/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package nextflow.cache.sqlite

import java.nio.file.Path
import java.sql.Connection
import java.sql.DriverManager
import java.sql.PreparedStatement
import java.sql.ResultSet
import java.sql.SQLException

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.cache.CacheStore
import nextflow.exception.AbortOperationException
import nextflow.util.CacheHelper

/**
 * SQLite-based implementation of the Nextflow cache store
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SQLiteCacheStore implements CacheStore {

    /** The underlying SQLite connection */
    private Connection connection

    /** The session UUID */
    private UUID uniqueId

    /** The unique run name associated to this cache instance */
    private String runName

    /** The base folder against which the cache is located. Default: current working directory  */
    private Path baseDir

    /** The actual path where DB file is located */
    private Path dataDir

    /** The path to the database file */
    private Path dbFile

    /** Database connection synchronization object */
    private final Object connectionLock = new Object()

    SQLiteCacheStore(UUID uniqueId, String runName, Path home=null) {
        this.uniqueId = uniqueId
        this.runName = runName
        this.baseDir = home ?: Const.appCacheDir.toAbsolutePath()
        this.dataDir = baseDir.resolve("cache/$uniqueId")
        this.dbFile = dataDir.resolve("cache.db")
    }

    private void openDb() {
        // make sure the db path exists
        dataDir.mkdirs()
        
        try {
            // Initialize SQLite JDBC driver
            Class.forName("org.sqlite.JDBC")
            
            // Open connection to SQLite database with better configuration
            def url = "jdbc:sqlite:${dbFile.toString()}"
            connection = DriverManager.getConnection(url)
            
            // Disable autocommit for better control
            connection.setAutoCommit(true)
            
            // Configure SQLite for better performance and concurrency
            def stmt = connection.createStatement()
            stmt.execute("""
                PRAGMA journal_mode=WAL;
                PRAGMA synchronous=NORMAL;
                PRAGMA cache_size=10000;
                PRAGMA temp_store=MEMORY;
                PRAGMA busy_timeout=30000;
            """)
            stmt.close()
            
            // Create tables if they don't exist
            createTables()
            
        } catch (Exception e) {
            String msg
            if (e.message?.contains('database is locked') || e.message?.contains('SQLITE_BUSY')) {
                msg = "Unable to acquire lock on session with ID $uniqueId"
                msg += "\n\n"
                msg += "Common reasons for this error are:"
                msg += "\n - You are trying to resume the execution of an already running pipeline"
                msg += "\n - A previous execution was abruptly interrupted, leaving the session open"
                msg += '\n'
                msg += '\nYou can check for running processes that might be holding the database lock.'
                throw new SQLException(msg)
            } else {
                msg = "Can't open cache DB: $dbFile"
                msg += '\n\n'
                msg += "Nextflow needs to be executed in a shared file system that supports file locks.\n"
                msg += "Alternatively, you can run it in a local directory and specify the shared work\n"
                msg += "directory by using the `-w` command line option."
                throw new SQLException(msg, e)
            }
        }
    }

    private void createTables() {
        connection.createStatement().execute("""
            CREATE TABLE IF NOT EXISTS cache_entries (
                key BLOB PRIMARY KEY,
                value BLOB
            )
        """)
        
        connection.createStatement().execute("""
            CREATE TABLE IF NOT EXISTS cache_index (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_name TEXT NOT NULL,
                key BLOB NOT NULL,
                cached BOOLEAN NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        connection.createStatement().execute("""
            CREATE INDEX IF NOT EXISTS idx_cache_index_run_name ON cache_index(run_name)
        """)
    }


    @Override
    SQLiteCacheStore open() {
        openDb()
        // Clear previous index entries for this run
        connection.createStatement().execute("DELETE FROM cache_index WHERE run_name = '${runName}'")
        return this
    }

    @Override
    synchronized SQLiteCacheStore openForRead() {
        openDb()
        // Check if index exists for this run
        PreparedStatement stmt = null
        ResultSet rs = null
        try {
            stmt = connection.prepareStatement("SELECT COUNT(*) FROM cache_index WHERE run_name = ?")
            stmt.setString(1, runName)
            rs = stmt.executeQuery()
            rs.next()
            if (rs.getInt(1) == 0) {
                throw new AbortOperationException("Missing cache index for run: $runName")
            }
            return this
        } finally {
            try { rs?.close() } catch (SQLException e) { /* ignore */ }
            try { stmt?.close() } catch (SQLException e) { /* ignore */ }
        }
    }

    @Override
    void drop() {
        // Close connection first if it exists
        close()
        // Delete the entire data directory
        dataDir.deleteDir()
    }

    @Override
    synchronized void close() {
        try {
            connection?.close()
        } catch (SQLException e) {
            log.warn("Error closing SQLite cache store connection", e)
        }
    }

    @Override
    synchronized void writeIndex(HashCode key, boolean cached) {
        PreparedStatement stmt = null
        try {
            stmt = connection.prepareStatement("INSERT INTO cache_index (run_name, key, cached) VALUES (?, ?, ?)")
            stmt.setString(1, runName)
            stmt.setBytes(2, key.asBytes())
            stmt.setBoolean(3, cached)
            stmt.executeUpdate()
        } catch (SQLException e) {
            log.error("Error writing cache index", e)
            throw new RuntimeException("Failed to write cache index", e)
        } finally {
            try { stmt?.close() } catch (SQLException e) { /* ignore */ }
        }
    }

    @Override
    synchronized void deleteIndex() {
        PreparedStatement stmt = null
        try {
            stmt = connection.prepareStatement("DELETE FROM cache_index WHERE run_name = ?")
            stmt.setString(1, runName)
            stmt.executeUpdate()
        } catch (SQLException e) {
            log.error("Error deleting cache index", e)
            throw new RuntimeException("Failed to delete cache index", e)
        } finally {
            try { stmt?.close() } catch (SQLException e) { /* ignore */ }
        }
    }

    @Override
    Iterator<Index> iterateIndex() {
        try {
            // Create a fresh statement for this iteration to avoid conflicts
            def stmt = connection.prepareStatement("SELECT key, cached FROM cache_index WHERE run_name = ? ORDER BY id")
            stmt.setString(1, runName)
            final ResultSet rs = stmt.executeQuery()
            
            return new Iterator<Index>() {
                private Index nextItem = null
                private boolean fetched = false

                @Override
                boolean hasNext() {
                    if (!fetched) {
                        nextItem = fetch()
                        fetched = true
                    }
                    return nextItem != null
                }

                @Override
                Index next() {
                    if (!fetched) {
                        nextItem = fetch()
                        fetched = true
                    }
                    final result = nextItem
                    nextItem = fetch()
                    fetched = true
                    return result
                }

                private Index fetch() {
                    try {
                        if (rs.next()) {
                            final key = HashCode.fromBytes(rs.getBytes("key"))
                            final cached = rs.getBoolean("cached")
                            return new Index(key, cached)
                        } else {
                            try {
                                rs.close()
                                stmt.close()
                            } catch (SQLException e) {
                                // ignore cleanup errors
                            }
                            return null
                        }
                    } catch (SQLException e) {
                        log.error("Error reading cache index", e)
                        try {
                            rs.close()
                            stmt.close()
                        } catch (SQLException cleanupE) {
                            // ignore cleanup errors
                        }
                        throw new RuntimeException("Failed to read cache index", e)
                    }
                }
            }
        } catch (SQLException e) {
            log.error("Error iterating cache index", e)
            throw new RuntimeException("Failed to iterate cache index", e)
        }
    }

    @Override
    synchronized byte[] getEntry(HashCode key) {
        PreparedStatement stmt = null
        ResultSet rs = null
        try {
            stmt = connection.prepareStatement("SELECT value FROM cache_entries WHERE key = ?")
            stmt.setBytes(1, key.asBytes())
            rs = stmt.executeQuery()
            if (rs.next()) {
                return rs.getBytes("value")
            }
            return null
        } catch (SQLException e) {
            log.error("Error getting cache entry", e)
            throw new RuntimeException("Failed to get cache entry", e)
        } finally {
            try { rs?.close() } catch (SQLException e) { /* ignore */ }
            try { stmt?.close() } catch (SQLException e) { /* ignore */ }
        }
    }

    @Override
    synchronized void putEntry(HashCode key, byte[] value) {
        PreparedStatement stmt = null
        try {
            stmt = connection.prepareStatement("INSERT OR REPLACE INTO cache_entries (key, value) VALUES (?, ?)")
            stmt.setBytes(1, key.asBytes())
            stmt.setBytes(2, value)
            stmt.executeUpdate()
        } catch (SQLException e) {
            log.error("Error putting cache entry", e)
            throw new RuntimeException("Failed to put cache entry", e)
        } finally {
            try { stmt?.close() } catch (SQLException e) { /* ignore */ }
        }
    }

    @Override
    synchronized void deleteEntry(HashCode key) {
        PreparedStatement stmt = null
        try {
            stmt = connection.prepareStatement("DELETE FROM cache_entries WHERE key = ?")
            stmt.setBytes(1, key.asBytes())
            stmt.executeUpdate()
        } catch (SQLException e) {
            log.error("Error deleting cache entry", e)
            throw new RuntimeException("Failed to delete cache entry", e)
        } finally {
            try { stmt?.close() } catch (SQLException e) { /* ignore */ }
        }
    }
}