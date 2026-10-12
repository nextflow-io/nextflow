/*
 * Copyright 2013-2026, Seqera Labs
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
 */

package nextflow.cache

import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.util.CacheHelper
import org.iq80.leveldb.DB
import org.iq80.leveldb.Options
import org.iq80.leveldb.impl.Iq80DBFactory
/**
 * Implement the default nextflow cache store that save the cache data
 * into the local directory using an embedded LevelDVB instance *
 *
 * KEY-SPACE CONTRACT
 * ------------------
 * This single LevelDB instance holds two disjoint key namespaces:
 *
 *   - task entries:        key = <hash>            (exactly KEY_SIZE bytes)
 *   - realized-hash idx:  key = 0x01 + <hash>     (KEY_SIZE + 1 bytes)
 *
 * The {@link #REALIZED_HASH_INDEX_PREFIX} byte makes index keys one byte longer than
 * entry keys, so the two namespaces can never collide (LevelDB compares the
 * full byte sequence). All key construction goes through {@link #realizedIndexKey}
 * for index keys and {@code key.asBytes()} for entry keys -- do not build raw
 * keys elsewhere.
 *
 * NOTE: the store is never iterated by LevelDB key (record enumeration uses
 * the separate flat index file via {@link #iterateIndex}). If a key scan is
 * ever added, it MUST filter by this contract (length / prefix) so it does
 * not treat index pointers as task-entry records.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DefaultCacheStore implements CacheStore {

    /** Prefix byte for realized-hash index keys -- see KEY-SPACE CONTRACT above. */
    private static final byte REALIZED_HASH_INDEX_PREFIX = 0x01

    /** The underlying Level DB instance */
    private DB db

    /** The session UUID */
    private UUID uniqueId

    /** The unique run name associated to this cache instance */
    private String runName

    /** The base folder against which the cache is located. Default: current working directory  */
    private Path baseDir

    /** The actual path where DB and index file are located */
    private Path dataDir

    private final int KEY_SIZE

    /** The path to the index file */
    private Path indexFile

    /** Index file read/write handle */
    private RandomAccessFile indexHandle


    DefaultCacheStore(UUID uniqueId, String runName, Path home=null) {
        this.KEY_SIZE = CacheHelper.hasher('x').hash().asBytes().size()
        this.uniqueId = uniqueId
        this.runName = runName
        this.baseDir = home ?: Const.appCacheDir.toAbsolutePath()
        this.dataDir = baseDir.resolve("cache/$uniqueId")
        this.indexFile = dataDir.resolve("index.$runName")
    }

    private void openDb() {
        // make sure the db path exists
        dataDir.mkdirs()
        // open a LevelDB instance
        final file=dataDir.resolve('db').toFile()
        try {
            db = Iq80DBFactory.@factory.open(file, new Options().createIfMissing(true))
        }
        catch( Exception e ) {
            if( e.message?.startsWith('Unable to acquire lock') ) {
                final msg = """\
                    Unable to acquire lock on session with ID $uniqueId

                    Common reasons for this error are:
                     - You are trying to resume the execution of an already running pipeline
                     - A previous execution was abruptly interrupted, leaving the session open

                    You can see which process is holding the lock file by using the following command:
                     - lsof $file/LOCK""".stripIndent()
                throw new IOException(msg)
            }
            else {
                final msg = """\
                    Can't open cache DB: $file

                    Cause: ${e.message}

                    This can happen when the cache DB is located on a file system that does not
                    support file locks, which is common for network file systems (NFS, Lustre, GPFS, …).
                    In that case you can either:
                      1. Set the NXF_CACHE_DIR environment variable to a lock-capable path.
                      2. Run Nextflow from a local directory and point the work directory to your
                         shared file system using the `-w` command line option.""".stripIndent()
                throw new IOException(msg, e)
            }
        }
    }

    @Override
    DefaultCacheStore open() {
        openDb()
        //
        indexFile.delete()
        indexHandle = new RandomAccessFile(indexFile.toFile(), 'rw')
        return this
    }

    @Override
    DefaultCacheStore openForRead() {
        openDb()
        if( !indexFile.exists() )
            throw new AbortOperationException("Missing cache index file: $indexFile")
        indexHandle = new RandomAccessFile(indexFile.toFile(), 'r')
        return this
    }

    @Override
    void drop() {
        dataDir.deleteDir()
    }

    @Override
    void close() {
        indexHandle.closeQuietly()
        db.closeQuietly()
    }

    @Override
    void writeIndex(HashCode key, boolean cached) {
        indexHandle.write(key.asBytes())
        indexHandle.writeBoolean(cached)
    }

    @Override
    void deleteIndex() {
        indexFile.delete()
    }

    @Override
    Iterator<Index> iterateIndex() {
        return new Iterator<Index>() {
            private Index next

            {
                next = fetch()
            }

            @Override
            boolean hasNext() {
                return next != null
            }

            @Override
            Index next() {
                final result = next
                next = fetch()
                return result
            }

            private Index fetch() {
                byte[] key = new byte[KEY_SIZE]
                if( indexHandle.read(key) == -1 )
                    return null
                final cached = indexHandle.readBoolean()
                return new Index (HashCode.fromBytes(key), cached)
            }
        }
    }

    @Override
    byte[] getEntry(HashCode key) {
        return db.get(key.asBytes())
    }

    @Override
    void putEntry(HashCode key, byte[] value) {
        db.put(key.asBytes(), value)
    }

    @Override
    void deleteEntry(HashCode key) {
        db.delete(key.asBytes())
    }

    private byte[] realizedIndexKey(HashCode contentHash) {
        final h = contentHash.asBytes()
        final k = new byte[h.length + 1]
        k[0] = REALIZED_HASH_INDEX_PREFIX
        System.arraycopy(h, 0, k, 1, h.length)
        return k
    }

    @Override
    HashCode getRealizedHash(HashCode contentHash) {
        final value = db.get(realizedIndexKey(contentHash))
        return value != null ? HashCode.fromBytes(value) : null
    }

    @Override
    void putRealizedHash(HashCode contentHash, HashCode finalHash) {
        db.put(realizedIndexKey(contentHash), finalHash.asBytes())
    }

    @Override
    void deleteRealizedHash(HashCode contentHash) {
        db.delete(realizedIndexKey(contentHash))
    }
}
