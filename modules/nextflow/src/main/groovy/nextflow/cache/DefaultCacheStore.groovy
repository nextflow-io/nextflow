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
 * into the local directory using an embedded LevelDVB instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DefaultCacheStore implements CacheStore {

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
}
