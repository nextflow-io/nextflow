/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cache

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx
import nextflow.util.CacheHelper
/**
 * Implements the cloud cache store
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class CloudCacheStore implements CacheStore {

    private final String LOCK_NAME = 'LOCK'

    private final int KEY_SIZE

    /** The session UUID */
    private UUID uniqueId

    /** The unique run name associated with this cache instance */
    private String runName

    /** The base path for the entire cache */
    private Path basePath

    /** The base path for this cache instance */
    private Path dataPath

    /** The lock file for this cache instance */
    private Path lock

    /** The path to the index file */
    private Path indexPath

    /** Index file input stream */
    private InputStream indexReader

    /** Index file output stream */
    private OutputStream indexWriter

    CloudCacheStore(UUID uniqueId, String runName, Path basePath=null) {
        this.KEY_SIZE = CacheHelper.hasher('x').hash().asBytes().size()
        this.uniqueId = uniqueId
        this.runName = runName
        this.basePath = basePath ?: defaultBasePath()
        this.dataPath = this.basePath.resolve("$uniqueId")
        this.lock = dataPath.resolve(LOCK_NAME)
        this.indexPath = dataPath.resolve("index.$runName")
    }

    private Path defaultBasePath() {
        final basePath = SysEnv.get('NXF_CLOUDCACHE_PATH')
        if( !basePath )
            throw new IllegalArgumentException("NXF_CLOUDCACHE_PATH must be defined when using the cloud cache store")

        return basePath as Path
    }

    @Override
    CloudCacheStore open() {
        acquireLock()
        indexWriter = new BufferedOutputStream(Files.newOutputStream(indexPath))
        return this
    }

    @Override
    CloudCacheStore openForRead() {
        if( !dataPath.exists() )
            throw new AbortOperationException("Missing cache directory: $dataPath")
        acquireLock()
        indexReader = Files.newInputStream(indexPath)
        return this
    }

    private void acquireLock() {
        if( lock.exists() ) {
            final msg = """
                Unable to acquire lock for session with ID ${uniqueId}

                Common reasons for this error are:
                - You are trying to resume the execution of an already running pipeline
                - A previous execution was abruptly interrupted, leaving the session open

                You can see the name of the conflicting run by inspecting the contents of the following path: ${lock}
                """
            throw new IOException(msg)
        }

        lock.text = runName
    }

    @Override
    void drop() {
        dataPath.deleteDir()
    }

    @Override
    void close() {
        FilesEx.closeQuietly(indexWriter)
        lock.delete()
    }

    @Override
    void writeIndex(HashCode key, boolean cached) {
        indexWriter.write(key.asBytes())
        indexWriter.write(cached ? 1 : 0)
    }

    @Override
    void deleteIndex() {
        indexPath.delete()
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
                if( indexReader.read(key) == -1 )
                    return null
                final cached = indexReader.read() == 1
                return new Index(HashCode.fromBytes(key), cached)
            }
        }
    }

    @Override
    byte[] getEntry(HashCode key) {
        try {
            return getCachePath(key).bytes
        }
        catch( NoSuchFileException e ) {
            return null
        }
    }

    @Override
    void putEntry(HashCode key, byte[] value) {
        getCachePath(key).bytes = value
    }

    @Override
    void deleteEntry(HashCode key) {
        getCachePath(key).delete()
    }

    private Path getCachePath(HashCode key) {
        dataPath.resolve(key.toString())
    }
}
