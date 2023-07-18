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

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import nextflow.cache.CacheStore
import nextflow.exception.AbortOperationException
import nextflow.util.CacheHelper
/**
 * Implements the path-based cache store
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PathCacheStore implements CacheStore {

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

    PathCacheStore(UUID uniqueId, String runName, Path basePath=null) {
        this.KEY_SIZE = CacheHelper.hasher('x').hash().asBytes().size()
        this.uniqueId = uniqueId
        this.runName = runName
        this.basePath = basePath ?: defaultBasePath()
        this.dataPath = this.basePath.resolve("$uniqueId")
        this.lock = dataPath.resolve(LOCK_NAME)
    }

    private Path defaultBasePath() {
        final basePath = System.getenv('NXF_CACHE_PATH')
        if( !basePath )
            throw new IllegalArgumentException("NXF_CACHE_PATH must be defined when using the path-based cache store")

        return basePath as Path
    }

    @Override
    PathCacheStore open() {
        acquireLock()
        return this
    }

    @Override
    PathCacheStore openForRead() {
        if( !dataPath.exists() )
            throw new AbortOperationException("Missing cache directory: $dataPath")
        acquireLock()
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
        lock.delete()
    }

    @Override
    void writeIndex(HashCode key, boolean cached) {
        getCachePath(key).bytes = [ cached as byte ] as byte[]
    }

    @Override
    void deleteIndex() {}

    @Override
    Iterator<Index> iterateIndex() {
        final entries = dataPath.list() as List<String>
        entries.remove(LOCK_NAME)

        return new Iterator<Index>() {
            private Iterator<String> delegate

            {
                delegate = entries.iterator()
            }

            @Override
            boolean hasNext() {
                return delegate.hasNext()
            }

            @Override
            Index next() {
                final key = HashCode.fromString(delegate.next())
                final cached = getCachePath(key).bytes[0] as Boolean
                return new Index(key, cached)
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
