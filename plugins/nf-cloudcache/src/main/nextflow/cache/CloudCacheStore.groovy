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

package nextflow.cache

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
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

    private final int KEY_SIZE

    /** The session UUID */
    private UUID uniqueId

    /** The unique run name associated with this cache instance */
    private String runName

    /** The base path for the entire cache */
    private Path basePath

    /** The base path for this cache instance */
    private Path dataPath

    /** The path to the index file */
    private Path indexPath

    /** Index file input stream */
    private InputStream indexReader

    /** Index file output stream */
    private OutputStream indexWriter

    CloudCacheStore(UUID uniqueId, String runName, Path basePath) {
        assert uniqueId, "Missing cloudcache 'uniqueId' argument"
        assert runName, "Missing cloudcache 'runName' argument"
        assert basePath, "Missing cloudcache 'basePath' argument"
        this.KEY_SIZE = CacheHelper.hasher('x').hash().asBytes().size()
        this.uniqueId = uniqueId
        this.runName = runName
        this.basePath = basePath
        this.dataPath = this.basePath.resolve("$uniqueId")
        this.indexPath = dataPath.resolve("index.$runName")
    }

    @Override
    CloudCacheStore open() {
        indexWriter = new BufferedOutputStream(Files.newOutputStream(indexPath))
        return this
    }

    @Override
    CloudCacheStore openForRead() {
        if( !dataPath.exists() )
            throw new AbortOperationException("Missing cache directory: $dataPath")
        indexReader = Files.newInputStream(indexPath)
        return this
    }

    @Override
    void drop() {
        dataPath.deleteDir()
    }

    @Override
    void close() {
        FilesEx.closeQuietly(indexWriter)
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
