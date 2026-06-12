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

import com.google.common.hash.HashCode
import groovy.transform.TupleConstructor

/**
 * Defines the contract for a pluggable cache storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface CacheStore {

    @TupleConstructor
    static class Index {
        final HashCode key
        final boolean cached
    }

    CacheStore open()
    CacheStore openForRead()
    void close()
    void drop()

    byte[] getEntry(HashCode key)
    void putEntry(HashCode key, byte[] value)
    void deleteEntry(HashCode key)

    /**
     * Retrieve the final hash of a successful execution for the given content hash,
     * or {@code null} if no successful-hash index entry exists.
     */
    HashCode getSuccessfulHash(HashCode contentHash)

    /**
     * Map a task content hash to the final hash of a successful execution.
     */
    void putSuccessfulHash(HashCode contentHash, HashCode finalHash)

    /**
     * Remove the successful-hash index entry for the given content hash.
     */
    void deleteSuccessfulHash(HashCode contentHash)

    void writeIndex(HashCode key, boolean cached)
    Iterator<Index> iterateIndex()
    void deleteIndex()

}
