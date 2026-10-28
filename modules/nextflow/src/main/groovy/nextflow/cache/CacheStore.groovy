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
     * Update an entry that already exists in this store, i.e. the read-modify-write of the
     * reference count / last-used stamp in {@link CacheDB#incTaskEntry}, whose value is a record
     * just read back through {@link #getEntry}.
     *
     * This is deliberately distinct from {@link #putEntry}, which stores a <b>new</b> entry: a
     * composite store has to send an update to the member that served the read — and may have to
     * drop it when that member is read-only — while a new entry must always go to the writable
     * one. Defaults to {@link #putEntry}, which is the correct behaviour for a single store.
     */
    default void updateEntry(HashCode key, byte[] value) { putEntry(key, value) }

    void writeIndex(HashCode key, boolean cached)
    Iterator<Index> iterateIndex()
    void deleteIndex()

}
