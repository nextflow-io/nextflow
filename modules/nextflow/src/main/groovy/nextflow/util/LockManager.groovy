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

package nextflow.util

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import java.util.function.BiFunction

import groovy.transform.CompileStatic
/**
 * A lock manager that allows the acquire of a lock on a unique key object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LockManager {

    /**
     * Associate a lock handle for each key. An entry lives as long as
     * there is at least one thread holding or waiting for its lock
     */
    private ConcurrentHashMap<Object, LockHandle> entries = new ConcurrentHashMap<>()


    /**
     * Acquire a lock. The lock needs to be released using the `release` method eg.
     *
     * def lock = lockManager.acquire(key)
     * try {
     *     // safe code
     * }
     * finally {
     *     lock.release()
     * }
     *
     * @param key A key object over which the lock needs to be acquired
     * @return The lock handler
     */
    LockHandle acquire(key) {
        final handle = entries.compute(key, (k, h) -> {
            h = h ?: new LockHandle(k)
            h.count++
            return h
        } as BiFunction<Object, LockHandle, LockHandle>)
        handle.sync.lock()
        return handle
    }

    class LockHandle {
        final Lock sync = new ReentrantLock()
        final Object key
        // updated only inside the map `compute` functions
        int count

        LockHandle(key) {
            this.key = key
        }

        void release() {
            sync.unlock()
            entries.computeIfPresent(key, (k, h) -> --h.count == 0 ? null : h)
        }
    }
}
