/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import java.util.function.Function

import groovy.transform.CompileStatic
/**
 * A lock manager that allows the acquire of a lock on a unique key object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LockManager {

    private int MAX_SIZE = 200

    /**
     * Maintain a pool of lock handler to reduce garbage collection
     */
    private List<LockHandle> pool = new ArrayList<>(MAX_SIZE)

    /**
     * Associate a lock handle for each key
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
        LockHandle result = entries.computeIfAbsent(key,newLock())
        result.sync.lock()
        result.count++
        return result
    }

    private Function<Object, LockHandle> newLock() {
        new Function<Object, LockHandle>() {
            @Override
            LockHandle apply(Object key) {
                return getOrCreate0(key)
            }
        }
    }


    private synchronized LockHandle getOrCreate0(key) {
        if( pool.size() ) {
            def handle = pool.remove(pool.size()-1)
            handle.key = key
            return handle
        }
        new LockHandle(key)
    }

    private synchronized void release0(LockHandle handle) {
        entries.remove(handle.key)
        handle.key = null
        if( pool.size()<MAX_SIZE )
            pool.add(handle)
    }


    class LockHandle {
        Lock sync
        volatile Object key
        volatile int count

        LockHandle(key) {
            this.key = key
            this.sync = new ReentrantLock()
        }

        void release() {
            if( count-- == 0 ) {
                release0(this)
            }
            sync.unlock()
        }
    }
}
