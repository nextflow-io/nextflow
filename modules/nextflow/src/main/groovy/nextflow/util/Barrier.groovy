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

import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic

/**
 * An easy implementation of a barrier that wait for the "arrival" of multiple
 * registered parties.
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Barrier {

    private IdentityHashMap<Object,Boolean> parties = new IdentityHashMap()

    /**
     * Define a critical section
     */
    private final Lock section

    private volatile boolean terminated

    /**
     * Signal condition
     */
    private final Condition condition

    Barrier() {
        section = new ReentrantLock()
        condition = section.newCondition()
    }

    /**
     * Register a "party" that is needed to wait the arrival
     *
     * @param item An object instance
     */
    void register( item ) {
        section.withLock {
            parties.put(item,false)
        }
    }

    /**
     * Signal the arrival of a "party"
     *
     * @see #register(java.lang.Object)
     *
     * @param item A object instance previously registered
     */
    void arrive( item ) {
        section.withLock {
            if( !parties.containsKey(item))
                throw new IllegalStateException("Barrier invalid state -- not a registered item: $item")

            parties.put(item,true)
            condition.signal()
        }
    }

    /**
     * Block until all registered parties "arrive" to the barrier condition
     */
    void awaitCompletion() {
        if( !parties.size() || terminated )
            return

        section.withLock {
            while( true ) {
                if( allArrived() || terminated ) break
                condition.await()
            }
        }
    }

    /**
     * @return {@code true} when all parties have arrived
     */
    boolean isCompleted() {
        section.withLock { allArrived() }
    }

    private boolean allArrived() {
        boolean result = true
        for( boolean complete : parties.values() )
            result = result && complete

        return result
    }

    /**
     * Force the barrier termination
     */
    void forceTermination() {
        terminated = true
        section.withLock {
            for( def key: parties.keySet() ) {
                parties.put(key, true)
            }
            condition.signal()
        }
    }

    int size() { parties.size() }
}
