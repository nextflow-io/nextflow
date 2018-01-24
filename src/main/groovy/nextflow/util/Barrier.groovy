/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
