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

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.Executors

import spock.lang.Specification

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class TrackingSemaphoreTest extends Specification {

    def 'should track total and available permits'() {
        given:
        def sem = new TrackingSemaphore(['0', '1', '2'])

        expect:
        sem.totalPermits() == 3
        sem.availablePermits() == 3
    }

    def 'should hand out distinct ids in insertion order'() {
        given:
        def sem = new TrackingSemaphore(['0', '1', '2', '3'])

        when:
        def a = sem.acquire(2)
        then:
        a == ['0', '1']
        sem.availablePermits() == 2

        when:
        def b = sem.acquire(2)
        then:
        b == ['2', '3']
        sem.availablePermits() == 0
    }

    def 'should reuse ids after release'() {
        given:
        def sem = new TrackingSemaphore(['0', '1', '2'])

        when:
        def a = sem.acquire(3)
        then:
        a == ['0', '1', '2']
        sem.availablePermits() == 0

        when:
        sem.release(['1'])
        def b = sem.acquire(1)
        then:
        b == ['1']
        sem.availablePermits() == 0
    }

    def 'should never assign the same id twice under concurrent acquire/release'() {
        given:
        def sem = new TrackingSemaphore(['0', '1', '2', '3'])
        def inUse = Collections.synchronizedSet(new HashSet<String>())
        def errors = new ConcurrentLinkedQueue<String>()

        when:
        def pool = Executors.newFixedThreadPool(4)
        def futures = (1..2000).collect {
            pool.submit({
                final ids = sem.acquire(1)
                for( final id : ids ) {
                    // each acquired id must be exclusively held while in use
                    if( !inUse.add(id) )
                        errors.add("id $id assigned twice".toString())
                    if( ids.size() != 1 )
                        errors.add("expected 1 id, got ${ids.size()}".toString())
                }
                for( final id : ids )
                    inUse.remove(id)
                sem.release(ids)
            } as Runnable)
        }
        futures*.get()
        pool.shutdown()

        then:
        errors.isEmpty()
        sem.availablePermits() == 4
        sem.totalPermits() == 4
    }

}
