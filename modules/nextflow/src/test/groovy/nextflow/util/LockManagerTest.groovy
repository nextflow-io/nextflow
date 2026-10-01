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

package test

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger

import nextflow.util.LockManager
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LockManagerTest extends Specification {

    def 'should lock on the same key' () {
        given:
        def manager = new LockManager()

        when:
        int counter=0
        def lock1 = manager.acquire(1)
        Thread.start { def l = manager.acquire(1); counter++; l.release() }
        and:
        sleep 100
        then:
        counter==0

        when:
        lock1.release()
        sleep 100
        then:
        counter ==1
    }

    def 'should not lock on different keys' () {
        given:
        def manager = new LockManager()

        when:
        int counter=0
        def lock1 = manager.acquire(1)
        Thread.start { def l = manager.acquire(2); counter++; l.release() }
        and:
        sleep 100
        then:
        counter==1

    }

    def 'should not retain a key once it is released' () {
        given:
        def manager = new LockManager()

        when:
        def lock1 = manager.acquire('a')
        def lock2 = manager.acquire('b')
        then:
        manager.@entries.size() == 2

        when:
        lock1.release()
        then:
        manager.@entries.size() == 1

        when:
        lock2.release()
        then:
        manager.@entries.isEmpty()
    }

    def 'should keep the key until a reentrant lock is fully released' () {
        given:
        def manager = new LockManager()

        when:
        def lock1 = manager.acquire('a')
        def lock2 = manager.acquire('a')
        then:
        lock1.is(lock2)

        when:
        lock2.release()
        then:
        manager.@entries.size() == 1

        when:
        lock1.release()
        then:
        manager.@entries.isEmpty()
    }

    @Timeout(30)
    def 'should provide mutual exclusion and not retain keys under contention' () {
        given:
        def manager = new LockManager()
        def inside = (0..<4).collect { new AtomicInteger() }
        def violations = new AtomicInteger()
        def errors = new ConcurrentLinkedQueue<Throwable>()

        when:
        def threads = (0..<8).collect {
            Thread.start {
                try {
                    for( int i=0; i<5_000; i++ ) {
                        def key = i % 4
                        def lock = manager.acquire(key)
                        try {
                            if( inside[key].incrementAndGet() != 1 )
                                violations.incrementAndGet()
                            inside[key].decrementAndGet()
                        }
                        finally {
                            lock.release()
                        }
                    }
                }
                catch( Throwable t ) { errors.add(t) }
            }
        }
        threads*.join()

        then:
        errors.isEmpty()
        violations.get() == 0
        manager.@entries.isEmpty()
    }
}
