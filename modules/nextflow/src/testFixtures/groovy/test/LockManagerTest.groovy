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

package test


import nextflow.util.LockManager
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LockManagerTest extends Specification {

    def 'should reuse the same instance from the pool' () {
        given:
        def manager = new LockManager()

        when:
        def lock1 = manager.acquire(1)
        and:
        lock1.release()

        and:
        def copy = manager.acquire(1)
        then:
        copy.is(lock1)

    }

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
}
