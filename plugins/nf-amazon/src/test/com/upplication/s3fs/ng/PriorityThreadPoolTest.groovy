/*
 * Copyright 2020-2022, Seqera Labs
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

package com.upplication.s3fs.ng

import com.upplication.s3fs.ng.PriorityThreadPool
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PriorityThreadPoolTest extends Specification {

    def task(int priority, Closure action) {
        new PriorityThreadPool.PriorityCallable(priority) {
            @Override
            Object call() throws Exception {
                return action.call()
            }
        }
    }
    def 'should order tasks executions' () {
        given:
        def pool = PriorityThreadPool.create('foo', 1, 100)
        def sequence = Collections.synchronizedList([])

        when:
        def f1 = pool.submit( task(0,  { sleep 100; sequence.add("A") }) )
        def f4 = pool.submit( task(50,  { sleep 100;  sequence.add("D") } ) )
        def f3 = pool.submit( task(30,  { sleep 100;  sequence.add("C") } ) )
        def f2 = pool.submit( task(20,  { sleep 100;  sequence.add("B") }) )
        and:
        f1.get()
        f2.get()
        f3.get()
        f4.get()

        then:
        sequence == ['A','B','C', 'D']
    }

}
