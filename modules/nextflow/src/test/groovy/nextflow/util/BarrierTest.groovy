/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import spock.lang.IgnoreIf
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BarrierTest extends Specification {

    def 'test is terminated' () {

        when:
        def barrier = new Barrier()
        then:
        barrier.isCompleted()

        when:
        barrier.register(1)
        then:
        !barrier.isCompleted()

        when:
        barrier.register(1)
        barrier.register(2)
        barrier.arrive(1)
        then:
        !barrier.isCompleted()

        when:
        barrier.arrive(2)
        then:
        barrier.isCompleted()

    }

    @IgnoreIf({ javaVersion == 1.7 })
    def 'test await termination' () {

        given:
        def barrier = new Barrier()
        def begin = System.currentTimeMillis()

        barrier.register(1)
        barrier.register(2)
        Thread.start { sleep 50; barrier.arrive(2) }
        Thread.start { sleep 100; barrier.arrive(1) }

        when:
        barrier.awaitCompletion()
        def elapsed = System.currentTimeMillis() - begin

        then:
        elapsed>=100 && elapsed<500
    }
}
