/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PhaseOpTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def testPhaseImpl() {

        setup:
        def result = null
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()

        when:
        def map = [ : ]
        result = PhaseOp.phaseImpl(map, 2, 0, 'a', { it })
        then:
        result == null
        map == [ a:[0: ['a']] ]

        when:
        map = [ : ]
        result = PhaseOp.phaseImpl(map, 2, 0, 'a', { it })
        result = PhaseOp.phaseImpl(map, 2, 1, 'a', { it })
        then:
        result == ['a','a']
        map == [ a:[:] ]


        when:
        def r1
        def r2
        def r3
        map = [ : ]
        r1 = PhaseOp.phaseImpl(map, 3, 0, 'a', { it })
        r1 = PhaseOp.phaseImpl(map, 3, 1, 'a', { it })
        r1 = PhaseOp.phaseImpl(map, 3, 2, 'a', { it })

        r2 = PhaseOp.phaseImpl(map, 3, 0, 'b', { it })
        r2 = PhaseOp.phaseImpl(map, 3, 1, 'b', { it })
        r2 = PhaseOp.phaseImpl(map, 3, 2, 'b', { it })

        r3 = PhaseOp.phaseImpl(map, 3, 0, 'z', { it })
        r3 = PhaseOp.phaseImpl(map, 3, 1, 'z', { it })
        r3 = PhaseOp.phaseImpl(map, 3, 1, 'z', { it })
        r3 = PhaseOp.phaseImpl(map, 3, 2, 'z', { it })

        then:
        r1 == ['a','a','a']
        r2 == ['b','b','b']
        r3 == ['z','z','z']
        map == [ a:[:], b:[:], z:[ 1:['z']] ]

    }

}
