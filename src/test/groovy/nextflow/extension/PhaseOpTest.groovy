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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
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

    def testPhase() {

        setup:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        when:
        def result = ch1.phase(ch2)
        then:
        result.val == [1,1]
        result.val == [2,2]
        result.val == [3,3]

        result.val == Channel.STOP


        when:
        ch1 = Channel.from( [sequence: 'aaaaaa', key: 1], [sequence: 'bbbbbb', key: 2] )
        ch2 = Channel.from( [val: 'zzzz', id: 3], [val: 'xxxxx', id: 1], [val: 'yyyyy', id: 2])
        result = ch1.phase(ch2) { Map it ->
            if( it.containsKey('key') ) {
                return it.key
            }
            else if( it.containsKey('id') ) {
                return it.id
            }
            return null
        }
        then:

        result.val == [ [sequence: 'aaaaaa', key: 1], [val: 'xxxxx', id: 1] ]
        result.val == [ [sequence: 'bbbbbb', key: 2], [val: 'yyyyy', id: 2] ]
        result.val == Channel.STOP

    }

    def testPhaseWithRemainder() {

        def ch1
        def ch2
        def result

        when:
        ch1 = Channel.from( 1,2,3 )
        ch2 = Channel.from( 1,0,0,2,7,8,9,3 )
        result = ch1.phase(ch2, remainder: true)

        then:
        result.val == [1,1]
        result.val == [2,2]
        result.val == [3,3]
        result.val == [null,0]
        result.val == [null,0]
        result.val == [null,7]
        result.val == [null,8]
        result.val == [null,9]
        result.val == Channel.STOP


        when:
        ch1 = Channel.from( 1,0,0,2,7,8,9,3 )
        ch2 = Channel.from( 1,2,3 )
        result = ch1.phase(ch2, remainder: true)

        then:
        result.val == [1,1]
        result.val == [2,2]
        result.val == [3,3]
        result.val == [0,null]
        result.val == [0,null]
        result.val == [7,null]
        result.val == [8,null]
        result.val == [9,null]
        result.val == Channel.STOP
    }

    def 'should phase entries' () {
        given:
        def ch1 = Channel.from(['sample1', 1], ['sample2', 2], ['sample3', 3])
        def ch2 = Channel.from(['sample1', 4], ['sample3', 6], ['sample2', 5])

        when:
        def result = ch1.phase(ch2).toList().getVal()
        then:
        result.size() == 3
        result.contains( [['sample1', 1], ['sample1', 4]] )
        result.contains( [['sample2', 2], ['sample2', 5]])
        result.contains( [['sample3', 3], ['sample3', 6]] )
    }


}
