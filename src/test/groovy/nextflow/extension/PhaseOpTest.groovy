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
