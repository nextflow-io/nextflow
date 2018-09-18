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

import nextflow.Channel

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GroupTupleOpTest extends Specification {



    def 'should reuse the same key' () {

        given:
        def key1 = ['a', 'b', 'c']
        def key2 = ['a', 'b', 'c']
        def key3 = ['a', 'B', 'd']

        def map = [:]
        def result

        when:
        result = map.getOrCreate( key1 ) { return 1 }
        then:
        result == 1
        map.get(key1) == 1
        map.get(key2) == 1
        map.get(key3) == null
        map.getOrCreate(key2) { return 10 } == 1

        map.containsKey(key1)
        map.containsKey(key2)
        !map.containsKey(key3)

        when:
        result = map.getOrCreate( key3 ) { return 3 }
        then:
        result == 3
        map.getOrCreate(key3) { return 0 } == 3
    }

    def 'should reuse same key with GroupSize' () {

        given:
        def key1 = [new GroupKey('A',1) ]
        def key2 = [new GroupKey('A',1) ]
        def key3 = [new GroupKey('B',1) ]

        def map = [:]
        def result

        when:
        result = map.getOrCreate( key1 ) { return 1 }
        then:
        result == 1
        map.get(key1) == 1
        map.get(key2) == 1
        map.get(key3) == null
        map.getOrCreate(key2) { return 10 } == 1

    }

    def 'should fetch groupsize' () {

        given:
        def kep = GroupTupleOp.sizeBy(VALUES)

        expect:
        kep == EXPECT

        where:
        EXPECT  | VALUES
        10      | [new GroupKey('A',10)]
        20      | [new GroupKey('A',20)]
        0       | [new GroupKey('A',20), 'B']
        0       | ['A']
        0       | []
    }


    def 'should group items using dyn group size' () {
        given:
        def k1 = new GroupKey("k1", 2)
        def k2 = new GroupKey('k2', 3)
        def k3 = new GroupKey('k3', 4)

        def tuples = [
                [k1,'a'], [k1,'b'], [k2,'x'], [k3, 'q'], [k1,'d'], [k1,'c'], [k2, 'y'], [k1,'f'], [k2, 'z']
        ]

        when:
        // here the size is defined as operator argument
        def result = tuples.channel().groupTuple(size: 2)
        then:
        result.val == [k1, ['a', 'b'] ]
        result.val == [k1, ['d', 'c'] ]
        result.val == [k2, ['x', 'y'] ]
        result.val == Channel.STOP

        when:
        // here the size is inferred by the key itself
        result = tuples.channel().groupTuple()
        then:
        result.val == [k1, ['a', 'b'] ]
        result.val == [k1, ['d', 'c'] ]
        result.val == [k2, ['x', 'y', 'z'] ]
        result.val == Channel.STOP


        when:
        result = tuples.channel().groupTuple(remainder: true)
        then:
        result.val == [k1, ['a', 'b'] ]
        result.val == [k1, ['d', 'c'] ]
        result.val == [k2, ['x', 'y', 'z'] ]
        result.val == [k3, ['q']]
        result.val == [k1, ['f']]
        result.val == Channel.STOP
    }
}
