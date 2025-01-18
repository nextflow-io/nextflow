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

import nextflow.Channel
import nextflow.Session
import nextflow.extension.op.OpDatum
import nextflow.prov.OperatorRun
import nextflow.util.ArrayBag
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GroupTupleOpTest extends Specification {

    def setup() {
        new Session()
    }

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

    def 'should unwrap tuple values' () {
        given:
        def r1 = new OperatorRun(new HashSet<Integer>([1,2]))
        def r2 = new OperatorRun(new HashSet<Integer>([2,3]))
        def d1 = OpDatum.of('A', r1)
        def d2 = OpDatum.of('B', r1)
        def d3 = OpDatum.of('C', r2)
        def d4 = OpDatum.of('D', r2)
        and:
        def bag1 = new ArrayBag(d1, d2)
        def bag2 = new ArrayBag(d3, d4)
        def key = Mock(GroupKey)
        and:
        def tuple = [key, bag1, bag2]

        when:
        def result = GroupTupleOp.unwrapValues(tuple)
        then:
        tuple[0] == key
        and:
        tuple[1] == new ArrayBag<>('A','B')
        and:
        tuple[2] == new ArrayBag<>('C', 'D')
        and:
        result.inputIds == [1,2,3] as Set
        
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
        result.unwrap() == [k1, ['a', 'b'] ]
        result.unwrap() == [k1, ['d', 'c'] ]
        result.unwrap() == [k2, ['x', 'y'] ]
        result.unwrap() == Channel.STOP

        when:
        // here the size is inferred by the key itself
        result = tuples.channel().groupTuple()
        then:
        result.unwrap() == [k1, ['a', 'b'] ]
        result.unwrap() == [k1, ['d', 'c'] ]
        result.unwrap() == [k2, ['x', 'y', 'z'] ]
        result.unwrap() == Channel.STOP


        when:
        result = tuples.channel().groupTuple(remainder: true)
        then:
        result.unwrap() == [k1, ['a', 'b'] ]
        result.unwrap() == [k1, ['d', 'c'] ]
        result.unwrap() == [k2, ['x', 'y', 'z'] ]
        result.unwrap() == [k3, ['q']]
        result.unwrap() == [k1, ['f']]
        result.unwrap() == Channel.STOP
    }
}
