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

package nextflow.extension

import nextflow.Channel
import spock.lang.Specification

import static test.ScriptHelper.runDataflow
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
        def result = runDataflow {
            tuples.channel().groupTuple(size: 2)
        }
        then:
        result.val == [k1, ['a', 'b'] ]
        result.val == [k1, ['d', 'c'] ]
        result.val == [k2, ['x', 'y'] ]
        result.val == Channel.STOP

        when:
        // here the size is inferred by the key itself
        result = runDataflow {
            tuples.channel().groupTuple()
        }
        then:
        result.val == [k1, ['a', 'b'] ]
        result.val == [k1, ['d', 'c'] ]
        result.val == [k2, ['x', 'y', 'z'] ]
        result.val == Channel.STOP


        when:
        result = runDataflow {
            tuples.channel().groupTuple(remainder: true)
        }
        then:
        result.val == [k1, ['a', 'b'] ]
        result.val == [k1, ['d', 'c'] ]
        result.val == [k2, ['x', 'y', 'z'] ]
        result.val == [k3, ['q']]
        result.val == [k1, ['f']]
        result.val == Channel.STOP
    }

    def 'should handle GString vs String keys correctly' () {
        given:
        def sampleName = "sample1"
        def gstringKey = "${sampleName}"  // GStringImpl
        def stringKey = "sample1"         // String

        def tuples = [
                [gstringKey, 'file1.txt'],
                [stringKey, 'file2.txt'],
                [gstringKey, 'file3.txt']
        ]

        when:
        def result = runDataflow {
            tuples.channel().groupTuple()
        }
        then:
        // Without fix, this would create separate groups for GString and String keys
        // With fix, they should be grouped together
        result.val == [gstringKey, ['file1.txt', 'file2.txt', 'file3.txt']]
        result.val == Channel.STOP
    }

    def 'should preserve non-string key types' () {
        given:
        def tuples = [
                [1, 'file1.txt'],           // Integer key
                [new File('/tmp'), 'file2.txt'],  // File key
                [[a: 1], 'file3.txt']       // Map key
        ]

        when:
        def result = runDataflow {
            tuples.channel().groupTuple()
        }
        then:
        result.val == [1, ['file1.txt']]
        result.val == [new File('/tmp'), ['file2.txt']]
        result.val == [[a: 1], ['file3.txt']]
        result.val == Channel.STOP
    }
}
