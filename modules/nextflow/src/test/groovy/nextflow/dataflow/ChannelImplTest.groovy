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

package nextflow.dataflow

import nextflow.dataflow.ChannelNamespace as channel
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.util.HashBag
import spock.lang.Specification
import spock.lang.Timeout

import static nextflow.Nextflow.record
import static nextflow.Nextflow.tuple
import static test.ScriptHelper.runDataflow
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(10)
class ChannelImplTest extends Specification {

    def testCollect() {

        when:
        def result = runDataflow {
            channel.of(1,2,3).collect()
        }
        then:
        result.val == new HashBag<>([1,2,3])

        when:
        result = runDataflow {
            channel.of(1).collect()
        }
        then:
        result.val == new HashBag<>([1])

        when:
        result = runDataflow {
            channel.empty().collect()
        }
        then:
        result.val == new HashBag<>([])
    }

    def testCross() {

        when:
        def result = runDataflow {
            def ch1 = channel.of( [1, 'x'], [2,'y'], [3,'z'] )
            def ch2 = channel.of( [1,11], [1,13], [2,21], [2,22], [2,23], [4,1], [4,2] )
            ch1.cross(ch2)
                .collect()
                .flatMap { vals ->
                    vals.toSorted { v -> [v[0][0], v[1][0]] }
                }
        }

        then:
        result.val == tuple([1, 'x'], [1,11])
        result.val == tuple([1, 'x'], [1,13])
        result.val == tuple([1, 'x'], [2,21])
        result.val == tuple([1, 'x'], [2,22])
        result.val == tuple([1, 'x'], [2,23])
        result.val == tuple([1, 'x'], [4,1])
        result.val == tuple([1, 'x'], [4,2])

        result.val == tuple([2, 'y'], [1,11])
        result.val == tuple([2, 'y'], [1,13])
        result.val == tuple([2, 'y'], [2,21])
        result.val == tuple([2, 'y'], [2,22])
        result.val == tuple([2, 'y'], [2,23])
        result.val == tuple([2, 'y'], [4,1])
        result.val == tuple([2, 'y'], [4,2])

        result.val == tuple([3, 'z'], [1,11])
        result.val == tuple([3, 'z'], [1,13])
        result.val == tuple([3, 'z'], [2,21])
        result.val == tuple([3, 'z'], [2,22])
        result.val == tuple([3, 'z'], [2,23])
        result.val == tuple([3, 'z'], [4,1])
        result.val == tuple([3, 'z'], [4,2])

        result.val == CH.stop()
    }

    def testFilter() {

        when:
        def c1 = runDataflow {
            channel.of(1,2,3,4,5).filter { v -> v > 3 }
        }
        then:
        c1.val == 4
        c1.val == 5
        c1.val == CH.stop()

        when:
        def c2 = runDataflow {
            channel.of('hola','hello','cioa','miao').filter { v -> v =~ /^h.*/ }
        }
        then:
        c2.val == 'hola'
        c2.val == 'hello'
        c2.val == CH.stop()

        when:
        def c3 = runDataflow {
            channel.of('hola','hello','cioa','miao').filter { v -> v ==~ /^h.*/ }
        }
        then:
        c3.val == 'hola'
        c3.val == 'hello'
        c3.val == CH.stop()
    }

    def testFlatMap() {

        when:
        def result = runDataflow {
            channel.of(1,2,3).flatMap { v -> [v, v*2] }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 2
        result.val == 4
        result.val == 3
        result.val == 6
        result.val == CH.stop()
    }

    def testFlatMapWithTuples() {

        when:
        def result = runDataflow {
            channel.of( [1,2], ['a','b'] ).flatMap { vals -> [vals, vals.reverse()] }
        }
        then:
        result.val == [1,2]
        result.val == [2,1]
        result.val == ['a','b']
        result.val == ['b','a']
        result.val == CH.stop()
    }

    def testFlatMapDefault() {

        when:
        def result = runDataflow {
            channel.of( [1,2], ['a',['b','c']] ).flatMap()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 'a'
        result.val == ['b','c']  // <-- nested list are preserved
        result.val == CH.stop()
    }

    def testGroupBy() {

        when:
        def result = runDataflow {
            channel.of([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'c'], [2, 'y'], [3, 'q'])
                .groupBy()
        }

        then:
        result.val == tuple(1, new HashBag<>(['a', 'b', 'c']))
        result.val == tuple(2, new HashBag<>(['x', 'y']))
        result.val == tuple(3, new HashBag<>(['q', 'q']))
        result.val == CH.stop()

        when:
        result = runDataflow {
            channel.of([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'d'], [1,'c'], [2, 'y'], [1,'f'])
                .groupBy()
        }

        then:
        result.val == tuple(1, new HashBag<>(['a', 'b', 'd', 'c', 'f']))
        result.val == tuple(2, new HashBag<>(['x', 'y']))
        result.val == tuple(3, new HashBag<>(['q']))
        result.val == CH.stop()
    }

    def testJoin() {
        when:
        def result = runDataflow {
            def ch1 = channel.of(record(id: 'X', a: 1), record(id: 'Y', a: 2), record(id: 'Z', a: 3), record(id: 'P', a: 7))
            def ch2 = channel.of(record(id: 'Z', b: 6), record(id: 'Y', b: 5), record(id: 'X', b: 4))
            ch1.join(ch2, by: 'id').collect()
        }.getVal()
        then:
        result.size() == 3
        result.contains( record(id: 'X', a: 1, b: 4) )
        result.contains( record(id: 'Y', a: 2, b: 5) )
        result.contains( record(id: 'Z', a: 3, b: 6) )
    }

    def testJoinDuplicates() {
        when:
        def result = runDataflow {
            def ch1 = channel.of(record(id: 'X', a: 1), record(id: 'X', a: 3))
            def ch2 = channel.of(record(id: 'X', b: 2), record(id: 'X', b: 4))
            ch1.join(ch2, by: 'id').collect()
        }.getVal()
        then:
        result.size() == 4
        result.contains( record(id: 'X', a: 1, b: 2) )
        result.contains( record(id: 'X', a: 1, b: 4) )
        result.contains( record(id: 'X', a: 3, b: 2) )
        result.contains( record(id: 'X', a: 3, b: 4) )
    }

    def testJoinWithRemainder() {

        when:
        def result = runDataflow {
            def ch1 = channel.of(record(id: 'X', a: 1), record(id: 'Y', a: 2), record(id: 'Z', a: 3), record(id: 'P', a: 7))
            def ch2 = channel.of(record(id: 'Z', b: 6), record(id: 'Y', b: 5), record(id: 'X', b: 4), record(id: 'Q', b: 8))
            ch1.join(ch2, by: 'id', remainder: true).collect()
        }.getVal()
        then:
        result.size() == 5
        result.contains( record(id: 'X', a: 1, b: 4) )
        result.contains( record(id: 'Y', a: 2, b: 5) )
        result.contains( record(id: 'Z', a: 3, b: 6) )
        result.contains( record(id: 'P', a: 7) )
        result.contains( record(id: 'Q', b: 8) )

        when:
        result = runDataflow {
            def ch1 = channel.of(record(id: 'X', a: 1), record(id: 'Y', a: 2), record(id: 'Z', a: 3), record(id: 'P', a: 7))
            def ch2 = channel.empty()
            ch1.join(ch2, by: 'id', remainder: true).collect()
        }.getVal()
        then:
        result.size() == 4
        result.contains( record(id: 'X', a: 1) )
        result.contains( record(id: 'Y', a: 2) )
        result.contains( record(id: 'Z', a: 3) )
        result.contains( record(id: 'P', a: 7) )
    }

    def testJoinError() {

        when:
        def result = runDataflow {
            def ch1 = channel.of( 1,2,3 )
            def ch2 = channel.of( 1,0,0,2,7,8,9,3 )
            ch1.join(ch2).collect()
        }.getVal()
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains 'Operator `join` requires the `by` option'
    }

    def testMap() {
        when:
        def result = runDataflow {
            channel.of(1,2,3).map { "Hello $it" }
        }
        then:
        result.val == 'Hello 1'
        result.val == 'Hello 2'
        result.val == 'Hello 3'
        result.val == CH.stop()
    }

    def testMapParamExpanding() {

        when:
        def result = runDataflow {
            channel.of(1,2,3).map { [it, it] }.map { x, y -> x+y }
        }
        then:
        result.val == 2
        result.val == 4
        result.val == 6
        result.val == CH.stop()
    }

    def testMix() {
        when:
        def result = runDataflow {
            def c1 = channel.of( 1,2,3 )
            def c2 = channel.of( 'a','b' )
            def c3 = channel.value( 'z' )
            c1.mix(c2,c3).collect()
        }.val

        then:
        1 in result
        2 in result
        3 in result
        'a' in result
        'b' in result
        'z' in result
        !('c' in result)
    }

    def testReduce() {

        when:
        def result = runDataflow {
            channel.of(1,2,3,4,5).reduce { a, e -> a += e }
        }
        then:
        result.getVal() == 15

        when:
        result = runDataflow {
            channel.of(99).reduce { a, e -> a += e }
        }
        then:
        result.getVal() == 99

        when:
        result = runDataflow {
            channel.empty().reduce { a, e -> a += e }
        }
        then:
        result.getVal() == null
    }


    def testReduceWithSeed() {

        when:
        def result = runDataflow {
            channel.of(1, 2, 3, 4, 5).reduce (1) { a, e -> a += e }
        }
        then:
        result.getVal() == 16

        when:
        result = runDataflow {
            channel.empty().reduce (10) { a, e -> a += e }
        }
        then:
        result.getVal() == 10
    }

    def testSubscribe() {

        when:
        def count = 0
        runDataflow {
            channel.of(1,2,3,4).subscribe { count++; }
        }
        sleep(100)
        then:
        count == 4
    }

    def testSubscribeComplete() {

        when:
        def count = 0
        def done = false
        runDataflow {
            channel.of(1,2,3).subscribe(onNext: { count++ }, onComplete: { done = true })
        }
        sleep 100
        then:
        done
        count == 3
    }

    def testSubscribeError() {

        when:
        int next=0
        int error=0
        int complete=0
        runDataflow {
            channel.of( 2,1,0,3,3 ).subscribe(
                onNext: { println it/it; next++ },
                onError: { error++ },
                onComplete: { complete++ }
            )
        }
        sleep 100

        then:
        // on the third iteration raise an exception
        // next have to be equals to two
        next == 2
        // error have to be invoked one time
        error == 1
        // complete never
        complete == 0
    }

    def testUnique() {

        when:
        def result = runDataflow {
            channel.of(1,1,1,5,7,7,7,3,3).unique()
        }
        then:
        result.val == 1
        result.val == 5
        result.val == 7
        result.val == 3
        result.val == CH.stop()

        when:
        result = runDataflow {
            channel.of(1,1,1,5,7,7,7,3,3).unique().collect()
        }
        then:
        result.val as Set == [1,5,7,3] as Set

        when:
        result = runDataflow {
            channel.of(1,3,4,5).unique { v -> v % 2 }.collect()
        }
        then:
        result.val as Set == [1,4] as Set

        when:
        result = runDataflow {
            channel.of(1).unique()
        }
        then:
        result.val == 1
        result.val == CH.stop()
    }

    def testUntil() {

        when:
        def result = runDataflow {
            channel.of(1,2,3,4).until { v -> v == 3 }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == CH.stop()

        when:
        result = runDataflow {
            channel.of(1,2,3).until { v -> v == 5 }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == CH.stop()
    }

}
