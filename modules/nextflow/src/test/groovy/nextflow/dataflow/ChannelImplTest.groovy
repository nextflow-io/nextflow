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

import nextflow.Global
import nextflow.Session
import nextflow.dataflow.ChannelNamespace as channel
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.util.HashBag
import spock.lang.Specification
import spock.lang.Timeout

import static nextflow.Nextflow.record
import static nextflow.Nextflow.tuple
import static test.ScriptHelper.runDataflow
import static test.ScriptHelper.runScript
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(5)
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

    def testCombine() {

        when:
        def result = runDataflow {
            def ch1 = channel.of( 1, 2, 3 )
            def ch2 = channel.of( 'x', 'y', 'z' )
            ch1.combine(ch2)
                .collect()
                .flatMap { vals -> vals.toSorted() }
        }
        then:
        result.val == tuple(1, 'x')
        result.val == tuple(1, 'y')
        result.val == tuple(1, 'z')
        result.val == tuple(2, 'x')
        result.val == tuple(2, 'y')
        result.val == tuple(2, 'z')
        result.val == tuple(3, 'x')
        result.val == tuple(3, 'y')
        result.val == tuple(3, 'z')
        result.val == CH.stop()

        when:
        result = runDataflow {
            def ch1 = channel.of( tuple(1, 'x'), tuple(2, 'y'), tuple(3, 'z') )
            def ch2 = channel.of( tuple(4, 13), tuple(4, 21) )
            ch1.combine(ch2)
                .collect()
                .flatMap { vals -> vals.toSorted() }
        }
        then:
        result.val == tuple(1, 'x', 4, 13)
        result.val == tuple(1, 'x', 4, 21)
        result.val == tuple(2, 'y', 4, 13)
        result.val == tuple(2, 'y', 4, 21)
        result.val == tuple(3, 'z', 4, 13)
        result.val == tuple(3, 'z', 4, 21)
        result.val == CH.stop()
    }

    def testCombineNamedArgs() {

        when:
        def result = runDataflow {
            def ch = channel.of( record(id: 1, a: 'x'), record(id: 2, a: 'y'), record(id: 3, a: 'z') )
            ch.combine(b: channel.value('foo'), c: 'bar')
        }
        then:
        result.val == record(id: 1, a: 'x', b: 'foo', c: 'bar')
        result.val == record(id: 2, a: 'y', b: 'foo', c: 'bar')
        result.val == record(id: 3, a: 'z', b: 'foo', c: 'bar')
        result.val == CH.stop()
    }

    def testCombineError() {

        when:
        runDataflow {
            channel.of(1, 2, 3).combine('not-a-channel')
        }
        then:
        thrown(ScriptRuntimeException)

        when:
        runDataflow {
            def ch = channel.of( record(id: 1) )
            ch.combine(b: channel.of(1, 2, 3))
        }
        then:
        thrown(ScriptRuntimeException)
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

        when:
        result = runDataflow {
            channel.of( [1,2], ['a','b'] ).flatMap { vals -> [vals, vals.reverse()] }
        }
        then:
        result.val == [1,2]
        result.val == [2,1]
        result.val == ['a','b']
        result.val == ['b','a']
        result.val == CH.stop()

        when:
        result = runDataflow {
            channel.of( [1,2], ['a',['b','c']] ).flatMap()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 'a'
        result.val == ['b','c']  // <-- nested list are preserved
        result.val == CH.stop()
    }

    def testFlatMapError() {

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of(tuple(1,2)).flatMap()
            }
            '''
        )
        def sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains 'Operator `flatMap` expected an Iterable but received a tuple'
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

        when:
        result = runDataflow {
            channel.of([1, 2, 'a'], [2, 1, 'x'], [1, 2, 'b'])
                .groupBy()
                .collect()
        }
        then:
        def bag = result.val
        bag.size() == 2
        bag.contains( tuple(1, new HashBag<>(['a', 'b'])) )
        bag.contains( tuple(2, new HashBag<>(['x'])) )
    }

    def testGroupByError() {

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of(1, 2, 3).groupBy()
            }
            '''
        )
        def sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains 'Operator `groupBy` expected a 3-tuple of (key, size, value) or a 2-tuple of (key, value)'

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of([1, 1, 'a'], [1, 1, 'b']).groupBy()
            }
            '''
        )
        sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains 'Operator `groupBy` received too many values for grouping key: 1'

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of([1, 2, 'a'], [1, 3, 'b']).groupBy()
            }
            '''
        )
        sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains 'Operator `groupBy` received inconsistent group size for key 1'

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of([1, 3, 'a'], [1, 3, 'b']).groupBy()
            }
            '''
        )
        sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains 'Operator `groupBy` received too few values for grouping keys: 1'
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
        e.message == 'Operator `join` requires the `by` option'

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                left = channel.of(1, 2, 3)
                right = channel.of(1, 2, 3)
                left.join(right, by: 'id')
            }
            '''
        )
        def sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains 'Operator `join` expected a record'
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

    def testMapWithTupleDestructuring() {

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
            c1.mix(c2).mix(c3).collect()
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

    def testMixError() {

        when:
        runDataflow {
            channel.of(1, 2, 3).mix('not-a-channel')
        }
        then:
        thrown(ScriptRuntimeException)
    }

    def testReduce() {

        when:
        def result = runDataflow {
            channel.of(1,2,3,4,5).reduce { acc, v -> acc + v }
        }
        then:
        result.getVal() == 15

        when:
        result = runDataflow {
            channel.of(99).reduce { acc, v -> acc + v }
        }
        then:
        result.getVal() == 99

        when:
        result = runDataflow {
            channel.of(1,2,3,4,5).reduce([]) { acc, v -> acc.add(v) ; acc }
        }
        then:
        result.getVal() == [1, 2, 3, 4, 5]
    }

    def testReduceError() {

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.empty().reduce { acc, v -> acc + v }
            }
            '''
        )
        def sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message.contains "Operator `reduce` received an empty channel with no initial value"
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
        result.val.sort() == [1,3,5,7]

        when:
        result = runDataflow {
            channel.of(1,3,4,5).unique { v -> v % 2 }.collect()
        }
        then:
        result.val.sort() == [1,4]

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

    def testView() {

        when:
        def result = runDataflow {
            channel.of(1,2,3).view()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == CH.stop()

        when:
        result = runDataflow {
            channel.of(1,2,3).view { v -> "value: $v" }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == CH.stop()
    }

    def 'should propagate errors to the session' () {
        given:
        def sess

        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of(1, 2, 3).subscribe { v ->
                    error('failed!')
                }
            }
            '''
        )
        sess = Global.session as Session
        then:
        sess.isAborted()
        sess.error.message == "failed!"
    }

    def 'should fall back to legacy operators' () {
        when:
        def result = runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.of(1, [2, 3], 4).flatten()
            }
            '''
        )
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == 4
        result.val == CH.stop()

        when:
        runDataflow {
            channel.of(1, 2, 3, 4).foo()
        }
        then:
        thrown(MissingMethodException)
    }

}
