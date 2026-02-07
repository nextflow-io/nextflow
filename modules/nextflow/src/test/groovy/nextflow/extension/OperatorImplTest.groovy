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

import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import spock.lang.Specification
import spock.lang.Timeout

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class OperatorImplTest extends Specification {

    def testFilter() {

        when:
        def c1 = runDataflow {
            Channel.of(1,2,3,4,5).filter { it > 3 }
        }
        then:
        c1.val == 4
        c1.val == 5
        c1.val == Channel.STOP

        when:
        def c2 = runDataflow {
            Channel.of('hola','hello','cioa','miao').filter { it =~ /^h.*/ }
        }
        then:
        c2.val == 'hola'
        c2.val == 'hello'
        c2.val == Channel.STOP

        when:
        def c3 = runDataflow {
            Channel.of('hola','hello','cioa','miao').filter { it ==~ /^h.*/ }
        }
        then:
        c3.val == 'hola'
        c3.val == 'hello'
        c3.val == Channel.STOP

        when:
        def c4 = runDataflow {
            Channel.of('hola','hello','cioa','miao').filter( ~/^h.*/ )
        }
        then:
        c4.val == 'hola'
        c4.val == 'hello'
        c4.val == Channel.STOP

        when:
        def c5 = runDataflow {
            Channel.of('hola',1,'cioa',2,3).filter( Number )
        }
        then:
        c5.val == 1
        c5.val == 2
        c5.val == 3
        c5.val == Channel.STOP

        when:
        def result = runDataflow {
            Channel.of(1,2,4,2,4,5,6,7,4).filter(1).count()
        }
        then:
        result.val == 1

        when:
        result = runDataflow {
            Channel.of(1,2,4,2,4,5,6,7,4).filter(2).count()
        }
        then:
        result.val == 2

        when:
        result = runDataflow {
            Channel.of(1,2,4,2,4,5,6,7,4).filter(4).count()
        }
        then:
        result.val == 3

    }

    def testFilterWithValue() {

        when:
        def result = runDataflow {
            Channel.value(3).filter { it>1 }
        }
        then:
        result.val == 3

        when:
        result = runDataflow {
            Channel.value(0).filter { it>1 }
        }
        then:
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.value(Channel.STOP).filter { it>1 }
        }
        then:
        result.val == Channel.STOP
    }

    def testSubscribe() {

        when:
        def count = 0
        runDataflow {
            Channel.of(1,2,3,4).subscribe { count++; }
        }
        sleep(100)
        then:
        count == 4

    }

    def testSubscribe1() {

        when:
        def count = 0
        def done = false
        runDataflow {
            Channel.of(1,2,3).subscribe(onNext: { count++ }, onComplete: { done = true })
        }
        sleep 100
        then:
        done
        count == 3

    }

    def testSubscribe2() {

        when:
        def count = 0
        def done = false
        runDataflow {
            Channel.value(1).subscribe(
                onNext: { count++ },
                onComplete: { done = true }
            )
        }
        sleep 100
        then:
        done
        count == 1

    }

    def testSubscribeError() {

        when:
        int next=0
        int error=0
        int complete=0
        runDataflow {
            Channel.of( 2,1,0,3,3 ).subscribe(
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


    def testMap() {
        when:
        def result = runDataflow {
            Channel.of(1,2,3).map { "Hello $it" }
        }
        then:
        result.val == 'Hello 1'
        result.val == 'Hello 2'
        result.val == 'Hello 3'
        result.val == Channel.STOP
    }

    def testMapWithVariable() {
        when:
        def result = runDataflow {
            Channel.value('Hello').map { it.reverse() }
        }
        then:
        result.val == 'olleH'
        result.val == 'olleH'
        result.val == 'olleH'
    }

    def testMapParamExpanding () {

        when:
        def result = runDataflow {
            Channel.of(1,2,3).map { [it, it] }.map { x, y -> x+y }
        }
        then:
        result.val == 2
        result.val == 4
        result.val == 6
        result.val == Channel.STOP
    }

    def testSkip() {

        when:
        def result = runDataflow {
            Channel.of(1,2,3).map { it == 2 ? Channel.VOID : "Hello $it" }
        }
        then:
        result.val == 'Hello 1'
        result.val == 'Hello 3'
        result.val == Channel.STOP

    }


    def testMapMany () {

        when:
        def result = runDataflow {
            Channel.of(1,2,3).flatMap { v -> [v, v*2] }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 2
        result.val == 4
        result.val == 3
        result.val == 6
        result.val == Channel.STOP
    }

    def testMapManyWithSingleton() {

        when:
        def result = runDataflow {
            Channel.value([1,2,3]).flatMap()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.empty().flatMap()
        }
        then:
        result.val == Channel.STOP

    }

    def testMapManyWithTuples () {

        when:
        def result = runDataflow {
            Channel.of( [1,2], ['a','b'] ).flatMap { vals -> [vals, vals.reverse()] }
        }
        then:
        result.val == [1,2]
        result.val == [2,1]
        result.val == ['a','b']
        result.val == ['b','a']
        result.val == Channel.STOP
    }

    def testMapManyDefault  () {

        when:
        def result = runDataflow {
            Channel.of( [1,2], ['a',['b','c']] ).flatMap()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 'a'
        result.val == ['b','c']  // <-- nested list are preserved
        result.val == Channel.STOP
    }

    def testMapManyWithHashArray () {

        when:
        def result = runDataflow {
            Channel.of(1,2,3).flatMap { n -> [ k: n, v: n*2] }
        }
        then:
        result.val == new MapEntry('k',1)
        result.val == new MapEntry('v',2)
        result.val == new MapEntry('k',2)
        result.val == new MapEntry('v',4)
        result.val == new MapEntry('k',3)
        result.val == new MapEntry('v',6)
        result.val == Channel.STOP

    }



    def testReduce() {

        when:
        def result = runDataflow {
            Channel.of(1,2,3,4,5).reduce { a, e -> a += e }
        }
        then:
        result.getVal() == 15

        when:
        result = runDataflow {
            Channel.of(99).reduce { a, e -> a += e }
        }
        then:
        result.getVal() == 99

        when:
        result = runDataflow {
            Channel.empty().reduce { a, e -> a += e }
        }
        then:
        result.getVal() == null

        when:
        result = runDataflow {
            Channel.of(6,5,4,3,2,1).reduce { a, e -> Channel.STOP }
        }
        then:
        result.val == 6

    }


    def testReduceWithSeed() {

        when:
        def result = runDataflow {
            Channel.of(1, 2, 3, 4, 5).reduce (1) { a, e -> a += e }
        }
        then:
        result.getVal() == 16

        when:
        result = runDataflow {
            Channel.empty().reduce (10) { a, e -> a += e }
        }
        then:
        result.getVal() == 10

        when:
        result = runDataflow {
            Channel.of(6,5,4,3,2,1).reduce(0) { a, e -> a < 3 ? a+1 : Channel.STOP }
        }
        then:
        result.val == 3

    }

    def testFirst() {

        when:
        def result = runDataflow {
            Channel.of(3,6,4,5,4,3,4).first()
        }
        then:
        result.val == 3
    }

    def testFirstWithCriteria() {

        when:
        def result = runDataflow {
            Channel.of(3,6,4,5,4,3,4).first{ it>4 }
        }
        then:
        result.val == 6
    }

    def testFirstWithValue() {

        when:
        def result = runDataflow {
            Channel.value(3).first()
        }
        then:
        result.val == 3

        when:
        result = runDataflow {
            Channel.value(3).first{ it>1 }
        }
        then:
        result.val == 3

        when:
        result = runDataflow {
            Channel.value(3).first{ it>3 }
        }
        then:
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.value(Channel.STOP).first { it>3 }
        }
        then:
        result.val == Channel.STOP
    }


    def testFirstWithCondition() {

        when:
        def result = runDataflow {
            Channel.of(3,6,4,5,4,3,4).first { it % 2 == 0  }
        }
        then:
        result.val == 6

        when:
        result = runDataflow {
            Channel.of( 'a', 'b', 'c', 1, 2 ).first( Number )
        }
        then:
        result.val == 1

        when:
        result = runDataflow {
            Channel.of( 'a', 'b', 1, 2, 'aaa', 'bbb' ).first( ~/aa.*/ )
        }
        then:
        result.val == 'aaa'

        when:
        result = runDataflow {
            Channel.of( 'a', 'b', 1, 2, 'aaa', 'bbb' ).first( 1 )
        }
        then:
        result.val == 1

    }


    def testTake() {

        when:
        def result = runDataflow {
            Channel.of(1,2,3,4,5,6).take(3)
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of(1).take(3)
        }
        then:
        result.val == 1
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of(1,2,3).take(-1)
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of(1,2,3).take(3)
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

    }

    def testLast() {

        when:
        def result = runDataflow {
            Channel.of(3,6,4,5,4,3,9).last()
        }
        then:
        result.val == 9

        when:
        result = runDataflow {
            Channel.value('x').last()
        }
        then:
        result.val == 'x'
    }




    def testCount() {

        when:
        def result = runDataflow {
            Channel.of(4,1,7,5).count()
        }
        then:
        result.val == 4

        when:
        result = runDataflow {
            Channel.of(4,1,7,1,1).count(1)
        }
        then:
        result.val == 3

        when:
        result = runDataflow {
            Channel.of('a','c','c','q','b').count ( ~/c/ ) 
        }
        then:
        result.val == 2

        when:
        result = runDataflow {
            Channel.value(5).count()
        }
        then:
        result.val == 1

        when:
        result = runDataflow {
            Channel.value(5).count(5)
        }
        then:
        result.val == 1

        when:
        result = runDataflow {
            Channel.value(5).count(6)
        }
        then:
        result.val == 0
    }

    def testToList() {

        when:
        def result = runDataflow {
            Channel.of(1,2,3).toList()
        }
        then:
        result.val == [1,2,3]

        when:
        result = runDataflow {
            Channel.value(1).toList()
        }
        then:
        result.val == [1]

        when:
        result = runDataflow {
            Channel.empty().toList()
        }
        then:
        result.val == []
    }

    def testToSortedList() {

        when:
        def result = runDataflow {
            Channel.of(3,1,4,2).toSortedList()
        }
        then:
        result.val == [1,2,3,4]

        when:
        result = runDataflow {
            Channel.of([1,'zeta'], [2,'gamma'], [3,'alpaha'], [4,'delta']).toSortedList { it[1] }
        }
        then:
        result.val == [[3,'alpaha'], [4,'delta'], [2,'gamma'], [1,'zeta'] ]

        when:
        result = runDataflow {
            Channel.value(1).toSortedList()
        }
        then:
        result.val == [1]

        when:
        result = runDataflow {
            Channel.empty().toSortedList()
        }
        then:
        result.val == []

    }


    def testUnique() {

        when:
        def result = runDataflow {
            Channel.of(1,1,1,5,7,7,7,3,3).unique().toList()
        }
        then:
        result.val == [1,5,7,3]

        when:
        result = runDataflow {
            Channel.of(1,3,4,5).unique { it%2 } .toList()
        }
        then:
        result.val == [1,4]

        when:
        result = runDataflow {
            Channel.of(1).unique()
        }
        then:
        result.val == 1

        when:
        result = runDataflow {
            Channel.value(1).unique()
        }
        then:
        result.val == 1
    }

    def testDistinct() {

        when:
        def result = runDataflow {
            Channel.of(1,1,2,2,2,3,1,1,2,2,3).distinct().toList()
        }
        then:
        result.val == [1,2,3,1,2,3]

        when:
        result = runDataflow {
            Channel.of(1,1,2,2,2,3,1,1,2,4,6).distinct { it%2 } .toList()
        }
        then:
        result.val == [1,2,3,2]
    }


    def testFlatten() {

        when:
        def r1 = runDataflow {
            Channel.of(1,2,3).flatten()
        }
        then:
        r1.val == 1
        r1.val == 2
        r1.val == 3
        r1.val == Channel.STOP

        when:
        def r2 = runDataflow {
            Channel.of([1,'a'], [2,'b']).flatten()
        }
        then:
        r2.val == 1
        r2.val == 'a'
        r2.val == 2
        r2.val == 'b'
        r2.val == Channel.STOP

        when:
        def r3 = runDataflow {
            Channel.of( [1,2] as Integer[], [3,4] as Integer[] ).flatten()
        }
        then:
        r3.val == 1
        r3.val == 2
        r3.val == 3
        r3.val == 4
        r3.val == Channel.STOP

        when:
        def r4 = runDataflow {
            Channel.of( [1,[2,3]], 4, [5,[6]] ).flatten()
        }
        then:
        r4.val == 1
        r4.val == 2
        r4.val == 3
        r4.val == 4
        r4.val == 5
        r4.val == 6
        r4.val == Channel.STOP
    }

    def testFlattenWithSingleton() {
        when:
        def result = runDataflow {
            Channel.value([3,2,1]).flatten()
        }
        then:
        result.val == 3
        result.val == 2
        result.val == 1
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.empty().flatten()
        }
        then:
        result.val == Channel.STOP
    }

    def testCollate() {

        when:
        def r1 = runDataflow {
            Channel.of(1,2,3,1,2,3,1).collate( 2, false )
        }
        then:
        r1.val == [1,2]
        r1.val == [3,1]
        r1.val == [2,3]
        r1.val == Channel.STOP

        when:
        def r2 = runDataflow {
            Channel.of(1,2,3,1,2,3,1).collate( 3 )
        }
        then:
        r2.val == [1,2,3]
        r2.val == [1,2,3]
        r2.val == [1]
        r2.val == Channel.STOP

    }

    def testCollateWithStep() {

        when:
        def r1 = runDataflow {
            Channel.of(1,2,3,4).collate( 3, 1, false )
        }
        then:
        r1.val == [1,2,3]
        r1.val == [2,3,4]
        r1.val == Channel.STOP

        when:
        def r2 = runDataflow {
            Channel.of(1,2,3,4).collate( 3, 1, true )
        }
        then:
        r2.val == [1,2,3]
        r2.val == [2,3,4]
        r2.val == [3,4]
        r2.val == [4]
        r2.val == Channel.STOP

        when:
        def r3 = runDataflow {
            Channel.of(1,2,3,4).collate( 3, 1  )
        }
        then:
        r3.val == [1,2,3]
        r3.val == [2,3,4]
        r3.val == [3,4]
        r3.val == [4]
        r3.val == Channel.STOP

        when:
        def r4 = runDataflow {
            Channel.of(1,2,3,4).collate( 4,4 )
        }
        then:
        r4.val == [1,2,3,4]
        r4.val == Channel.STOP

        when:
        def r5 = runDataflow {
            Channel.of(1,2,3,4).collate( 6,6 )
        }
        then:
        r5.val == [1,2,3,4]
        r5.val == Channel.STOP

        when:
        def r6 = runDataflow {
            Channel.of(1,2,3,4).collate( 6,6,false )
        }
        then:
        r6.val == Channel.STOP

    }

    def testCollateIllegalArgs() {
        when:
        Channel.empty().collate(0)
        then:
        thrown(IllegalArgumentException)

        when:
        Channel.empty().collate(-1)
        then:
        thrown(IllegalArgumentException)

        when:
        Channel.empty().collate(0,1)
        then:
        thrown(IllegalArgumentException)

        when:
        Channel.empty().collate(1,0)
        then:
        thrown(IllegalArgumentException)

    }

    def testCollateWithValueChannel() {
        when:
        def result = runDataflow {
            Channel.value(1).collate(1)
        }
        then:
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.value(1).collate(10)
        }
        then:
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.value(1).collate(10, true)
        }
        then:
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.value(1).collate(10, false)
        }
        then:
        result.val == Channel.STOP
    }

    def testMix() {
        when:
        def result = runDataflow {
            def c1 = Channel.of( 1,2,3 )
            def c2 = Channel.of( 'a','b' )
            def c3 = Channel.value( 'z' )
            c1.mix(c2,c3).toList()
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

    def testMixWithSingleton() {
        when:
        def result = runDataflow {
            Channel.value(1).mix( Channel.of(2,3) ).toList()
        }
        then:
        result.val.sort() == [1,2,3]
    }



    def testDefaultMappingClosure() {

        expect:
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9] ) == 7
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9], 2 ) == 9
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [] ) == null
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [], 2 ) == null

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9] as Object[] ) == 7
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9] as Object[], 1 ) == 8
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( ['alpha', 'beta'] as String[] ) == 'alpha'
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( ['alpha', 'beta'] as String[], 1 ) == 'beta'

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet ) == 6
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet, 1 ) == 7
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet, 2 ) == 8
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet, 5 ) == null

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9] ) == 1
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9], 1 ) == 2
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9], 2 ) == 9
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9], 3 ) == null

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(0) ) == 'a'
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(0), 1 ) == 1
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(0), 2 ) == null

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(1) ) == 'b'
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(1), 1 ) == 2
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(1), 2 ) == null

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( [:] ) == null

        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( 99 ) == 99
        OperatorImpl.DEFAULT_MAPPING_CLOSURE.call( 99, 2 ) == null

    }


    def testCross() {

        when:
        def result = runDataflow {
            def ch1 = Channel.of( [1, 'x'], [2,'y'], [3,'z'] )
            def ch2 = Channel.of( [1,11], [1,13], [2,21],[2,22], [2,23], [4,1], [4,2]  )
            ch1.cross(ch2)
        }

        then:
        result.val == [ [1, 'x'], [1,11] ]
        result.val == [ [1, 'x'], [1,13] ]
        result.val == [ [2, 'y'], [2,21] ]
        result.val == [ [2, 'y'], [2,22] ]
        result.val == [ [2, 'y'], [2,23] ]
        result.val == Channel.STOP

    }

    def testCross2() {

        when:
        def result = runDataflow {
            def ch1 = Channel.of(['PF00006', 'PF00006.sp_lib'])
            def ch2 = Channel.of(['PF00006', 'PF00006_mafft.aln'], ['PF00006', 'PF00006_clustalo.aln'])
            ch1.cross(ch2)
        }

        then:
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_mafft.aln'] ]
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_clustalo.aln'] ]
        result.val == Channel.STOP

    }


    def testCross3() {

        when:
        def result = runDataflow {
            def ch1 = Channel.of(['PF00006', 'PF00006.sp_lib'])
            def ch2 = Channel.of(['PF00006', 'PF00006_mafft.aln'], ['PF00006', 'PF00006_clustalo.aln'])
            ch1.cross(ch2)
        }

        then:
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_mafft.aln'] ]
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_clustalo.aln'] ]
        result.val == Channel.STOP

    }


    def testConcat() {

        when:
        def all = runDataflow {
            def c1 = Channel.of(1,2,3)
            def c2 = Channel.of('a','b','c')
            c1.concat(c2)
        }
        then:
        all.val == 1
        all.val == 2
        all.val == 3
        all.val == 'a'
        all.val == 'b'
        all.val == 'c'
        all.val == Channel.STOP

        when:
        def result = runDataflow {
            def d1 = Channel.of([1, 2]).flatMap { vals -> sleep 20; vals }
            def d2 = Channel.of('a','b','c')
            def d3 = Channel.of(['p', 'q']).flatMap { vals -> sleep 100; vals }
            d1.concat(d2,d3)
        }

        then:
        result.val == 1
        result.val == 2
        result.val == 'a'
        result.val == 'b'
        result.val == 'c'
        result.val == 'p'
        result.val == 'q'
        result.val == Channel.STOP

    }

    def testContactWithSingleton() {
        when:
        def result = runDataflow {
            Channel.value(1).concat( Channel.of(2,3) )
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
    }


    def testGroupTuple() {

        when:
        def result = runDataflow {
            Channel.of([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'c'], [2, 'y'], [3, 'q'])
                .groupTuple()
        }

        then:
        result.val == [1, ['a', 'b','c'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['q', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithCount() {

        when:
        def result = runDataflow {
            Channel.of([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'d'], [1,'c'], [2, 'y'], [1,'f'])
                .groupTuple(size: 2)
        }

        then:
        result.val == [1, ['a', 'b'] ]
        result.val == [1, ['d', 'c'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'d'], [1,'c'], [2, 'y'], [1,'f'])
                .groupTuple(size: 2, remainder: true)
        }

        then:
        result.val == [1, ['a', 'b'] ]
        result.val == [1, ['d', 'c'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['q']]
        result.val == [1, ['f']]
        result.val == Channel.STOP

    }

    def testGroupTupleWithSortNatural() {

        when:
        def result = runDataflow {
            Channel.of([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: true)
        }

        then:
        result.val == [1, ['a', 'b','c','w','z'] ]
        result.val == [2, ['x','y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: 'natural')
        }

        then:
        result.val == [1, ['a', 'b','c','w','z'] ]
        result.val == [2, ['x','y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }


    def testGroupTupleWithSortHash() {

        when:
        def result = runDataflow {
            Channel.of([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: 'hash')
        }

        then:
        result.val == [1, ['a', 'c','z','b','w'] ]
        result.val == [2, ['y','x'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithComparator() {

        when:
        def result = runDataflow {
            Channel.of([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: { o1, o2 -> o2<=>o1 } as Comparator )
        }

        then:
        result.val == [1, ['z','w','c','b','a'] ]
        result.val == [2, ['y','x'] ]
        result.val == [3, ['q','p'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithClosureWithSingle() {

        when:
        def result = runDataflow {
            Channel.of([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: { it } )
        }

        then:
        result.val == [1, ['a', 'b','c','w','z'] ]
        result.val == [2, ['x','y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithComparatorWithPair() {

        when:
        def result = runDataflow {
            Channel.of([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: { o1, o2 -> o2<=>o1 } )
        }

        then:
        result.val == [1, ['z','w','c','b','a'] ]
        result.val == [2, ['y','x'] ]
        result.val == [3, ['q','p'] ]
        result.val == Channel.STOP

    }


    def testGroupTupleWithIndex () {

        given:
        def file1 = Paths.get('/path/file_1')
        def file2 = Paths.get('/path/file_2')
        def file3 = Paths.get('/path/file_3')

        when:
        def result = runDataflow {
            Channel.of([1,'a', file1], [1,'b',file2], [2,'x',file2], [3, 'q',file1], [1,'c',file3], [2, 'y',file3], [3, 'q',file1])
                .groupTuple(by: 2)
        }

        then:
        result.val == [ [1,3,3], ['a','q','q'], file1 ]
        result.val == [ [1,2], ['b','x'], file2 ]
        result.val == [ [1,2], ['c','y'], file3 ]
        result.val == Channel.STOP


        when:
        result = runDataflow {
            Channel.of([1,'a', file1], [1,'b',file2], [2,'x',file2], [3, 'q',file1], [1,'c',file3], [2, 'y',file3], [3, 'q',file1])
                .groupTuple(by: [2])
        }

        then:
        result.val == [ [1,3,3], ['a','q','q'], file1 ]
        result.val == [ [1,2], ['b','x'], file2 ]
        result.val == [ [1,2], ['c','y'], file3 ]
        result.val == Channel.STOP


        when:
        result = runDataflow {
            Channel.of([1,'a', file1], [1,'b',file2], [2,'x',file2], [1, 'q',file1], [3, 'y', file3], [1,'c',file2], [2, 'y',file2], [3, 'q',file1], [1, 'z', file2], [3, 'c', file3])
                .groupTuple(by: [0,2])
        }

        then:
        result.val == [ 1, ['a','q'], file1 ]
        result.val == [ 1, ['b','c','z'], file2 ]
        result.val == [ 2, ['x','y'], file2 ]
        result.val == [ 3, ['y','c'], file3 ]
        result.val == [ 3, ['q'], file1 ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithNotMatchingCardinality() {

        when:
        def result = runDataflow {
            Channel.of(
                    [1,'a'],
                    [1,'b'],
                    [2,'x'],
                    [3,'p'],
                    [1,'c','d'],
                    [2,'y'],
                    [3,'q']
                )
                .groupTuple()
        }

        then:
        result.val == [1, ['a', 'b', 'c'], ['d'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithNotMatchingCardinalityAndFixedSize() {

        when:
        def result = runDataflow {
            Channel.of(
                    [1,'a'],
                    [1,'b'],
                    [2,'x'],
                    [3,'p'],
                    [1,'c','d'],
                    [2,'y'],
                    [3,'q']
                )
                .groupTuple(size:2)
        }

        then:
        result.val == [1, ['a', 'b'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP
    }

    def testGroupTupleWithNotMatchingCardinalityAndFixedSizeAndRemainder() {

        when:
        def result = runDataflow {
            Channel.of(
                    [1,'a'],
                    [1,'b'],
                    [2,'x'],
                    [3,'p'],
                    [1,'c','d'],
                    [2, 'y'],
                    [3, 'q']
                )
                .groupTuple(size:2, remainder: true)
        }

        then:
        result.val == [1, ['a', 'b'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == [1, ['c'], ['d']]
        result.val == Channel.STOP
    }

    def testChannelIfEmpty() {

        def result

        when:
        result = runDataflow {
            Channel.of(1,2,3).ifEmpty(100)
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.empty().ifEmpty(100)
        }
        then:
        result.val == 100
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.empty().ifEmpty { 1+2 }
        }
        then:
        result.val == 3
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.value(1).ifEmpty(100)
        }
        then:
        result instanceof DataflowVariable
        result.val == 1

        when:
        result = runDataflow {
            Channel.empty().ifEmpty(100)
        }
        then:
        !(result instanceof DataflowVariable)
        result.val == 100

    }

    def 'should create a channel given a list'() {

        when:
        def result = runDataflow {
            [10,20,30].channel()
        }
        then:
        result.val == 10
        result.val == 20
        result.val == 30
        result.val == Channel.STOP

    }


    def 'should assign a channel to new variable' () {

        when:
        def result = runDataflow {
            Channel.of(10,20,30)
                .map { it +2 }
                .set { ch_result }
            ch_result
        }

        then:
        result.val == 12
        result.val == 22
        result.val == 32
        result.val == Channel.STOP
    }

    def 'should always the same value' () {

        when:
        def x = runDataflow {
            Channel.value('Hello')
        }
        then:
        x.val == 'Hello'
        x.val == 'Hello'
        x.val == 'Hello'
    }

    def 'should emit channel items until the condition is verified' () {

        when:
        def result = runDataflow {
            Channel.of(1,2,3,4).until { it == 3 }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of(1,2,3).until { it == 5 }
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

    }


    def 'should assign singleton channel to a new variable' () {
        when:
        def result = runDataflow {
            Channel.value('Hello').set { ch_result }
            ch_result
        }

        then:
        result.val == 'Hello'
        result.val == 'Hello'
        result.val == 'Hello'

    }

    def 'should assign queue channel to a new variable' () {
        when:
        def result = runDataflow {
            Channel.of(1,2,3).set { ch_result }
            ch_result
        }

        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
    }

}
