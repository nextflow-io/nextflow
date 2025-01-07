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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import nextflow.prov.Tracker
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class OperatorImplTest extends Specification {

    def setupSpec() {
        new Session()
    }

    private mval(Object obj) {
        if( obj instanceof DataflowReadChannel ) {
            def result = obj.val
            return result instanceof Tracker.Msg ? result.value : result
        }
        else
            return obj
    }

    def testFilter() {

        when:
        def c1 = Channel.of(1,2,3,4,5).filter { it > 3 }
        then:
        mval(c1) == 4
        mval(c1) == 5
        mval(c1) == Channel.STOP

        when:
        def c2 = Channel.of('hola','hello','cioa','miao').filter { it =~ /^h.*/ }
        then:
        mval(c2) == 'hola'
        mval(c2) == 'hello'
        mval(c2) == Channel.STOP

        when:
        def c3 = Channel.of('hola','hello','cioa','miao').filter { it ==~ /^h.*/ }
        then:
        mval(c3) == 'hola'
        mval(c3) == 'hello'
        mval(c3) == Channel.STOP

        when:
        def c4 = Channel.of('hola','hello','cioa','miao').filter( ~/^h.*/ )
        then:
        mval(c4) == 'hola'
        mval(c4) == 'hello'
        mval(c4) == Channel.STOP

        when:
        def c5 = Channel.of('hola',1,'cioa',2,3).filter( Number )
        then:
        mval(c5) == 1
        mval(c5) == 2
        mval(c5) == 3
        mval(c5) == Channel.STOP

        expect:
        mval( Channel.of(1,2,4,2,4,5,6,7,4).filter(1) .count() ) == 1
        mval( Channel.of(1,2,4,2,4,5,6,7,4).filter(2) .count() ) == 2
        mval( Channel.of(1,2,4,2,4,5,6,7,4).filter(4) .count() ) == 3

    }

    def testFilterWithValue() {
        expect:
        mval(Channel.value(3).filter { it>1 }) == 3
        mval(Channel.value(0).filter { it>1 }) == Channel.STOP
        mval(Channel.value(Channel.STOP).filter { it>1 }) == Channel.STOP
    }

    def testSubscribe() {

        when:
        def channel = Channel.create()
        int count = 0
        channel.subscribe { count++; } << 1 << 2 << 3
        sleep(100)
        then:
        count == 3

        when:
        count = 0
        channel = Channel.of(1,2,3,4)
        channel.subscribe { count++; }
        sleep(100)
        then:
        count == 4


    }

    def testSubscribe1() {

        when:
        def count = 0
        def done = false
        Channel.of(1,2,3).subscribe onNext:  { count++ }, onComplete: { done = true }
        sleep 100
        then:
        done
        count == 3

    }

    def testSubscribe2() {

        when:
        def count = 0
        def done = false
        Channel.value(1).subscribe onNext:  { count++ }, onComplete: { done = true }
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
        Channel
                .from( 2,1,0,3,3 )
                .subscribe onNext: { println it/it; next++ }, onError: { error++ }, onComplete: { complete++ }
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
        def result = Channel.of(1,2,3).map { "Hello $it" }
        then:
        mval(result) == 'Hello 1'
        mval(result) == 'Hello 2'
        mval(result) == 'Hello 3'
        mval(result) == Channel.STOP
    }

    def testMapWithVariable() {
        given:
        def variable = Channel.value('Hello')
        when:
        def result = variable.map { it.reverse() }
        then:
        mval(result) == 'olleH'
        mval(result) == 'olleH'
        mval(result) == 'olleH'
    }

    def testMapParamExpanding () {

        when:
        def result = Channel.of(1,2,3).map { [it, it] }.map { x, y -> x+y }
        then:
        mval(result) == 2
        mval(result) == 4
        mval(result) == 6
        mval(result) == Channel.STOP
    }

    def testSkip() {

        when:
        def result = Channel.of(1,2,3).map { it == 2 ? Channel.VOID : "Hello $it" }
        then:
        mval(result) == 'Hello 1'
        mval(result) == 'Hello 3'
        mval(result) == Channel.STOP

    }


    def testMapMany () {

        when:
        def result = Channel.of(1,2,3).flatMap { it -> [it, it*2] }
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 2
        mval(result) == 4
        mval(result) == 3
        mval(result) == 6
        mval(result) == Channel.STOP
    }

    def testMapManyWithSingleton() {

        when:
        def result = Channel.value([1,2,3]).flatMap()
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 3
        mval(result) == Channel.STOP

        when:
        result = Channel.empty().flatMap()
        then:
        mval(result) == Channel.STOP

    }

    def testMapManyWithTuples () {

        when:
        def result = Channel.of( [1,2], ['a','b'] ).flatMap { it -> [it, it.reverse()] }
        then:
        mval(result) == [1, 2]
        mval(result) == [2, 1]
        mval(result) == ['a', 'b']
        mval(result) == ['b', 'a']
        mval(result) == Channel.STOP
    }

    def testMapManyDefault  () {

        when:
        def result = Channel.of( [1,2], ['a',['b','c']] ).flatMap()
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 'a'
        mval(result) == ['b', 'c']  // <-- nested list are preserved
        mval(result) == Channel.STOP
    }

    def testMapManyWithHashArray () {

        when:
        def result = Channel.of(1,2,3).flatMap { it -> [ k: it, v: it*2] }
        then:
        mval(result) == new MapEntry('k',1)
        mval(result) == new MapEntry('v',2)
        mval(result) == new MapEntry('k',2)
        mval(result) == new MapEntry('v',4)
        mval(result) == new MapEntry('k',3)
        mval(result) == new MapEntry('v',6)
        mval(result) == Channel.STOP

    }



    def testReduce() {

        when:
        def channel = Channel.create()
        def result = channel.reduce { a, e -> a += e }
        channel << 1 << 2 << 3 << 4 << 5 << Channel.STOP
        then:
        mval(result) == 15


        when:
        channel = Channel.of(1,2,3,4,5)
        result = channel.reduce { a, e -> a += e }
        then:
        mval(result) == 15

        when:
        channel = Channel.create()
        result = channel.reduce { a, e -> a += e }
        channel << 99 << Channel.STOP
        then:
        mval(result) == 99

        when:
        channel = Channel.create()
        result = channel.reduce { a, e -> a += e }
        channel << Channel.STOP
        then:
        mval(result) == null

        when:
        result = Channel.of(6,5,4,3,2,1).reduce { a, e -> Channel.STOP }
        then:
        mval(result) == 6

    }


    def testReduceWithSeed() {

        when:
        def channel = Channel.create()
        def result = channel.reduce (1) { a, e -> a += e }
        channel << 1 << 2 << 3 << 4 << 5 << Channel.STOP
        then:
        mval(result) == 16

        when:
        channel = Channel.create()
        result = channel.reduce (10) { a, e -> a += e }
        channel << Channel.STOP
        then:
        mval(result) == 10

        when:
        result = Channel.of(6,5,4,3,2,1).reduce(0) { a, e -> a < 3 ? a+1 : Channel.STOP }
        then:
        mval(result) == 3

    }

    def testFirst() {

        expect:
        mval(Channel.of(3,6,4,5,4,3,4).first()) == 3
    }

    def testFirstWithCriteria() {
        expect:
        mval(Channel.of(3,6,4,5,4,3,4).first{ it>4 }) == 6
    }

    def testFirstWithValue() {

        expect:
        mval(Channel.value(3).first()) == 3
        mval(Channel.value(3).first{ it>1 }) == 3
        mval(Channel.value(3).first{ it>3 }) == Channel.STOP
        mval(Channel.value(Channel.STOP).first { it>3 }) == Channel.STOP
    }


    def testFirstWithCondition() {

        expect:
        mval(Channel.of(3,6,4,5,4,3,4).first { it % 2 == 0  }) == 6
        mval(Channel.of( 'a', 'b', 'c', 1, 2 ).first( Number )) == 1
        mval(Channel.of( 'a', 'b', 1, 2, 'aaa', 'bbb' ).first( ~/aa.*/ )) == 'aaa'
        mval(Channel.of( 'a', 'b', 1, 2, 'aaa', 'bbb' ).first( 1 )) == 1

    }


    def testTake() {

        when:
        def result = Channel.of(1,2,3,4,5,6).take(3)
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 3
        mval(result) == Channel.STOP

        when:
        result = Channel.of(1).take(3)
        then:
        mval(result) == 1
        mval(result) == Channel.STOP

        when:
        result = Channel.of(1,2,3).take(0)
        then:
        mval(result) == Channel.STOP

        when:
        result = Channel.of(1,2,3).take(-1)
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 3
        mval(result) == Channel.STOP

        when:
        result = Channel.of(1,2,3).take(3)
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 3
        mval(result) == Channel.STOP

    }

    def testLast() {

        expect:
        mval(Channel.of(3,6,4,5,4,3,9).last()) == 9
        mval(Channel.value('x').last()) == 'x'
    }




    def testCount() {
        expect:
        mval(Channel.of(4,1,7,5).count()) == 4
        mval(Channel.of(4,1,7,1,1).count(1)) == 3
        mval(Channel.of('a','c','c','q','b').count ( ~/c/ )) == 2
        mval(Channel.value(5).count()) == 1
        mval(Channel.value(5).count(5)) == 1
        mval(Channel.value(5).count(6)) == 0
    }

    def testToList() {

        when:
        def channel = Channel.of(1,2,3)
        then:
        mval(channel.toList()) == [1, 2, 3]

        when:
        channel = Channel.create()
        channel << Channel.STOP
        then:
        mval(channel.toList()) == []

        when:
        channel = Channel.value(1)
        then:
        mval(channel.toList()) == [1]

        when:
        channel = Channel.empty()
        then:
        mval(channel.toList()) == []
    }

    def testToSortedList() {

        when:
        def channel = Channel.of(3,1,4,2)
        then:
        mval(channel.toSortedList()) == [1, 2, 3, 4]

        when:
        channel = Channel.empty()
        then:
        mval(channel.toSortedList()) == []

        when:
        channel = Channel.of([1,'zeta'], [2,'gamma'], [3,'alpaha'], [4,'delta'])
        then:
        mval(channel.toSortedList { it[1] }) == [[3, 'alpaha'], [4, 'delta'], [2, 'gamma'], [1, 'zeta'] ]

        when:
        channel = Channel.value(1)
        then:
        mval(channel.toSortedList()) == [1]

        when:
        channel = Channel.empty()
        then:
        mval(channel.toSortedList()) == []

    }

    def testUnique() {
        expect:
        mval(Channel.of(1,1,1,5,7,7,7,3,3).unique().toList()) == [1, 5, 7, 3]
        mval(Channel.of(1,3,4,5).unique { it%2 } .toList()) == [1, 4]
        and:
        mval(Channel.of(1).unique()) == 1
        mval(Channel.value(1).unique()) == 1
    }

    def testDistinct() {
        expect:
        mval(Channel.of(1,1,2,2,2,3,1,1,2,2,3).distinct().toList()) == [1, 2, 3, 1, 2, 3]
        mval(Channel.of(1,1,2,2,2,3,1,1,2,4,6).distinct { it%2 } .toList()) == [1, 2, 3, 2]
    }


    def testFlatten() {

        when:
        def r1 = Channel.of(1,2,3).flatten()
        then:
        mval(r1) == 1
        mval(r1) == 2
        mval(r1) == 3
        mval(r1) == Channel.STOP

        when:
        def r2 = Channel.of([1,'a'], [2,'b']).flatten()
        then:
        mval(r2) == 1
        mval(r2) == 'a'
        mval(r2) == 2
        mval(r2) == 'b'
        mval(r2) == Channel.STOP

        when:
        def r3 = Channel.of( [1,2] as Integer[], [3,4] as Integer[] ).flatten()
        then:
        mval(r3) == 1
        mval(r3) == 2
        mval(r3) == 3
        mval(r3) == 4
        mval(r3) == Channel.STOP

        when:
        def r4 = Channel.of( [1,[2,3]], 4, [5,[6]] ).flatten()
        then:
        mval(r4) == 1
        mval(r4) == 2
        mval(r4) == 3
        mval(r4) == 4
        mval(r4) == 5
        mval(r4) == 6
        mval(r4) == Channel.STOP
    }

    def testFlattenWithSingleton() {
        when:
        def result = Channel.value([3,2,1]).flatten()
        then:
        mval(result) == 3
        mval(result) == 2
        mval(result) == 1
        mval(result) == Channel.STOP

        when:
        result = Channel.empty().flatten()
        then:
        mval(result) ==  Channel.STOP
    }

    def testCollate() {

        when:
        def r1 = Channel.of(1,2,3,1,2,3,1).collate( 2, false )
        then:
        mval(r1) == [1, 2]
        mval(r1) == [3, 1]
        mval(r1) == [2, 3]
        mval(r1) == Channel.STOP

        when:
        def r2 = Channel.of(1,2,3,1,2,3,1).collate( 3 )
        then:
        mval(r2) == [1, 2, 3]
        mval(r2) == [1, 2, 3]
        mval(r2) == [1]
        mval(r2) == Channel.STOP

    }

    def testCollateWithStep() {

        when:
        def r1 = Channel.of(1,2,3,4).collate( 3, 1, false )
        then:
        mval(r1) == [1, 2, 3]
        mval(r1) == [2, 3, 4]
        mval(r1) == Channel.STOP

        when:
        def r2 = Channel.of(1,2,3,4).collate( 3, 1, true )
        then:
        mval(r2) == [1, 2, 3]
        mval(r2) == [2, 3, 4]
        mval(r2) == [3, 4]
        mval(r2) == [4]
        mval(r2) == Channel.STOP

        when:
        def r3 = Channel.of(1,2,3,4).collate( 3, 1  )
        then:
        mval(r3) == [1, 2, 3]
        mval(r3) == [2, 3, 4]
        mval(r3) == [3, 4]
        mval(r3) == [4]
        mval(r3) == Channel.STOP

        when:
        def r4 = Channel.of(1,2,3,4).collate( 4,4 )
        then:
        mval(r4) == [1, 2, 3, 4]
        mval(r4) == Channel.STOP

        when:
        def r5 = Channel.of(1,2,3,4).collate( 6,6 )
        then:
        mval(r5) == [1, 2, 3, 4]
        mval(r5) == Channel.STOP

        when:
        def r6 = Channel.of(1,2,3,4).collate( 6,6,false )
        then:
        mval(r6) == Channel.STOP

    }

    def testCollateIllegalArgs() {
        when:
        Channel.create().collate(0)
        then:
        thrown(IllegalArgumentException)

        when:
        Channel.create().collate(-1)
        then:
        thrown(IllegalArgumentException)

        when:
        Channel.create().collate(0,1)
        then:
        thrown(IllegalArgumentException)

        when:
        Channel.create().collate(1,0)
        then:
        thrown(IllegalArgumentException)

    }

    def testCollateWithValueChannel() {
        when:
        def result = Channel.value(1).collate(1)
        then:
        mval(result) == [1]
        mval(result) == Channel.STOP

        when:
        result = Channel.value(1).collate(10)
        then:
        mval(result) == [1]
        mval(result) == Channel.STOP

        when:
        result = Channel.value(1).collate(10, true)
        then:
        mval(result) == [1]
        mval(result) == Channel.STOP

        when:
        result = Channel.value(1).collate(10, false)
        then:
        mval(result) == Channel.STOP
    }

    def testMix() {
        when:
        def c1 = Channel.of( 1,2,3 )
        def c2 = Channel.of( 'a','b' )
        def c3 = Channel.value( 'z' )
        def result = c1.mix(c2,c3).toList().val

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
        def result = Channel.value(1).mix( Channel.of(2,3)  )
        then:
        mval(result.toList().sort()) == [1, 2, 3]
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

        setup:
        def ch1 = Channel.of( [1, 'x'], [2,'y'], [3,'z'] )
        def ch2 = Channel.of( [1,11], [1,13], [2,21],[2,22], [2,23], [4,1], [4,2]  )

        when:
        def result = ch1.cross(ch2)

        then:
        mval(result) == [[1, 'x'], [1, 11] ]
        mval(result) == [[1, 'x'], [1, 13] ]
        mval(result) == [[2, 'y'], [2, 21] ]
        mval(result) == [[2, 'y'], [2, 22] ]
        mval(result) == [[2, 'y'], [2, 23] ]
        mval(result) == Channel.STOP

    }

    def testCross2() {

        setup:
        def ch1 = Channel.create()
        def ch2 = Channel.of ( ['PF00006', 'PF00006_mafft.aln'], ['PF00006', 'PF00006_clustalo.aln'])

        when:
        Thread.start {  sleep 100;   ch1 << ['PF00006', 'PF00006.sp_lib'] << Channel.STOP }
        def result = ch1.cross(ch2)

        then:
        mval(result) == [['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_mafft.aln'] ]
        mval(result) == [['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_clustalo.aln'] ]
        mval(result) == Channel.STOP

    }


    def testCross3() {

        setup:
        def ch1 = Channel.of(['PF00006', 'PF00006.sp_lib'])
        def ch2 = Channel.create ( )

        when:
        Thread.start {  sleep 100;  ch2 << ['PF00006', 'PF00006_mafft.aln'] <<  ['PF00006', 'PF00006_clustalo.aln']<< Channel.STOP }
        def result = ch1.cross(ch2)

        then:
        mval(result) == [['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_mafft.aln'] ]
        mval(result) == [['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_clustalo.aln'] ]
        mval(result) == Channel.STOP

    }


    def testConcat() {

        when:
        def c1 = Channel.of(1,2,3)
        def c2 = Channel.of('a','b','c')
        def all = c1.concat(c2)
        then:
        mval(all) == 1
        mval(all) == 2
        mval(all) == 3
        mval(all) == 'a'
        mval(all) == 'b'
        mval(all) == 'c'
        mval(all) == Channel.STOP

        when:
        def d1 = Channel.create()
        def d2 = Channel.of('a','b','c')
        def d3 = Channel.create()
        def result = d1.concat(d2,d3)

        Thread.start { sleep 20; d3 << 'p' << 'q' << Channel.STOP }
        Thread.start { sleep 100; d1 << 1 << 2 << Channel.STOP }

        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 'a'
        mval(result) == 'b'
        mval(result) == 'c'
        mval(result) == 'p'
        mval(result) == 'q'
        mval(result) == Channel.STOP

    }

    def testContactWithSingleton() {
        when:
        def result = Channel.value(1).concat( Channel.of(2,3) )
        then:
        mval(result) == 1
        mval(result) == 2
        mval(result) == 3
        mval(result) == Channel.STOP
    }


    def testGroupTuple() {

        when:
        def result = Channel
                        .from([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'c'], [2, 'y'], [3, 'q'])
                        .groupTuple()

        then:
        result.val == [1, ['a', 'b','c'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['q', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithCount() {

        when:
        def result = Channel
                .from([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'d'], [1,'c'], [2, 'y'], [1,'f'])
                .groupTuple(size: 2)

        then:
        result.val == [1, ['a', 'b'] ]
        result.val == [1, ['d', 'c'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == Channel.STOP

        when:
        result = Channel
                .from([1,'a'], [1,'b'], [2,'x'], [3, 'q'], [1,'d'], [1,'c'], [2, 'y'], [1,'f'])
                .groupTuple(size: 2, remainder: true)

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
        def result = Channel
                .from([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: true)

        then:
        result.val == [1, ['a', 'b','c','w','z'] ]
        result.val == [2, ['x','y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

        when:
        result = Channel
                .from([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: 'natural')

        then:
        result.val == [1, ['a', 'b','c','w','z'] ]
        result.val == [2, ['x','y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }


    def testGroupTupleWithSortHash() {

        when:
        def result = Channel
                .from([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: 'hash')

        then:
        result.val == [1, ['a', 'c','z','b','w'] ]
        result.val == [2, ['y','x'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithComparator() {

        when:
        def result = Channel
                .from([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: { o1, o2 -> o2<=>o1 } as Comparator )

        then:
        result.val == [1, ['z','w','c','b','a'] ]
        result.val == [2, ['y','x'] ]
        result.val == [3, ['q','p'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithClosureWithSingle() {

        when:
        def result = Channel
                .from([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: { it -> it } )

        then:
        result.val == [1, ['a', 'b','c','w','z'] ]
        result.val == [2, ['x','y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithComparatorWithPair() {

        when:
        def result = Channel
                .from([1,'z'], [1,'w'], [1,'a'], [1,'b'], [2, 'y'], [2,'x'], [3, 'q'], [1,'c'], [3, 'p'])
                .groupTuple(sort: { o1, o2 -> o2<=>o1 } )

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
        def result = Channel
                .from([1,'a', file1], [1,'b',file2], [2,'x',file2], [3, 'q',file1], [1,'c',file3], [2, 'y',file3], [3, 'q',file1])
                .groupTuple(by: 2)

        then:
        result.val == [ [1,3,3], ['a','q','q'], file1 ]
        result.val == [ [1,2], ['b','x'], file2 ]
        result.val == [ [1,2], ['c','y'], file3 ]
        result.val == Channel.STOP


        when:
        result = Channel
                .from([1,'a', file1], [1,'b',file2], [2,'x',file2], [3, 'q',file1], [1,'c',file3], [2, 'y',file3], [3, 'q',file1])
                .groupTuple(by: [2])

        then:
        result.val == [ [1,3,3], ['a','q','q'], file1 ]
        result.val == [ [1,2], ['b','x'], file2 ]
        result.val == [ [1,2], ['c','y'], file3 ]
        result.val == Channel.STOP


        when:
        result = Channel
                .from([1,'a', file1], [1,'b',file2], [2,'x',file2], [1, 'q',file1], [3, 'y', file3], [1,'c',file2], [2, 'y',file2], [3, 'q',file1], [1, 'z', file2], [3, 'c', file3])
                .groupTuple(by: [0,2])

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
        def result = Channel
                .of([1,'a'],
                    [1,'b'],
                    [2,'x'],
                    [3,'p'],
                    [1,'c','d'],
                    [2,'y'],
                    [3,'q'])
                .groupTuple()

        then:
        result.val == [1, ['a', 'b', 'c'], ['d'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP

    }

    def testGroupTupleWithNotMatchingCardinalityAndFixedSize() {

        when:
        def result = Channel
                .of([1,'a'],
                    [1,'b'],
                    [2,'x'],
                    [3,'p'],
                    [1,'c','d'],
                    [2,'y'],
                    [3,'q'])
                .groupTuple(size:2)

        then:
        result.val == [1, ['a', 'b'] ]
        result.val == [2, ['x', 'y'] ]
        result.val == [3, ['p', 'q'] ]
        result.val == Channel.STOP
    }

    def testGroupTupleWithNotMatchingCardinalityAndFixedSizeAndRemainder() {

        when:
        def result = Channel
                .of([1,'a'],
                    [1,'b'],
                    [2,'x'],
                    [3,'p'],
                    [1,'c','d'],
                    [2, 'y'],
                    [3, 'q'])
                .groupTuple(size:2, remainder: true)

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
        result = Channel.of(1,2,3).ifEmpty(100)
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = Channel.empty().ifEmpty(100)
        then:
        result.val == 100
        result.val == Channel.STOP

        when:
        result = Channel.empty().ifEmpty { 1+2  }
        then:
        result.val == 3
        result.val == Channel.STOP

        when:
        result = Channel.value(1).ifEmpty(100)
        then:
        result instanceof DataflowVariable
        result.val == 1

        when:
        result = Channel.empty().ifEmpty(100)
        then:
        result instanceof DataflowQueue
        result.val == 100

    }

    def 'should create a channel given a list'() {

        when:
        def result = [10,20,30].channel()
        then:
        result.val == 10
        result.val == 20
        result.val == 30
        result.val == Channel.STOP

    }


    def 'should assign a channel to new variable' () {
        given:
        def session = new Session()

        when:
        Channel.of(10,20,30)
                .map { it +2 }
                .set { result }

        then:
        session.binding.result.val == 12
        session.binding.result.val == 22
        session.binding.result.val == 32
        session.binding.result.val == Channel.STOP

    }

    def 'should always the same value' () {

        when:
        def x = Channel.value('Hello')
        then:
        x.val == 'Hello'
        x.val == 'Hello'
        x.val == 'Hello'
    }

    def 'should emit channel items until the condition is verified' () {

        when:
        def result = Channel.of(1,2,3,4).until { it == 3 }
        then:
        result.val == 1
        result.val == 2
        result.val == Channel.STOP

        when:
        result = Channel.of(1,2,3).until { it == 5 }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

    }


    def 'should assign singleton channel to a new variable' () {
        given:
        def session = new Session()

        when:
        Channel.value('Hello').set { result }

        then:
        session.binding.result.val == 'Hello'
        session.binding.result.val == 'Hello'
        session.binding.result.val == 'Hello'

    }

    def 'should assign queue channel to a new variable' () {
        given:
        def session = new Session()

        when:
        Channel.of(1,2,3).set { result }

        then:
        session.binding.result.val == 1
        session.binding.result.val == 2
        session.binding.result.val == 3
        session.binding.result.val == Channel.STOP
    }

}
