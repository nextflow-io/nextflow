/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class OperatorExTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should check dataflow read channel' () {
        expect:
        OperatorEx.isReadChannel(DataflowVariable.class)
        OperatorEx.isReadChannel(DataflowQueue.class)
        !OperatorEx.isReadChannel(DataflowBroadcast.class)
    }

    def 'should check extension method' () {
        given:
        def ext = new OperatorEx()
        expect:
        ext.isExtensionMethod(new DataflowVariable(), 'map')
        ext.isExtensionMethod(new DataflowVariable(), 'flatMap')
        !ext.isExtensionMethod(new DataflowVariable(), 'foo')
    }

    def 'should invoke ext method' () {
        given:
        def ext = new OperatorEx()
        def ch = new DataflowQueue(); ch<<1<<2<<3

        when:
        def result = ext.invokeExtensionMethod(ch, 'map', { it -> it * it })
        then:
        result instanceof DataflowReadChannel
        result.val == 1
        result.val == 4
        result.val == 9 
    }

    def testFilter() {

        when:
        def c1 = Channel.from(1,2,3,4,5).filter { it > 3 }
        then:
        c1.val == 4
        c1.val == 5
        c1.val == Channel.STOP

        when:
        def c2 = Channel.from('hola','hello','cioa','miao').filter { it =~ /^h.*/ }
        then:
        c2.val == 'hola'
        c2.val == 'hello'
        c2.val == Channel.STOP

        when:
        def c3 = Channel.from('hola','hello','cioa','miao').filter { it ==~ /^h.*/ }
        then:
        c3.val == 'hola'
        c3.val == 'hello'
        c3.val == Channel.STOP

        when:
        def c4 = Channel.from('hola','hello','cioa','miao').filter( ~/^h.*/ )
        then:
        c4.val == 'hola'
        c4.val == 'hello'
        c4.val == Channel.STOP

        when:
        def c5 = Channel.from('hola',1,'cioa',2,3).filter( Number )
        then:
        c5.val == 1
        c5.val == 2
        c5.val == 3
        c5.val == Channel.STOP

        expect:
        Channel.from(1,2,4,2,4,5,6,7,4).filter(1) .count().val == 1
        Channel.from(1,2,4,2,4,5,6,7,4).filter(2) .count().val == 2
        Channel.from(1,2,4,2,4,5,6,7,4).filter(4) .count().val == 3

    }

    def testFilterWithValue() {
        expect:
        Channel.value(3).filter { it>1 }.val == 3
        Channel.value(0).filter { it>1 }.val == Channel.STOP
        Channel.value(Channel.STOP).filter { it>1 }.val == Channel.STOP
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
        channel = Channel.from(1,2,3,4)
        channel.subscribe { count++; }
        sleep(100)
        then:
        count == 4


    }

    def testSubscribe1() {

        when:
        def count = 0
        def done = false
        Channel.from(1,2,3).subscribe onNext:  { count++ }, onComplete: { done = true }
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
        def result = Channel.from(1,2,3).map { "Hello $it" }
        then:
        result.val == 'Hello 1'
        result.val == 'Hello 2'
        result.val == 'Hello 3'
        result.val == Channel.STOP
    }

    def testMapWithVariable() {
        given:
        def variable = Channel.value('Hello')
        when:
        def result = variable.map { it.reverse() }
        then:
        result.val == 'olleH'
        result.val == 'olleH'
        result.val == 'olleH'
    }

    def testMapParamExpanding () {

        when:
        def result = Channel.from(1,2,3).map { [it, it] }.map { x, y -> x+y }
        then:
        result.val == 2
        result.val == 4
        result.val == 6
        result.val == Channel.STOP
    }

    def testSkip() {

        when:
        def result = Channel.from(1,2,3).map { it == 2 ? Channel.VOID : "Hello $it" }
        then:
        result.val == 'Hello 1'
        result.val == 'Hello 3'
        result.val == Channel.STOP

    }


    def testMapMany () {

        when:
        def result = Channel.from(1,2,3).flatMap { it -> [it, it*2] }
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
        def result = Channel.value([1,2,3]).flatMap()
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = Channel.value().close().flatMap()
        then:
        result.val == Channel.STOP

    }

    def testMapManyWithTuples () {

        when:
        def result = Channel.from( [1,2], ['a','b'] ).flatMap { it -> [it, it.reverse()] }
        then:
        result.val == [1,2]
        result.val == [2,1]
        result.val == ['a','b']
        result.val == ['b','a']
        result.val == Channel.STOP
    }

    def testMapManyDefault  () {

        when:
        def result = Channel.from( [1,2], ['a',['b','c']] ).flatMap()
        then:
        result.val == 1
        result.val == 2
        result.val == 'a'
        result.val == ['b','c']  // <-- nested list are preserved
        result.val == Channel.STOP
    }

    def testMapManyWithHashArray () {

        when:
        def result = Channel.from(1,2,3).flatMap { it -> [ k: it, v: it*2] }
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
        def channel = Channel.create()
        def result = channel.reduce { a, e -> a += e }
        channel << 1 << 2 << 3 << 4 << 5 << Channel.STOP
        then:
        result.getVal() == 15


        when:
        channel = Channel.from(1,2,3,4,5)
        result = channel.reduce { a, e -> a += e }
        then:
        result.getVal() == 15

        when:
        channel = Channel.create()
        result = channel.reduce { a, e -> a += e }
        channel << 99 << Channel.STOP
        then:
        result.getVal() == 99

        when:
        channel = Channel.create()
        result = channel.reduce { a, e -> a += e }
        channel << Channel.STOP
        then:
        result.getVal() == null

        when:
        result = Channel.from(6,5,4,3,2,1).reduce { a, e -> Channel.STOP }
        then:
        result.val == 6

    }


    def testReduceWithSeed() {

        when:
        def channel = Channel.create()
        def result = channel.reduce (1) { a, e -> a += e }
        channel << 1 << 2 << 3 << 4 << 5 << Channel.STOP
        then:
        result.getVal() == 16

        when:
        channel = Channel.create()
        result = channel.reduce (10) { a, e -> a += e }
        channel << Channel.STOP
        then:
        result.getVal() == 10

        when:
        result = Channel.from(6,5,4,3,2,1).reduce(0) { a, e -> a < 3 ? a+1 : Channel.STOP }
        then:
        result.val == 3

    }

    def testFirst() {

        expect:
        Channel.from(3,6,4,5,4,3,4).first().val == 3
    }

    def testFirstWithCriteria() {
        expect:
        Channel.from(3,6,4,5,4,3,4).first{ it>4 } .val == 6
    }

    def testFirstWithValue() {

        expect:
        Channel.value(3).first().val == 3
        Channel.value(3).first{ it>1 }.val == 3
        Channel.value(3).first{ it>3 }.val == Channel.STOP
        Channel.value(Channel.STOP).first { it>3 }.val == Channel.STOP
    }


    def testFirstWithCondition() {

        expect:
        Channel.from(3,6,4,5,4,3,4).first { it % 2 == 0  } .val == 6
        Channel.from( 'a', 'b', 'c', 1, 2 ).first( Number ) .val == 1
        Channel.from( 'a', 'b', 1, 2, 'aaa', 'bbb' ).first( ~/aa.*/ ) .val == 'aaa'
        Channel.from( 'a', 'b', 1, 2, 'aaa', 'bbb' ).first( 1 ) .val == 1

    }


    def testTake() {

        when:
        def result = Channel.from(1,2,3,4,5,6).take(3)
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = Channel.from(1).take(3)
        then:
        result.val == 1
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2,3).take(0)
        then:
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2,3).take(-1)
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2,3).take(3)
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

    }

    def testLast() {

        expect:
        Channel.from(3,6,4,5,4,3,9).last().val == 9
        Channel.value('x').last().val == 'x'
    }




    def testCount() {
        expect:
        Channel.from(4,1,7,5).count().val == 4
        Channel.from(4,1,7,1,1).count(1).val == 3
        Channel.from('a','c','c','q','b').count ( ~/c/ ) .val == 2
        Channel.value(5).count().val == 1
        Channel.value(5).count(5).val == 1
        Channel.value(5).count(6).val == 0
    }

    def testCountBy() {
        expect:
        Channel.from('hello','ciao','hola', 'hi', 'bonjour').countBy { it[0] } .val == [c:1, b:1, h:3]
    }


    def testGroupBy() {

        def result

        when:
        result = Channel.from('hello','ciao','hola', 'hi', 'bonjour').groupBy { String str -> str[0] }
        then:
        result.val == [c:['ciao'], b:['bonjour'], h:['hello','hola','hi']]

        when:
        result = Channel.from( [id: 1, str:'one'], [id: 2, str:'two'], [id: 2, str:'dos'] ).groupBy()
        then:
        result.val == [ 1: [[id: 1, str:'one']], 2: [[id: 2, str:'two'], [id: 2, str:'dos']] ]

        when:
        result = Channel.from( [1, 'a' ], [2, 'b'], [1, 'c'], [1, 'd'], [3, 'z'] ).groupBy()
        then:
        result.val == [1: [[1,'a'], [1,'c'], [1,'d']], 2: [[2,'b']], 3: [[3,'z']]]
        result instanceof DataflowVariable

        when:
        result = Channel.from( [1, 'a' ], [2, 'b'], [3, 'a'], [4, 'c'], [5, 'a'] ).groupBy(1)
        then:
        result.val == [a: [[1,'a'], [3,'a'], [5,'a']], b: [[2,'b']], c: [[4,'c']]]
        result instanceof DataflowVariable

        when:
        result = Channel.from( [id: 1, str:'one'], [id: 2, str:'two'], [id: 2, str:'dos'] ).groupBy()
        then:
        result.val == [ 1: [[id: 1, str:'one']], 2: [[id: 2, str:'two'], [id: 2, str:'dos']] ]
        result instanceof DataflowVariable

        when:
        result = Channel.from('hello','ciao','hola', 'hi', 'bonjour').groupBy { String str -> str[0] }
        then:
        result.val == [c:['ciao'], b:['bonjour'], h:['hello','hola','hi']]
        result instanceof DataflowVariable

    }

    def testToList() {

        when:
        def channel = Channel.from(1,2,3)
        then:
        channel.toList().val == [1,2,3]

        when:
        channel = Channel.create()
        channel << Channel.STOP
        then:
        channel.toList().val == []

        when:
        channel = Channel.value(1)
        then:
        channel.toList().val == [1]

        when:
        channel = Channel.value().close()
        then:
        channel.toList().val == []
    }

    def testToSortedList() {

        when:
        def channel = Channel.from(3,1,4,2)
        then:
        channel.toSortedList().val == [1,2,3,4]

        when:
        channel = Channel.create()
        channel << Channel.STOP
        then:
        channel.toSortedList().val == []

        when:
        channel = Channel.from([1,'zeta'], [2,'gamma'], [3,'alpaha'], [4,'delta'])
        then:
        channel.toSortedList { it[1] } .val == [[3,'alpaha'], [4,'delta'], [2,'gamma'], [1,'zeta'] ]

        when:
        channel = Channel.value(1)
        then:
        channel.toSortedList().val == [1]

        when:
        channel = Channel.value().close()
        then:
        channel.toSortedList().val == []

    }


    def testUnique() {
        expect:
        Channel.from(1,1,1,5,7,7,7,3,3).unique().toList().val == [1,5,7,3]
        Channel.from(1,3,4,5).unique { it%2 } .toList().val == [1,4]
    }

    def testDistinct() {
        expect:
        Channel.from(1,1,2,2,2,3,1,1,2,2,3).distinct().toList().val == [1,2,3,1,2,3]
        Channel.from(1,1,2,2,2,3,1,1,2,4,6).distinct { it%2 } .toList().val == [1,2,3,2]
    }



    def testSpread() {

        when:
        def left = Channel.from(1,2,3)
        def right = ['aa','bb']
        def r1 = left.spread(right)
        then:
        r1.val == [1, 'aa']
        r1.val == [1, 'bb']
        r1.val == [2, 'aa']
        r1.val == [2, 'bb']
        r1.val == [3, 'aa']
        r1.val == [3, 'bb']
        r1.val == Channel.STOP

        when:
        left = Channel.from(1,2)
        right = Channel.from('a','bb','ccc')
        def r2 = left.spread(right)
        then:
        r2.val == [1, 'a']
        r2.val == [1, 'bb']
        r2.val == [1, 'ccc']
        r2.val == [2, 'a']
        r2.val == [2, 'bb']
        r2.val == [2, 'ccc']
        r2.val == Channel.STOP

    }

    def testSpreadChained() {

        when:
        def str1 = Channel.from('a','b','c')
        def str2 = Channel.from('x','y')
        def result = Channel.from(1,2).spread(str1).spread(str2)
        then:
        result.val == [1,'a','x']
        result.val == [1,'a','y']
        result.val == [1,'b','x']
        result.val == [1,'b','y']
        result.val == [1,'c','x']
        result.val == [1,'c','y']
        result.val == [2,'a','x']
        result.val == [2,'a','y']
        result.val == [2,'b','x']
        result.val == [2,'b','y']
        result.val == [2,'c','x']
        result.val == [2,'c','y']
        result.val == Channel.STOP

    }


    def testSpreadTuple() {

        when:
        def result = Channel.from([1, 'x'], [2,'y'], [3, 'z']).spread( ['alpha','beta','gamma'] )

        then:
        result.val == [1, 'x', 'alpha']
        result.val == [1, 'x', 'beta']
        result.val == [1, 'x', 'gamma']

        result.val == [2, 'y', 'alpha']
        result.val == [2, 'y', 'beta']
        result.val == [2, 'y', 'gamma']

        result.val == [3, 'z', 'alpha']
        result.val == [3, 'z', 'beta']
        result.val == [3, 'z', 'gamma']

        result.val == Channel.STOP
    }

    def testSpreadMap() {

        when:
        def result = Channel.from([id:1, val:'x'], [id:2,val:'y'], [id:3, val:'z']).spread( ['alpha','beta','gamma'] )

        then:
        result.val == [[id:1, val:'x'], 'alpha']
        result.val == [[id:1, val:'x'], 'beta']
        result.val == [[id:1, val:'x'], 'gamma']

        result.val == [[id:2,val:'y'], 'alpha']
        result.val == [[id:2,val:'y'], 'beta']
        result.val == [[id:2,val:'y'], 'gamma']

        result.val == [[id:3, val:'z'], 'alpha']
        result.val == [[id:3, val:'z'], 'beta']
        result.val == [[id:3, val:'z'], 'gamma']

        result.val == Channel.STOP
    }

    def testSpreadWithSingleton() {
        when:
        def result = Channel.value(7).spread(['a','b','c'])
        then:
        result.val == [7, 'a']
        result.val == [7, 'b']
        result.val == [7, 'c']
        result.val == Channel.STOP
    }

    def testSpreadWithPath() {
        given:
        def path1 = Paths.get('/some/data/file1')
        def path2 = Paths.get('/some/data/file2')

        when:
        def result = Channel.from(1,2,3).spread( Channel.value(path1) )
        then:
        result.val == [1, path1]
        result.val == [2, path1]
        result.val == [3, path1]
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2,3).spread( ['abc'] )
        then:
        result.val == [1, 'abc']
        result.val == [2, 'abc']
        result.val == [3, 'abc']
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2,3).spread( Channel.value('abc') )
        then:
        result.val == [1, 'abc']
        result.val == [2, 'abc']
        result.val == [3, 'abc']
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2).spread( Channel.from(path1,path2) )
        then:
        result.val == [1, path1]
        result.val == [1, path2]
        result.val == [2, path1]
        result.val == [2, path2]
        result.val == Channel.STOP

    }

    def testFlatten() {

        when:
        def r1 = Channel.from(1,2,3).flatten()
        then:
        r1.val == 1
        r1.val == 2
        r1.val == 3
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from([1,'a'], [2,'b']).flatten()
        then:
        r2.val == 1
        r2.val == 'a'
        r2.val == 2
        r2.val == 'b'
        r2.val == Channel.STOP

        when:
        def r3 = Channel.from( [1,2] as Integer[], [3,4] as Integer[] ).flatten()
        then:
        r3.val == 1
        r3.val == 2
        r3.val == 3
        r3.val == 4
        r3.val == Channel.STOP

        when:
        def r4 = Channel.from( [1,[2,3]], 4, [5,[6]] ).flatten()
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
        def result = Channel.value([3,2,1]).flatten()
        then:
        result.val == 3
        result.val == 2
        result.val == 1
        result.val == Channel.STOP

        when:
        result = Channel.value().close().flatten()
        then:
        result.val ==  Channel.STOP
    }

    def testCollate() {

        when:
        def r1 = Channel.from(1,2,3,1,2,3,1).collate( 2, false )
        then:
        r1.val == [1,2]
        r1.val == [3,1]
        r1.val == [2,3]
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from(1,2,3,1,2,3,1).collate( 3 )
        then:
        r2.val == [1,2,3]
        r2.val == [1,2,3]
        r2.val == [1]
        r2.val == Channel.STOP

    }

    def testCollateWithStep() {

        when:
        def r1 = Channel.from(1,2,3,4).collate( 3, 1, false )
        then:
        r1.val == [1,2,3]
        r1.val == [2,3,4]
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from(1,2,3,4).collate( 3, 1, true )
        then:
        r2.val == [1,2,3]
        r2.val == [2,3,4]
        r2.val == [3,4]
        r2.val == [4]
        r2.val == Channel.STOP

        when:
        def r3 = Channel.from(1,2,3,4).collate( 3, 1  )
        then:
        r3.val == [1,2,3]
        r3.val == [2,3,4]
        r3.val == [3,4]
        r3.val == [4]
        r3.val == Channel.STOP

        when:
        def r4 = Channel.from(1,2,3,4).collate( 4,4 )
        then:
        r4.val == [1,2,3,4]
        r4.val == Channel.STOP

        when:
        def r5 = Channel.from(1,2,3,4).collate( 6,6 )
        then:
        r5.val == [1,2,3,4]
        r5.val == Channel.STOP

        when:
        def r6 = Channel.from(1,2,3,4).collate( 6,6,false )
        then:
        r6.val == Channel.STOP

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
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = Channel.value(1).collate(10)
        then:
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = Channel.value(1).collate(10, true)
        then:
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = Channel.value(1).collate(10, false)
        then:
        result.val == Channel.STOP
    }

    def testMix() {
        when:
        def c1 = Channel.from( 1,2,3 )
        def c2 = Channel.from( 'a','b' )
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
        def result = Channel.value(1).mix( Channel.from([2,3])  )
        then:
        result.toList().val.sort() == [1,2,3]
    }



    def testDefaultMappingClosure() {

        expect:
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9] ) == 7
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9], 2 ) == 9
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [] ) == null
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [], 2 ) == null

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9] as Object[] ) == 7
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [7, 8, 9] as Object[], 1 ) == 8
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( ['alpha', 'beta'] as String[] ) == 'alpha'
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( ['alpha', 'beta'] as String[], 1 ) == 'beta'

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet ) == 6
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet, 1 ) == 7
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet, 2 ) == 8
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [6, 7, 8, 9 ] as LinkedHashSet, 5 ) == null

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9] ) == 1
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9], 1 ) == 2
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9], 2 ) == 9
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9], 3 ) == null

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(0) ) == 'a'
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(0), 1 ) == 1
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(0), 2 ) == null

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(1) ) == 'b'
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(1), 1 ) == 2
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9].entrySet().getAt(1), 2 ) == null

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( [:] ) == null

        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( 99 ) == 99
        OperatorEx.DEFAULT_MAPPING_CLOSURE.call( 99, 2 ) == null

    }


    def testCross() {

        setup:
        def ch1 = Channel.from(  [1, 'x'], [2,'y'], [3,'z'] )
        def ch2 = Channel.from( [1,11], [1,13], [2,21],[2,22], [2,23], [4,1], [4,2]  )

        when:
        def result = ch1.cross(ch2)

        then:
        result.val == [ [1, 'x'], [1,11] ]
        result.val == [ [1, 'x'], [1,13] ]
        result.val == [ [2, 'y'], [2,21] ]
        result.val == [ [2, 'y'], [2,22] ]
        result.val == [ [2, 'y'], [2,23] ]
        result.val == Channel.STOP

    }

    def testCross2() {

        setup:
        def ch1 = Channel.create()
        def ch2 = Channel.from ( ['PF00006', 'PF00006_mafft.aln'], ['PF00006', 'PF00006_clustalo.aln'])

        when:
        Thread.start {  sleep 100;   ch1 << ['PF00006', 'PF00006.sp_lib'] << Channel.STOP }
        def result = ch1.cross(ch2)

        then:
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_mafft.aln'] ]
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_clustalo.aln'] ]
        result.val == Channel.STOP

    }


    def testCross3() {

        setup:
        def ch1 = Channel.from([['PF00006', 'PF00006.sp_lib'] ])
        def ch2 = Channel.create ( )

        when:
        Thread.start {  sleep 100;  ch2 << ['PF00006', 'PF00006_mafft.aln'] <<  ['PF00006', 'PF00006_clustalo.aln']<< Channel.STOP }
        def result = ch1.cross(ch2)

        then:
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_mafft.aln'] ]
        result.val == [ ['PF00006', 'PF00006.sp_lib'], ['PF00006', 'PF00006_clustalo.aln'] ]
        result.val == Channel.STOP

    }


    def testConcat() {

        when:
        def c1 = Channel.from(1,2,3)
        def c2 = Channel.from('a','b','c')
        def all = c1.concat(c2)
        then:
        all.val == 1
        all.val == 2
        all.val == 3
        all.val == 'a'
        all.val == 'b'
        all.val == 'c'
        all.val == Channel.STOP

        when:
        def d1 = Channel.create()
        def d2 = Channel.from('a','b','c')
        def d3 = Channel.create()
        def result = d1.concat(d2,d3)

        Thread.start { sleep 20; d3 << 'p' << 'q' << Channel.STOP }
        Thread.start { sleep 100; d1 << 1 << 2 << Channel.STOP }

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
        def result = Channel.value(1).concat( Channel.from(2,3) )
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
    }


    def testDataflowSeparateWithOpenArray() {

        when:
        def s1 = Channel.create()
        def s2 = Channel.create()
        def s3 = Channel.create()

        Channel.from(1,2,3,4)
                .separate(s1,s2,s3) { item -> [item+1, item*item, item-1] }

        then:
        s1.val == 2
        s1.val == 3
        s1.val == 4
        s1.val == 5
        s1.val == Channel.STOP
        s2.val == 1
        s2.val == 4
        s2.val == 9
        s2.val == 16
        s2.val == Channel.STOP
        s3.val == 0
        s3.val == 1
        s3.val == 2
        s3.val == 3
        s3.val == Channel.STOP

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

    def testGroupTupleWithClosure() {

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

    def testChannelIfEmpty() {

        def result

        when:
        result = Channel.from(1,2,3).ifEmpty(100)
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
        result = Channel.value().close().ifEmpty(100)
        then:
        result instanceof DataflowVariable
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
        Channel.from(10,20,30)
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
        def result = Channel.from(1,2,3,4).until { it == 3 }
        then:
        result.val == 1
        result.val == 2
        result.val == Channel.STOP

        when:
        result = Channel.from(1,2,3).until { it == 5 }
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
        Channel.from(1,2,3).set { result }

        then:
        session.binding.result.val == 1
        session.binding.result.val == 2
        session.binding.result.val == 3
        session.binding.result.val == Channel.STOP
    }

}
