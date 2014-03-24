/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import spock.lang.Specification


/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowExtensionsTest extends Specification {

    def testHandlerNames() {

        when:
        DataflowExtensions.checkSubscribeHandlers( [:] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowExtensions.checkSubscribeHandlers( [ onNext:{}] )
        then:
        true

        when:
        DataflowExtensions.checkSubscribeHandlers( [ onNext:{}, xxx:{}] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowExtensions.checkSubscribeHandlers( [ xxx:{}] )
        then:
        thrown(IllegalArgumentException)


    }

    def testGrep() {

        when:
        def channel = new DataflowQueue()
        channel << 'Hello' << 'Halo' << 'Ciao' << 'Hola' << 'Bonjour'
        def result = channel.grep( ~/^H.*/ )
        then:
        result.val == 'Hello'
        result.val == 'Halo'
        result.val == 'Hola'


        when:
        def channel2 = new DataflowQueue()
        channel2 << 'Hello' << 1 << 'Ciao' << 2 << 'Bonjour'
        def result2 = channel2.grep(Number)
        then:
        result2.val == 1
        result2.val == 2

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

//
//        when:
//        int a=0
//        int b=0
//        def channel3 = Channel.from(1,2,3,4,5)
//        channel3.subscribe { if(it%2 != 0) a++ }
//        channel3.subscribe { if(it%2 == 0) b++; println ">> $it"  }
//        //channel3 << 1 << 2 << 3 << 4 << 5
//        sleep(100)
//        then:
//        a == 3
//        b == 2


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
        Channel.just(1).subscribe onNext:  { count++ }, onComplete: { done = true }
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


    def testDoFinally() {

        when:
        def flag
        Channel.from(1,2,3).doFinally { flag = 1 }
        sleep 100
        then:
        flag == 1

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

    }

    def testLast() {

        expect:
        Channel.from(3,6,4,5,4,3,9).last().val == 9
        Channel.just('x').last().val == 'x'
    }

    def testSplit() {
        when:
        def result = Channel.from(1,2,3,4)
        def (ch1, ch2) = result.split(2)

        then:
        ch1.val == 1
        ch1.val == 2
        ch1.val == 3
        ch1.val == 4
        ch1.val == Channel.STOP

        ch2.val == 1
        ch2.val == 2
        ch2.val == 3
        ch2.val == 4
        ch2.val == Channel.STOP

    }



    def testMin() {

        expect:
        Channel.from(4,1,7,5).min().val == 1
        Channel.from("hello","hi","hey").min { it.size() } .val == "hi"
        Channel.from("hello","hi","hey").min { a,b -> a.size()<=>b.size() } .val == "hi"
        Channel.from("hello","hi","hey").min { a,b -> a.size()<=>b.size() } .val == "hi"
        Channel.from("hello","hi","hey").min ({ a,b -> a.size()<=>b.size() } as Comparator) .val == "hi"

    }

    def testMax() {
        expect:
        Channel.from(4,1,7,5).max().val == 7
        Channel.from("hello","hi","hey").max { it.size() } .val == "hello"
        Channel.from("hello","hi","hey").max { a,b -> a.size()<=>b.size() } .val == "hello"
        Channel.from("hello","hi","hey").max { a,b -> a.size()<=>b.size() } .val == "hello"
        Channel.from("hello","hi","hey").max ({ a,b -> a.size()<=>b.size() } as Comparator) .val == "hello"

    }

    def testSum() {
        expect:
        Channel.from(4,1,7,5).sum().val == 17
        Channel.from(4,1,7,5).sum { it * 2 } .val == 34
    }

    def testCount() {
        expect:
        Channel.from(4,1,7,5).count().val == 4
        Channel.from(4,1,7,1,1).count(1).val == 3
        Channel.from('a','c','c','q','b').count ( ~/c/ ) .val == 2
    }

    def testCountBy() {
        expect:
        Channel.from('hello','ciao','hola', 'hi', 'bonjour').countBy { it[0] } .val == [c:1, b:1, h:3]
    }

    def testGroupBy() {

        expect:
        Channel.from('hello','ciao','hola', 'hi', 'bonjour').groupBy { String str -> str[0] } .val == [c:['ciao'], b:['bonjour'], h:['hello','hola','hi']]

        Channel.from( [id: 1, str:'one'], [id: 2, str:'two'], [id: 2, str:'dos'] ).groupBy().val == [ 1: [[id: 1, str:'one']], 2: [[id: 2, str:'two'], [id: 2, str:'dos']] ]


    }

    def testRouteBy() {

        when:
        def result = Channel.from('hello','ciao','hola', 'hi', 'bonjour').route { it[0] }
        result.subscribe {
          def (key, channel) = it
          channel.subscribe { println "group: $key >> $it" }
        }

        sleep 1000
        then:
        true

    }

    def testRouteByMap() {

        setup:
        def r1 = Channel.create()
        def r2 = Channel.create()
        def r3 = Channel.create()

        when:
        Channel.from('hello','ciao','hola', 'hi', 'x', 'bonjour').route ( b: r1, c: r2, h: r3 ) { it[0] }

        then:
        r1.val == 'bonjour'
        r1.val == Channel.STOP

        r2.val == 'ciao'
        r2.val == Channel.STOP

        r3.val == 'hello'
        r3.val == 'hola'
        r3.val == 'hi'
        r3.val == Channel.STOP

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


    def testSeparate() {

        when:
        def str = 'abcdef'
        def (ch1, ch2) = Channel.from(0..3).separate(2) { [it, str[it]] }
        then:
        ch1.val == 0
        ch1.val == 1
        ch1.val == 2
        ch1.val == 3
        ch1.val == Channel.STOP

        ch2.val == 'a'
        ch2.val == 'b'
        ch2.val == 'c'
        ch2.val == 'd'
        ch2.val == Channel.STOP

    }

    def testSeparate2() {

        when:
        def s1 = Channel.create()
        def s2 = Channel.create()
        Channel.from(1,2,3,4)
                .map { it }
                .separate(s1,s2) { item -> [item+1, item*item] }

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





    }

    def testSpread() {

        when:
        def r1 = Channel.from(1,2,3).spread(['a','b'])
        then:
        r1.val == [1, 'a']
        r1.val == [1, 'b']
        r1.val == [2, 'a']
        r1.val == [2, 'b']
        r1.val == [3, 'a']
        r1.val == [3, 'b']
        r1.val == Channel.STOP

        when:
        def str = Channel.from('a','b','c')
        def r2 = Channel.from(1,2).spread(str)
        then:
        r2.val == [1, 'a']
        r2.val == [1, 'b']
        r2.val == [1, 'c']
        r2.val == [2, 'a']
        r2.val == [2, 'b']
        r2.val == [2, 'c']
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

        when:
        def r5 = Channel.just( [1,2,3] ).flatten()
        then:
        r5.val == 1
        r5.val == 2
        r5.val == 3
        r5.val == Channel.STOP

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





    def testBufferClose() {

        when:
        def r1 = Channel.from(1,2,3,1,2,3).buffer({ it == 2 })
        then:
        r1.val == [1,2]
        r1.val == [3,1,2]
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from('a','b','c','a','b','z').buffer(~/b/)
        then:
        r2.val == ['a','b']
        r2.val == ['c','a','b']
        r2.val == Channel.STOP

    }

    def testBufferWithCount() {

        when:
        def r1 = Channel.from(1,2,3,1,2,3,1).buffer( size:2 )
        then:
        r1.val == [1,2]
        r1.val == [3,1]
        r1.val == [2,3]
        r1.val == Channel.STOP

        when:
        r1 = Channel.from(1,2,3,1,2,3,1).buffer( size:2, remainder: true )
        then:
        r1.val == [1,2]
        r1.val == [3,1]
        r1.val == [2,3]
        r1.val == [1]
        r1.val == Channel.STOP


        when:
        def r2 = Channel.from(1,2,3,4,5,1,2,3,4,5,1,2,9).buffer( size:3, skip:2 )
        then:
        r2.val == [3,4,5]
        r2.val == [3,4,5]
        r2.val == Channel.STOP

        when:
        r2 = Channel.from(1,2,3,4,5,1,2,3,4,5,1,2,9).buffer( size:3, skip:2, remainder: true )
        then:
        r2.val == [3,4,5]
        r2.val == [3,4,5]
        r2.val == [9]
        r2.val == Channel.STOP

    }

    def testBufferInvalidArg() {

        when:
        Channel.create().buffer( xxx: true )

        then:
        IllegalArgumentException e = thrown()
        e.message == "Unknown argument argument 'xxx' for operator 'buffer' -- Possible arguments: size, skip, remainder"


    }


    def testBufferOpenClose() {

        when:
        def r1 = Channel.from(1,2,3,4,5,1,2,3,4,5,1,2).buffer( 2, 4 )
        then:
        r1.val == [2,3,4]
        r1.val == [2,3,4]
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from('a','b','c','a','b','z').buffer(~/a/,~/b/)
        then:
        r2.val == ['a','b']
        r2.val == ['a','b']
        r2.val == Channel.STOP

    }

    def testMix() {
        when:
        def c1 = Channel.from( 1,2,3 )
        def c2 = Channel.from( 'a','b' )
        def c3 = Channel.just( 'z' )
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


    def testPhaseImpl() {

        setup:
        def result = null
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()

        when:
        def map = [ : ]
        result = DataflowExtensions.phaseImpl(map, 2, 0, 'a', { it })
        then:
        result == null
        map == [ a:[0: ['a']] ]

        when:
        map = [ : ]
        result = DataflowExtensions.phaseImpl(map, 2, 0, 'a', { it })
        result = DataflowExtensions.phaseImpl(map, 2, 1, 'a', { it })
        then:
        result == ['a','a']
        map == [ a:[:] ]


        when:
        def r1
        def r2
        def r3
        map = [ : ]
        r1 = DataflowExtensions.phaseImpl(map, 3, 0, 'a', { it })
        r1 = DataflowExtensions.phaseImpl(map, 3, 1, 'a', { it })
        r1 = DataflowExtensions.phaseImpl(map, 3, 2, 'a', { it })

        r2 = DataflowExtensions.phaseImpl(map, 3, 0, 'b', { it })
        r2 = DataflowExtensions.phaseImpl(map, 3, 1, 'b', { it })
        r2 = DataflowExtensions.phaseImpl(map, 3, 2, 'b', { it })

        r3 = DataflowExtensions.phaseImpl(map, 3, 0, 'z', { it })
        r3 = DataflowExtensions.phaseImpl(map, 3, 1, 'z', { it })
        r3 = DataflowExtensions.phaseImpl(map, 3, 1, 'z', { it })
        r3 = DataflowExtensions.phaseImpl(map, 3, 2, 'z', { it })

        then:
        r1 == ['a','a','a']
        r2 == ['b','b','b']
        r3 == ['z','z','z']
        map == [ a:[:], b:[:], z:[ 1:['z']] ]

    }

    def testDefaultMappingClosure() {

        expect:
        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( [a:1, b:2, z:9] ) == 1
        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( [:] ) == null

        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( [3,2,1] ) == 3
        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( [] ) == null

        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( [1,2,3] as Object[] ) == 1
        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( ['alpha','beta'] as String[] ) == 'alpha'
        DataflowExtensions.DEFAULT_MAPPING_CLOSURE.call( 99 ) == 99

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

    def testChopFasta() {

        setup:
        def fasta = """\
                >1aboA
                NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
                NYITPVN
                >1ycsB
                KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
                VPRNLLGLYP
                ; comment
                >1pht
                GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
                WLNGYNETTGERGDFPGTYVEYIGRKKISP
                """.stripIndent()

        when:
        def records = Channel.from(fasta).chopFasta()
        then:
        records.val == '>1aboA\nNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS\nNYITPVN\n'
        records.val == '>1ycsB\nKGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY\nVPRNLLGLYP\n'
        records.val == '>1pht\nGYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG\nWLNGYNETTGERGDFPGTYVEYIGRKKISP\n'
        records.val == Channel.STOP

        when:
        def fasta2 = '''
            >alpha123
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            VPRNLLGLYP
            '''
        records = Channel.from(fasta, fasta2).chopFasta(record:[id:true])
        then:
        records.val == [id:'1aboA']
        records.val == [id:'1ycsB']
        records.val == [id:'1pht']
        records.val == [id:'alpha123']
        records.val == Channel.STOP


        when:
        def ids = []
        def list = []
        Channel.from(fasta).chopFasta(record:[id:true]) { item, int index ->
            ids << item.id
            list << index
            return item
        }

        sleep 100
        then:
        ids == ['1aboA','1ycsB','1pht']
        list == [0,1,2]

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

    def testInto() {
        when:
        def x = Channel.create()
        def y = Channel.create()
        def source = Channel.from([1,2], ['a','b'], ['p','q'])
        source.into(x,y)
        then:
        x.val == 1
        x.val == 'a'
        x.val == 'p'
        x.val == Channel.STOP
        y.val == 2
        y.val == 'b'
        y.val == 'q'
        y.val == Channel.STOP

        when:
        def x2 = Channel.create()
        def y2 = Channel.create()
        def source2 = Channel.from([1,2], ['a','c','b'], 'z')
        source2.into(x2,y2)
        then:
        x2.val == 1
        x2.val == 'a'
        x2.val == 'z'
        x2.val == Channel.STOP
        y2.val == 2
        y2.val == 'c'
        y2.val == null
        y2.val == Channel.STOP
    }


}