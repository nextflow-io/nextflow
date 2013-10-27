package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import spock.lang.Specification


/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowExtensionsTest extends Specification {

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

    def testMapMany () {

        when:
        def result = Channel.from(1,2,3).mapMany { it -> [it, it*2] }
        then:
        result.val == 1
        result.val == 2
        result.val == 2
        result.val == 4
        result.val == 3
        result.val == 6
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
        Channel.from('a','c','c','q','b').count { it == 'c' } .val == 2
    }

    def testCountBy() {
        expect:
        Channel.from('hello','ciao','hola', 'hi', 'bonjour').countBy { it[0] } .val == [c:1, b:1, h:3]
    }

    def testGroupBy() {

        expect:
        Channel.from('hello','ciao','hola', 'hi', 'bonjour').groupBy { it[0] } .val == [c:['ciao'], b:['bonjour'], h:['hello','hola','hi']]

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

}