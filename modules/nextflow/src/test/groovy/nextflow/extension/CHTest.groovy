package nextflow.extension

import spock.lang.Specification

import java.util.concurrent.TimeUnit

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.NextflowMeta
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CHTest extends Specification {

    def 'should create dataflow variable or queue' () {

        expect:
        CH.create() instanceof DataflowQueue
        CH.create(false) instanceof DataflowQueue
        CH.create(true) instanceof DataflowVariable

        CH.createBy(new DataflowVariable()) instanceof DataflowVariable
        CH.createBy(new DataflowQueue()) instanceof DataflowQueue


        when:
        NextflowMeta.instance.enableDsl2()
        then:
        CH.create() instanceof DataflowBroadcast
        CH.create(false) instanceof DataflowBroadcast
        CH.create(true) instanceof DataflowVariable

        CH.createBy(new DataflowVariable()) instanceof DataflowVariable
        CH.createBy(new DataflowQueue()) instanceof DataflowBroadcast

        cleanup:
        NextflowMeta.instance.disableDsl2()

    }


    def 'should check queue channel' () {
        expect:
        CH.isChannelQueue(new DataflowQueue())
        CH.isChannelQueue(new DataflowBroadcast().createReadChannel())
        !CH.isChannelQueue(new DataflowVariable())
        !CH.isChannelQueue('hello')
    }

    def 'should validate allScalar method' () {

        expect:
        CH.allScalar([1])
        CH.allScalar([1, 2, 3])
        !CH.allScalar([new DataflowVariable(), new DataflowQueue()])
        !CH.allScalar([new DataflowQueue(), new DataflowQueue()])
    }

    def 'should check value' () {
        when:
        def v1 = CH.value()
        then:
        v1 instanceof DataflowVariable
        !v1.isBound()

        when:
        def v2 = CH.value('hello')
        then:
        v2 instanceof DataflowVariable
        v2.val == 'hello'
    }

    def 'should bind a value' () {
        given:
        def ch = new DataflowQueue()
        when:
        CH.emit(ch, 'hello')
        then:
        ch.val == 'hello'
    }

    def 'should emit a list of values' () {
        given:
        def ch = new DataflowQueue()
        when:
        CH.emitValues(ch, [1, 2, 3])
        then:
        ch.val == 1
        ch.val == 2
        ch.val == 3
        and:
        // no stop is emitted
        ch.getVal(500, TimeUnit.MILLISECONDS) == null
    }


    def 'should create dataflow empty queue' () {
        when:
        def ch = CH.queue()
        then:
        ch instanceof DataflowQueue
        and:
        ch.getVal(200, TimeUnit.MILLISECONDS) == null
    }

    def 'should create dataflow queue w/o stop' () {
        when:
        def ch = CH.emitValues(CH.queue(), [1,2,3])
        then:
        ch instanceof DataflowQueue
        ch.val == 1
        ch.val == 2
        ch.val == 3
        and:
        ch.getVal(200, TimeUnit.MILLISECONDS) == null
    }

    def 'should create dataflow queue with stop' () {
        when:
        def ch = CH.emitAndClose(CH.queue(), [1,2,3])
        then:
        ch instanceof DataflowQueue
        ch.val == 1
        ch.val == 2
        ch.val == 3
        and:
        ch.val == Channel.STOP
    }

}
