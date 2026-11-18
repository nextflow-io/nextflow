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

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Channel
import spock.lang.Specification
import spock.lang.Timeout

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CHTest extends Specification {

    static final PoisonPill STOP = PoisonPill.getInstance()

    def 'should create dataflow variable or queue' () {

        expect:
        CH.create() instanceof DataflowBroadcast
        CH.create(false) instanceof DataflowBroadcast
        CH.create(true) instanceof DataflowVariable

        CH.createBy(new DataflowVariable()) instanceof DataflowVariable
        CH.createBy(new DataflowQueue()) instanceof DataflowBroadcast
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
        def v2 = runDataflow {
            CH.value('hello')
        }
        then:
        v2 instanceof DataflowVariable
        v2.val == 'hello'
    }

    def 'should bind a value' () {
        given:
        def ch = new DataflowQueue()
        when:
        runDataflow {
            CH.emit(ch, 'hello')
        }
        then:
        ch.val == 'hello'
    }

    def 'should emit a list of values' () {
        given:
        def ch = new DataflowQueue()
        when:
        runDataflow {
            CH.emitValues(ch, [1, 2, 3])
        }
        then:
        ch.val == 1
        ch.val == 2
        ch.val == 3
        and:
        // no stop is emitted
        ch.length() == 0
    }


    def 'should create dataflow empty queue' () {
        when:
        def ch = CH.queue()
        then:
        ch instanceof DataflowQueue
        and:
        ch.length() == 0
    }

    def 'should create dataflow queue w/o stop' () {
        given:
        def ch = CH.queue()
        when:
        runDataflow {
            CH.emitValues(ch, [1,2,3])
        }
        then:
        ch instanceof DataflowQueue
        ch.val == 1
        ch.val == 2
        ch.val == 3
        and:
        ch.length() == 0
    }

    def 'should create dataflow queue with stop' () {
        given:
        def ch = CH.queue()
        when:
        runDataflow {
            CH.emitAndClose(ch, [1,2,3])
        }
        then:
        ch instanceof DataflowQueue
        ch.val == 1
        ch.val == 2
        ch.val == 3
        and:
        ch.val == Channel.STOP
    }

    static List read(DataflowReadChannel channel, int n) {
        final result = []
        for( int i = 0; i < n; i++ )
            result.add(channel.getVal())
        return result
    }

    @Timeout(10)
    def 'should forward every item and the poison pill of a queue' () {
        given:
        def queue = new DataflowQueue()
        def target = new DataflowBroadcast()
        def reader = target.createReadChannel()

        when:
        CH.forward(queue, target)
        for( def item : [1, 2, 3, STOP] )
            queue.bind(item)
        then:
        read(reader, 4) == [1, 2, 3, STOP]
    }

    @Timeout(10)
    def 'should forward every item of a broadcast reader without consuming the other readers' () {
        given:
        def broadcast = new DataflowBroadcast()
        def source = broadcast.createReadChannel()
        def other = broadcast.createReadChannel()
        def target = new DataflowQueue()

        when:
        CH.forward(source, target)
        for( def item : [1, 2, STOP] )
            broadcast.bind(item)
        then:
        read(target, 3) == [1, 2, STOP]
        read(other, 3) == [1, 2, STOP]
    }

    @Timeout(10)
    def 'should forward a value once bound' () {
        given:
        def value = new DataflowVariable()
        def target = new DataflowQueue()

        when:
        CH.forward(value, target)
        value.bind('x')
        then:
        target.getVal() == 'x'
    }

    @Timeout(10)
    def 'should forward a bound value in the caller thread' () {
        given:
        def value = new DataflowVariable()
        def target = new DataflowQueue()
        value.bind('x')

        when:
        CH.forward(value, target)
        then:
        target.length() == 1
    }

    @Timeout(10)
    def 'should forward the stop of a value and its end of channel' () {
        // a value bound to STOP, as by CH.close0 on an unbound value
        when:
        def result = [true, false].collect { bound ->
            def value = new DataflowVariable()
            def target = new DataflowVariable()
            if( bound )
                value.bind(STOP)
            CH.forward(value, target)
            if( !bound )
                CH.close0(value)
            return target.getVal()
        }
        then:
        result[0].is(STOP)
        result[1].is(STOP)
    }

    @Timeout(10)
    def 'should bind the result of the wrapped closure to every output' () {
        given:
        def input = new DataflowQueue()
        def o1 = new DataflowQueue()
        def o2 = new DataflowQueue()
        def triple = { Object x -> (x as Integer) * 3 }
        def copy = new DataflowQueue()

        when:
        new DataflowProcessor([input], [o1, o2], null, 1, new ChainWithClosure(triple)).start()
        new DataflowProcessor([o2], [copy], null, 1, new ChainWithClosure(new CopyChannelsClosure())).start()
        for( def item : [1, 2, STOP] )
            input.bind(item)
        then:
        read(o1, 3) == [3, 6, STOP]
        read(copy, 3) == [3, 6, STOP]
    }

}
