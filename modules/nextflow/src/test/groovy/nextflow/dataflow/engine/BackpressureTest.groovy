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
package nextflow.dataflow.engine

import java.util.concurrent.TimeUnit

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.SysEnv
import nextflow.extension.CH
import spock.lang.Specification
import spock.lang.Timeout

import static nextflow.dataflow.engine.BackpressureInternals.await

/**
 * Unit tests of the source gates.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(value = 20, unit = TimeUnit.SECONDS)
class BackpressureTest extends Specification {

    /**
     * A backpressure with a source channel feeding one process.
     */
    static class Gated {
        final Backpressure bp
        final DataflowBroadcast source
        final ProcessBacklog backlog

        Gated(Backpressure bp, DataflowBroadcast source, ProcessBacklog backlog) {
            this.bp = bp
            this.source = source
            this.backlog = backlog
        }

        Backpressure.Gate gate() {
            return bp.gate(source)
        }
    }

    private static Gated gated(int highWater) {
        final bp = new Backpressure(null, highWater)
        final source = new DataflowBroadcast()
        final backlog = new ProcessBacklog(bp, 'foo', [source.createReadChannel()], null)
        bp.markSource('channel.of', source)
        bp.gate(source).setFeeds([backlog])
        bp.addProcess(backlog)
        return new Gated(bp, source, backlog)
    }

    private static Thread sourceThread(Backpressure bp, Object source) {
        final thread = new Thread({ Backpressure.runAsSource({ bp.awaitOpen(source) } as Runnable) } as Runnable)
        thread.start()
        return thread
    }

    private static List range(int from, int to) {
        return (from..to).toList()
    }

    private static <T> T withEnv(Map<String,String> env, Closure<T> action) {
        SysEnv.push(env)
        try {
            return action.call()
        }
        finally {
            SysEnv.pop()
        }
    }

    def 'should get the high-water mark' () {
        expect:
        withEnv([:]) { Backpressure.highWater(4) } == 64
        withEnv([:]) { Backpressure.highWater(16) } == 64
        withEnv([:]) { Backpressure.highWater(50) } == 100

        and:
        withEnv(NXF_DATAFLOW_HIGH_WATER: '16') { Backpressure.highWater(50) } == 16
        withEnv(NXF_DATAFLOW_HIGH_WATER: '16') { new Backpressure(null, 16).getLowWater() } == 8

        and:
        // zero disables backpressure
        withEnv(NXF_DATAFLOW_HIGH_WATER: '0') { Backpressure.highWater(4) } == 0
    }

    def 'should reject an invalid high-water mark' () {
        given:
        SysEnv.push(NXF_DATAFLOW_HIGH_WATER: value)

        when:
        Backpressure.highWater(4)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Invalid NXF_DATAFLOW_HIGH_WATER value: '" + value + "'")

        cleanup:
        SysEnv.pop()

        where:
        value << ['-1', 'abc', '1.5', 'true']
    }

    def 'should never block a thread not flagged as a source thread' () {
        given:
        def ctx = gated(4)
        4.times { ctx.backlog.taskPending() }

        expect:
        !ctx.bp.isOpen(ctx.source)

        when:
        // returns immediately, otherwise the test times out
        ctx.bp.awaitOpen(ctx.source)
        Backpressure.awaitSource(ctx.source)
        then:
        noExceptionThrown()
    }

    def 'should block a source thread from the high-water mark down to the low-water mark' () {
        given:
        def ctx = gated(8)
        8.times { ctx.backlog.taskPending() }

        when:
        def thread = sourceThread(ctx.bp, ctx.source)
        await { ctx.gate().isWaiting() }
        then:
        ctx.gate().getBlockCount() == 1

        when:
        // the backlog drops below the high-water mark but stays above the low-water mark
        3.times { ctx.backlog.taskSubmitted() }
        // a short wait to catch an early wake-up: it can only miss a bug, never cause a false failure
        Thread.sleep(50)
        then:
        thread.isAlive()

        when:
        ctx.backlog.taskSubmitted()
        thread.join(5_000)
        then:
        !thread.isAlive()
        !ctx.gate().isWaiting()
    }

    def 'should count the threads waiting on a gate' () {
        given:
        def ctx = gated(4)
        4.times { ctx.backlog.taskPending() }

        when:
        def first = sourceThread(ctx.bp, ctx.source)
        def second = sourceThread(ctx.bp, ctx.source)
        await { ctx.gate().getWaiters() == 2 }
        then:
        ctx.gate().getBlockCount() == 2

        when:
        // an interrupted waiter leaves while the other keeps waiting
        first.interrupt()
        first.join(5_000)
        then:
        !first.isAlive()
        ctx.gate().getWaiters() == 1
        ctx.gate().isWaiting()
        second.isAlive()

        when:
        ctx.bp.release()
        second.join(5_000)
        then:
        !second.isAlive()
        ctx.gate().getWaiters() == 0

        cleanup:
        ctx.bp.release()
    }

    def 'should deregister a closed process' () {
        given:
        def ctx = gated(4)

        expect:
        ctx.bp.processes() == [ctx.backlog]

        when:
        ctx.backlog.close()
        then:
        ctx.bp.processes().isEmpty()
    }

    def 'should not gate value channels nor topics' () {
        given:
        def bp = new Backpressure(null, 4)
        def value = new DataflowVariable()
        def topic = CH.topic('nums')

        when:
        bp.markSource('channel.value', value)
        bp.markSource('channel.topic', topic)
        then:
        bp.gate(value) == null
        bp.gate(topic) == null

        cleanup:
        CH.init()
    }

    def 'should keep the gate of a source marked again' () {
        given:
        def ctx = gated(4)
        def gate = ctx.gate()

        when:
        ctx.bp.markSource('channel.of', ctx.source)
        then:
        ctx.gate().is(gate)
        ctx.gate().getFeeds() == [ctx.backlog]
    }

    def 'should count the queued input sets for the source gate' () {
        given:
        def ctx = gated(4)

        when:
        for( int i = 0; i < 4; i++ )
            ctx.source.bind(i)
        then:
        !ctx.bp.isOpen(ctx.source)

        when:
        ctx.backlog.close()
        then:
        ctx.bp.isOpen(ctx.source)
    }

    def 'should release the blocked sources' () {
        given:
        def ctx = gated(4)
        4.times { ctx.backlog.taskPending() }

        when:
        def thread = sourceThread(ctx.bp, ctx.source)
        await { ctx.gate().isWaiting() }
        then:
        thread.isAlive()

        when:
        ctx.bp.release()
        thread.join(5_000)
        then:
        !thread.isAlive()
        ctx.bp.isOpen(ctx.source)
    }

    def 'should stop emitting the remaining items once released' () {
        given:
        def ctx = gated(8)
        def reader = ctx.source.createReadChannel()

        when:
        ctx.bp.emit(ctx.source, range(1, 100), true)
        await { ctx.gate().isWaiting() }
        ctx.bp.release()
        ctx.bp.shutdown()
        then:
        ctx.bp.pump.awaitTermination(5, TimeUnit.SECONDS)
        reader.length() == 8

        cleanup:
        ctx.bp.shutdown()
    }

    def 'should emit on the calling thread while open and hand the rest to a source thread' () {
        given:
        def ctx = gated(8)
        def reader = ctx.source.createReadChannel()

        when:
        ctx.bp.emit(ctx.source, range(1, 100), true)
        then:
        // the caller returns once the process backlog reaches the high-water mark
        ctx.backlog.ready(1000) == 8

        when:
        // a short wait to catch an item emitted past the mark: it can only miss a bug, never cause a false failure
        Thread.sleep(50)
        then:
        ctx.backlog.ready(1000) == 8

        when:
        ctx.backlog.close()
        def result = []
        for( int i = 0; i < 100; i++ )
            result.add(reader.getVal())
        then:
        result == range(1, 100)
        reader.getVal().is(Channel.STOP)

        cleanup:
        ctx.bp.shutdown()
    }

    def 'should emit nothing for a null iterable' () {
        given:
        def bp = new Backpressure(null, 2)
        def ch = new DataflowBroadcast()
        def reader = ch.createReadChannel()

        when:
        bp.emit(ch, null, true)
        then:
        reader.getVal().is(Channel.STOP)
    }

    def 'should emit everything on the calling thread for an ungated channel' () {
        given:
        def bp = new Backpressure(null, 2)
        def ch = new DataflowBroadcast()
        def reader = ch.createReadChannel()

        when:
        bp.emit(ch, range(1, 10), true)
        then:
        reader.length() == 11
    }
}
