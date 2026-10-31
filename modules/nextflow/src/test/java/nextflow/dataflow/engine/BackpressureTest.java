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
package nextflow.dataflow.engine;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import nextflow.Channel;
import nextflow.SysEnv;
import nextflow.extension.CH;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static nextflow.dataflow.engine.BackpressureInternals.await;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Unit tests of the source gates.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(value = 20, unit = TimeUnit.SECONDS)
class BackpressureTest {

    /**
     * A backpressure with a source channel feeding one process.
     */
    record Gated(Backpressure bp, DataflowBroadcast<Object> source, ProcessBacklog backlog) {

        Backpressure.Gate gate() {
            return bp.gate(source);
        }
    }

    @SuppressWarnings("rawtypes")
    private static Gated gated(int highWater) {
        final var bp = new Backpressure(null, highWater);
        final var source = new DataflowBroadcast<Object>();
        final var backlog = new ProcessBacklog(bp, "foo", List.<DataflowReadChannel>of(source.createReadChannel()), null);
        bp.markSource("channel.of", source);
        bp.gate(source).setFeeds(List.of(backlog));
        bp.addProcess(backlog);
        return new Gated(bp, source, backlog);
    }

    private static Thread sourceThread(Backpressure bp, Object source) {
        final var thread = new Thread(() -> Backpressure.runAsSource(() -> bp.awaitOpen(source)));
        thread.start();
        return thread;
    }

    private static void times(int n, Runnable action) {
        for( int i = 0; i < n; i++ )
            action.run();
    }

    private static List<Object> range(int from, int to) {
        return new ArrayList<>(IntStream.rangeClosed(from, to).boxed().toList());
    }

    @Test
    void shouldGetTheHighWaterMark() {
        SysEnv.push(Map.of());
        try {
            assertThat(Backpressure.highWater(4)).isEqualTo(64);
            assertThat(Backpressure.highWater(16)).isEqualTo(64);
            assertThat(Backpressure.highWater(50)).isEqualTo(100);
        }
        finally {
            SysEnv.pop();
        }

        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", "16"));
        try {
            assertThat(Backpressure.highWater(50)).isEqualTo(16);
            assertThat(new Backpressure(null, 16).getLowWater()).isEqualTo(8);
        }
        finally {
            SysEnv.pop();
        }

        // zero disables backpressure
        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", "0"));
        try {
            assertThat(Backpressure.highWater(4)).isZero();
        }
        finally {
            SysEnv.pop();
        }
    }

    @ParameterizedTest
    @ValueSource(strings = {"-1", "abc", "1.5", "true"})
    void shouldRejectAnInvalidHighWaterMark(String value) {
        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", value));
        try {
            assertThatThrownBy(() -> Backpressure.highWater(4))
                .isInstanceOf(IllegalArgumentException.class)
                .hasMessageContaining("Invalid NXF_DATAFLOW_HIGH_WATER value: '" + value + "'");
        }
        finally {
            SysEnv.pop();
        }
    }

    @Test
    void shouldNeverBlockAThreadNotFlaggedAsASourceThread() {
        var ctx = gated(4);
        times(4, ctx.backlog()::taskPending);
        assertThat(ctx.bp().isOpen(ctx.source())).isFalse();

        // returns immediately, otherwise the test times out
        ctx.bp().awaitOpen(ctx.source());
        Backpressure.awaitSource(ctx.source());
    }

    @Test
    void shouldBlockASourceThreadFromTheHighWaterMarkDownToTheLowWaterMark() throws Exception {
        var ctx = gated(8);
        times(8, ctx.backlog()::taskPending);

        var thread = sourceThread(ctx.bp(), ctx.source());
        await(() -> ctx.gate().isWaiting());
        assertThat(ctx.gate().getBlockCount()).isEqualTo(1);

        // the backlog drops below the high-water mark but stays above the low-water mark
        times(3, ctx.backlog()::taskSubmitted);
        // a short wait to catch an early wake-up: it can only miss a bug, never cause a false failure
        Thread.sleep(50);
        assertThat(thread.isAlive()).isTrue();

        ctx.backlog().taskSubmitted();
        thread.join(5_000);
        assertThat(thread.isAlive()).isFalse();
        assertThat(ctx.gate().isWaiting()).isFalse();
    }

    @Test
    void shouldCountTheThreadsWaitingOnAGate() throws Exception {
        var ctx = gated(4);
        times(4, ctx.backlog()::taskPending);

        var first = sourceThread(ctx.bp(), ctx.source());
        var second = sourceThread(ctx.bp(), ctx.source());
        try {
            await(() -> ctx.gate().getWaiters() == 2);
            assertThat(ctx.gate().getBlockCount()).isEqualTo(2);

            // an interrupted waiter leaves while the other keeps waiting
            first.interrupt();
            first.join(5_000);
            assertThat(first.isAlive()).isFalse();
            assertThat(ctx.gate().getWaiters()).isEqualTo(1);
            assertThat(ctx.gate().isWaiting()).isTrue();
            assertThat(second.isAlive()).isTrue();
        }
        finally {
            ctx.bp().release();
        }
        second.join(5_000);
        assertThat(second.isAlive()).isFalse();
        assertThat(ctx.gate().getWaiters()).isZero();
    }

    @Test
    void shouldDeregisterAClosedProcess() {
        var ctx = gated(4);
        assertThat(ctx.bp().processes()).containsExactly(ctx.backlog());

        ctx.backlog().close();
        assertThat(ctx.bp().processes()).isEmpty();
    }

    @Test
    void shouldNotGateValueChannelsNorTopics() {
        var bp = new Backpressure(null, 4);
        var value = new DataflowVariable<Object>();
        var topic = CH.topic("nums");
        try {
            bp.markSource("channel.value", value);
            bp.markSource("channel.topic", topic);
            assertThat(bp.gate(value)).isNull();
            assertThat(bp.gate(topic)).isNull();
        }
        finally {
            CH.init();
        }
    }

    @Test
    void shouldKeepTheGateOfASourceMarkedAgain() {
        var ctx = gated(4);
        var gate = ctx.gate();

        ctx.bp().markSource("channel.of", ctx.source());
        assertThat(ctx.gate()).isSameAs(gate);
        assertThat(ctx.gate().getFeeds()).containsExactly(ctx.backlog());
    }

    @Test
    void shouldCountTheQueuedInputSetsForTheSourceGate() {
        var ctx = gated(4);

        for( int i = 0; i < 4; i++ )
            ctx.source().bind(i);
        assertThat(ctx.bp().isOpen(ctx.source())).isFalse();

        ctx.backlog().close();
        assertThat(ctx.bp().isOpen(ctx.source())).isTrue();
    }

    @Test
    void shouldReleaseTheBlockedSources() throws Exception {
        var ctx = gated(4);
        times(4, ctx.backlog()::taskPending);

        var thread = sourceThread(ctx.bp(), ctx.source());
        await(() -> ctx.gate().isWaiting());
        assertThat(thread.isAlive()).isTrue();

        ctx.bp().release();
        thread.join(5_000);
        assertThat(thread.isAlive()).isFalse();
        assertThat(ctx.bp().isOpen(ctx.source())).isTrue();
    }

    @Test
    void shouldEmitOnTheCallingThreadWhileOpenAndHandTheRestToASourceThread() throws Exception {
        var ctx = gated(8);
        var reader = ctx.source().createReadChannel();
        try {
            ctx.bp().emit(ctx.source(), range(1, 100), true);
            // the caller returns once the process backlog reaches the high-water mark
            assertThat(ctx.backlog().ready(1000)).isEqualTo(8);
            // a short wait to catch an item emitted past the mark: it can only miss a bug, never cause a false failure
            Thread.sleep(50);
            assertThat(ctx.backlog().ready(1000)).isEqualTo(8);

            ctx.backlog().close();
            var result = new ArrayList<Object>();
            for( int i = 0; i < 100; i++ )
                result.add(reader.getVal());
            assertThat(result).isEqualTo(range(1, 100));
            assertThat(reader.getVal()).isSameAs(Channel.STOP);
        }
        finally {
            ctx.bp().shutdown();
        }
    }

    @Test
    void shouldEmitNothingForANullIterable() throws Exception {
        var bp = new Backpressure(null, 2);
        var ch = new DataflowBroadcast<Object>();
        var reader = ch.createReadChannel();

        bp.emit(ch, null, true);
        assertThat(reader.getVal()).isSameAs(Channel.STOP);
    }

    @Test
    void shouldEmitEverythingOnTheCallingThreadForAnUngatedChannel() {
        var bp = new Backpressure(null, 2);
        var ch = new DataflowBroadcast<Object>();
        DataflowReadChannel<Object> reader = ch.createReadChannel();

        bp.emit(ch, range(1, 10), true);
        assertThat(reader.length()).isEqualTo(11);
    }
}
