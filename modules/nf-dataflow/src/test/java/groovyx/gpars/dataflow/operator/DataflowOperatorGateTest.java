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
package groovyx.gpars.dataflow.operator;

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BooleanSupplier;

import groovyx.gpars.Closures;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowVariable;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static groovyx.gpars.TestThreads.poll;
import static groovyx.gpars.TestThreads.await;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * The backpressure gate of a dataflow operator ({@link DataflowProcessor#setGate}, {@link DataflowProcessor#resume}).
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class DataflowOperatorGateTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    /**
     * An operator with a queue input and a value input, which records the queue items it fires on.
     */
    private static DataflowProcessor operator(DataflowQueue<Object> input, List<Object> fired, int maxForks, DataflowVariable<Object> stopped, AtomicInteger started) {
        final var listener = new DataflowEventListener() {
            @Override
            public List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
                if( started != null )
                    started.incrementAndGet();
                return messages;
            }

            @Override
            public void afterStop(DataflowProcessor processor) {
                if( stopped != null )
                    stopped.bind(true);
            }
        };
        final var value = new DataflowVariable<Object>();
        value.bind(0);
        return new DataflowProcessor(List.of(input, value), List.of(), List.of(listener), maxForks, Closures.action((x, y) -> fired.add(x)));
    }

    /**
     * Busy-wait until a condition holds, so as to race with the operator thread.
     */
    private static void spinUntil(BooleanSupplier condition) {
        final long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(5);
        while( !condition.getAsBoolean() ) {
            if( System.nanoTime() > deadline )
                throw new AssertionError("Condition not met in time");
            Thread.onSpinWait();
        }
    }

    private static void send(DataflowQueue<Object> input, int n) {
        for( int i = 0; i < n; i++ )
            input.bind(i);
    }

    @ParameterizedTest
    @ValueSource(ints = { 1, 4 })
    void shouldPauseAfterAFiringWhileTheGateIsClosedAndResumeOnceItOpens(int maxForks) throws Exception {
        var fired = new CopyOnWriteArrayList<Object>();
        var closed = new AtomicBoolean(true);
        var started = new AtomicInteger();
        var input = new DataflowQueue<Object>();
        var op = operator(input, fired, maxForks, null, started);
        // like a process backlog, the gate counts the input sets taken in beforeRun
        op.setGate(() -> started.get() >= 3 && closed.get());
        op.start();
        try {
            send(input, 10);
            await(() -> fired.size() == 3 && op.isPaused());
            assertThat(op.getPauseCount()).isEqualTo(1);
            // a short wait to catch an extra firing: it can only miss a bug, never cause a false failure
            Thread.sleep(50);
            assertThat(fired).hasSize(3);

            // a resume while the gate is still closed has no effect
            op.resume();
            Thread.sleep(50);
            assertThat(fired).hasSize(3);
            assertThat(op.isPaused()).isTrue();

            closed.set(false);
            op.resume();
            await(() -> fired.size() == 10);
            assertThat(op.isPaused()).isFalse();
        }
        finally {
            op.terminate();
        }
    }

    @Test
    void shouldTakeThePoisonPillOnceResumed() throws Exception {
        var fired = new CopyOnWriteArrayList<Object>();
        var closed = new AtomicBoolean(true);
        var stopped = new DataflowVariable<Object>();
        var input = new DataflowQueue<Object>();
        var op = operator(input, fired, 1, stopped, null);
        op.setGate(closed::get);
        op.start();

        input.bind(1);
        input.bind(2);
        input.bind(STOP);
        await(() -> fired.size() == 1 && op.isPaused());
        assertThat(stopped.isBound()).isFalse();

        closed.set(false);
        op.resume();
        assertThat(poll(stopped, 5, TimeUnit.SECONDS)).isEqualTo(true);
        assertThat(fired).containsExactly(1, 2);
    }

    @Test
    void shouldNotPauseWithoutAGate() throws Exception {
        var fired = new CopyOnWriteArrayList<Object>();
        var input = new DataflowQueue<Object>();
        var op = operator(input, fired, 1, null, null);
        op.start();
        try {
            send(input, 5);
            op.resume();
            await(() -> fired.size() == 5);
            assertThat(op.isPaused()).isFalse();
            assertThat(op.getPauseCount()).isZero();
        }
        finally {
            op.terminate();
        }
    }

    @ParameterizedTest
    @ValueSource(ints = { 1, 4 })
    void shouldNotLoseAResumeWhileTheGateIsToggledConcurrently(int maxForks) throws Exception {
        final int n = 5_000;
        var fired = new CopyOnWriteArrayList<Object>();
        var started = new AtomicInteger();
        var allowed = new AtomicInteger(1);
        var stopped = new DataflowVariable<Object>();
        var input = new DataflowQueue<Object>();
        var op = operator(input, fired, maxForks, stopped, started);
        // like a process backlog, the gate closes once the allowed input sets are taken
        op.setGate(() -> started.get() >= allowed.get());
        op.start();
        send(input, n);
        input.bind(STOP);

        // open the gate by one input set at a time, racing with the operator pausing after a
        // firing: a single resume per opening, so a lost one stalls the operator. The last
        // opening lets the operator take the poison pill
        for( int i = 1; i <= n; i++ ) {
            final int expected = i;
            spinUntil(() -> started.get() >= expected);
            allowed.incrementAndGet();
            op.resume();
        }

        assertThat(poll(stopped, 5, TimeUnit.SECONDS)).isEqualTo(true);
        assertThat(fired).hasSize(n);
        assertThat(op.getPauseCount()).isPositive();
    }

    @Test
    void shouldCountTheItemsBufferedForABroadcastReader() throws Exception {
        var broadcast = new DataflowBroadcast<Object>();
        var reader = broadcast.createReadChannel();

        assertThat(reader.length()).isZero();

        for( int i = 0; i < 5; i++ )
            broadcast.bind(i);
        assertThat(reader.length()).isEqualTo(5);

        reader.getVal();
        reader.getVal();
        assertThat(reader.length()).isEqualTo(3);
    }
}
