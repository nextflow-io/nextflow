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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import groovyx.gpars.Closures;
import groovyx.gpars.dataflow.Dataflow;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.group.DefaultPGroup;
import groovyx.gpars.scheduler.DefaultPool;
import groovyx.gpars.scheduler.Pool;
import groovyx.gpars.scheduler.Scheduler;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.TestThreads.await;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Operator engine behaviors beyond the GPars conformance suite. They hold on both engines.
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class OperatorEngineTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    static class StopCounter extends DataflowEventAdapter {
        final AtomicInteger stops = new AtomicInteger();

        @Override
        public void afterStop(DataflowProcessor processor) {
            stops.incrementAndGet();
        }
    }

    @Test
    void shouldRunOnTheEngineSelectedByTheEnvironment() {
        var op = new DataflowOperator(Dataflow.retrieveCurrentDFPGroup(), Map.of("inputs", List.of(new DataflowQueue<>())), Closures.action(x -> {}));
        assertThat(op.engine).isInstanceOf(Scheduler.isLegacyEngine() ? ActorEngine.class : MailboxEngine.class);
    }

    @Test
    void shouldBeActiveFromStartUntilItStops() throws Exception {
        var input = new DataflowQueue<Object>();
        var op = new DataflowOperator(Dataflow.retrieveCurrentDFPGroup(), Map.of("inputs", List.of(input)), Closures.action(x -> {}));
        assertThat(op.isActive()).isFalse();

        op.start();
        assertThat(op.isActive()).isTrue();

        input.bind(STOP);
        op.join();
        assertThat(op.isActive()).isFalse();
    }

    @Test
    void shouldStopTakingInputsWhileAllForksAreBusy() throws Exception {
        var input = new DataflowQueue<Object>();
        var release = new CountDownLatch(1);
        var started = new AtomicInteger();
        var stopped = new StopCounter();

        Dataflow.operator(Map.of("inputs", List.of(input), "maxForks", 2, "listeners", List.of(stopped)), Closures.action(x -> {
            started.incrementAndGet();
            release.await();
        }));
        for( int i = 0; i < 5; i++ )
            input.bind(i);

        // two bodies run, a third input set waits for one of them, the rest stays in the channel
        await(() -> started.get() == 2 && input.length() == 2);
        Thread.sleep(50);
        assertThat(started.get()).isEqualTo(2);
        assertThat(input.length()).isEqualTo(2);

        release.countDown();
        input.bind(STOP);
        await(() -> stopped.stops.get() == 1);
        assertThat(started.get()).isEqualTo(5);
    }

    @Test
    void shouldFireAfterStopOnceAndReleaseJoinWhenTerminatedConcurrently() throws Exception {
        for( int round = 0; round < 50; round++ ) {
            var input = new DataflowQueue<Object>();
            var stopped = new StopCounter();
            var op = Dataflow.operator(Map.of("inputs", List.of(input), "listeners", List.of(stopped)), Closures.action(x -> {}));
            var go = new CountDownLatch(1);
            var threads = new ArrayList<Thread>();
            for( int i = 0; i < 4; i++ )
                threads.add(start(() -> { go.await(); op.terminate(); }));
            input.bind(1);
            go.countDown();
            for( var thread : threads )
                thread.join();

            op.join();
            assertThat(op.isActive()).isFalse();
            assertThat(stopped.stops.get()).isEqualTo(1);
        }
    }

    @Test
    void shouldNotLeakTheInterruptOfAnExternalTerminateToTheNextTaskOnThePoolThread() throws Exception {
        // a single thread, so the next operator runs on the thread interrupted for the first one
        var group = new DefaultPGroup(new DefaultPool(true, 1));
        try {
            var started = new CountDownLatch(1);
            var sleeper = group.operator(Map.of("inputs", List.of(constant(1))), Closures.action(x -> {
                started.countDown();
                Thread.sleep(30_000);
            }));
            assertThat(started.await(5, TimeUnit.SECONDS)).isTrue();
            sleeper.terminate();
            sleeper.join();

            var interrupted = new DataflowVariable<Boolean>();
            group.operator(Map.of("inputs", List.of(constant(2))), Closures.action(x -> interrupted.bind(Thread.currentThread().isInterrupted())));
            assertThat(interrupted.getVal(5, TimeUnit.SECONDS)).isFalse();
        }
        finally {
            group.shutdown();
        }
    }

    @Test
    void shouldNotScheduleAnUnstartedOperator() throws Exception {
        var executions = new AtomicInteger();
        var pool = new DefaultPool(true, 2);
        var counting = new Pool() {
            public void resize(int poolSize) { pool.resize(poolSize); }
            public void resetDefaultSize() { pool.resetDefaultSize(); }
            public int getPoolSize() { return pool.getPoolSize(); }
            public void shutdown() { pool.shutdown(); }
            public void execute(Runnable task) {
                executions.incrementAndGet();
                pool.execute(task);
            }
        };
        var group = new DefaultPGroup(counting);
        try {
            var input = new DataflowQueue<Object>();
            var fired = new CopyOnWriteArrayList<Object>();
            var op = new DataflowOperator(group, Map.of("inputs", List.of(input)), Closures.action(x -> fired.add(x)));
            op.terminateAfterNextRun();
            Thread.sleep(100);
            assertThat(executions.get()).isZero();

            op.start();
            input.bind(1);
            await(() -> fired.size() == 1);
            op.terminate();
            op.join();
        }
        finally {
            group.shutdown();
        }
    }

    @Test
    void shouldClearTheInterruptOfAnExternalTerminateBeforeAfterStop() throws Exception {
        var started = new CountDownLatch(1);
        var interrupted = new DataflowVariable<Boolean>();
        var listener = new DataflowEventAdapter() {
            @Override
            public void afterStop(DataflowProcessor processor) {
                interrupted.bind(Thread.currentThread().isInterrupted());
            }
        };
        var op = Dataflow.operator(Map.of("inputs", List.of(constant(1)), "listeners", List.of(listener)), Closures.action(x -> {
            started.countDown();
            // wait for the interrupt of the terminate, leaving the flag set
            final long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(5);
            while( !Thread.currentThread().isInterrupted() && System.nanoTime() < deadline )
                Thread.onSpinWait();
        }));
        assertThat(started.await(5, TimeUnit.SECONDS)).isTrue();
        op.terminate();
        assertThat(interrupted.getVal(5, TimeUnit.SECONDS)).isFalse();
    }

    @Test
    void shouldForwardThePoisonPillWhenTerminatedWhileWaitingForTheBodiesInFlight() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var started = new CountDownLatch(1);
        var poisoned = new CountDownLatch(1);
        var release = new CountDownLatch(1);
        var listener = new DataflowEventAdapter() {
            @Override
            public Object controlMessageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
                poisoned.countDown();
                return message;
            }
        };
        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "maxForks", 2, "listeners", List.of(listener)), Closures.action(x -> {
            started.countDown();
            release.await();
        }));
        try {
            input.bind(1);
            assertThat(started.await(5, TimeUnit.SECONDS)).isTrue();
            input.bind(STOP);
            assertThat(poisoned.await(5, TimeUnit.SECONDS)).isTrue();

            op.terminate();
            op.join();
            assertThat(output.getVal(5, TimeUnit.SECONDS)).isInstanceOf(PoisonPill.class);
        }
        finally {
            release.countDown();
        }
    }

    @Test
    void shouldPauseWhenAParkedInputSetForksWhileTheGateIsClosed() throws Exception {
        var input = new DataflowQueue<Object>();
        var release = new CountDownLatch(1);
        var started = new AtomicInteger();
        var closed = new AtomicBoolean();
        var stopped = new StopCounter();

        var op = new DataflowOperator(Dataflow.retrieveCurrentDFPGroup(), Map.of("inputs", List.of(input), "maxForks", 2, "listeners", List.of(stopped)), Closures.action(x -> {
            started.incrementAndGet();
            release.await();
        }));
        op.setGate(closed::get);
        op.start();
        for( int i = 0; i < 5; i++ )
            input.bind(i);

        // two bodies run and a third input set is parked: the gate is checked again once it forks
        await(() -> started.get() == 2 && input.length() == 2);
        closed.set(true);
        release.countDown();
        await(() -> started.get() == 3 && op.isPaused());
        Thread.sleep(50);
        assertThat(started.get()).isEqualTo(3);
        assertThat(input.length()).isEqualTo(2);
        assertThat(op.getPauseCount()).isEqualTo(1);

        closed.set(false);
        op.resume();
        input.bind(STOP);
        await(() -> stopped.stops.get() == 1);
        assertThat(started.get()).isEqualTo(5);
    }

    private static DataflowQueue<Object> constant(Object value) {
        var result = new DataflowQueue<Object>();
        result.bind(value);
        return result;
    }
}
