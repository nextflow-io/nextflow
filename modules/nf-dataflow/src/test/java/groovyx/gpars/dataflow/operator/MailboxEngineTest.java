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
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import groovyx.gpars.group.DefaultPGroup;
import groovyx.gpars.scheduler.DefaultPool;
import groovyx.gpars.scheduler.Pool;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.api.condition.DisabledIfEnvironmentVariable;

import static groovyx.gpars.TestThreads.await;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Operator engine behaviors beyond the GPars conformance suite.
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class MailboxEngineTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    static class StopCounter extends DataflowEventAdapter {
        final AtomicInteger stops = new AtomicInteger();

        @Override
        public void afterStop(DataflowProcessor processor) {
            stops.incrementAndGet();
        }
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
    void shouldForwardASinglePoisonPillWhenSeveralInputsHandOverOneInTheSameDrain() throws Exception {
        // a manual pool, so that both pills are handed over before the drain that takes them
        var tasks = new ArrayList<Runnable>();
        var manual = new Pool() {
            public void resize(int poolSize) { }
            public void resetDefaultSize() { }
            public int getPoolSize() { return 1; }
            public void shutdown() { }
            public void execute(Runnable task) { tasks.add(task); }
        };
        var a = new DataflowQueue<Object>();
        var b = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var counter = new StopCounter();
        var op = new DataflowOperator(new DefaultPGroup(manual), Map.of("inputs", List.of(a, b), "outputs", List.of(output), "listeners", List.of(counter)), Closures.action((x, y) -> {}));

        op.start();
        // subscribe both inputs
        runAll(tasks);
        a.bind(STOP);
        b.bind(STOP);
        runAll(tasks);

        assertThat(op.isActive()).isFalse();
        assertThat(counter.stops.get()).isEqualTo(1);
        assertThat(output.getVal()).isSameAs(STOP);
        assertThat(output.length()).isZero();
    }

    private static void runAll(List<Runnable> tasks) {
        while( !tasks.isEmpty() )
            tasks.remove(0).run();
    }

    @Test
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
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

        try {
            // two bodies run, a third input set waits for one of them, the rest stays in the channel
            await(() -> started.get() == 2 && input.length() == 2);
            Thread.sleep(50);
            assertThat(started.get()).isEqualTo(2);
            assertThat(input.length()).isEqualTo(2);
        }
        finally {
            release.countDown();
        }
        input.bind(STOP);
        await(() -> stopped.stops.get() == 1);
        assertThat(started.get()).isEqualTo(5);
    }

    @Test
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
    void shouldPairBeforeRunAndAfterRunWhenAParkedInputSetIsDroppedOnStop() throws Exception {
        var input = new DataflowQueue<Object>();
        var release = new CountDownLatch(1);
        var started = new AtomicInteger();
        var before = new AtomicInteger();
        var after = new AtomicInteger();
        var stopped = new StopCounter();
        var listener = new DataflowEventAdapter() {
            @Override
            public List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
                before.incrementAndGet();
                return messages;
            }

            @Override
            public void afterRun(DataflowProcessor processor, List<Object> messages) {
                after.incrementAndGet();
            }
        };

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "maxForks", 2, "listeners", List.of(listener, stopped)), Closures.action(x -> {
            started.incrementAndGet();
            release.await();
        }));
        for( int i = 0; i < 5; i++ )
            input.bind(i);

        try {
            // two bodies run and a third input set is parked, then the operator is stopped
            await(() -> started.get() == 2 && input.length() == 2);
            op.terminate();
            await(() -> stopped.stops.get() == 1);
        }
        finally {
            release.countDown();
        }
        // the parked input set is dropped without a beforeRun, which would expect an afterRun
        await(() -> after.get() == 2);
        Thread.sleep(50);
        assertThat(before.get()).isEqualTo(2);
        assertThat(after.get()).isEqualTo(2);
        assertThat(started.get()).isEqualTo(2);
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
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
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
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
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

        try {
            // two bodies run and a third input set is parked: the gate is checked again once it forks
            await(() -> started.get() == 2 && input.length() == 2);
            closed.set(true);
        }
        finally {
            release.countDown();
        }
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

    @Test
    void shouldFireOnceWhenAllInputsAreValues() throws Exception {
        // a firing needs a value that arrived: the value inputs are re-read, which does not fire again
        var a = new DataflowVariable<Object>();
        var b = new DataflowVariable<Object>();
        a.bind(1);
        b.bind(2);
        var fired = new AtomicInteger();
        var op = Dataflow.operator(Map.of("inputs", List.of(a, b)), Closures.action((x, y) -> fired.incrementAndGet()));
        try {
            await(() -> fired.get() == 1);
            Thread.sleep(100);
            assertThat(fired.get()).isEqualTo(1);
            assertThat(op.isActive()).isTrue();
        }
        finally {
            op.terminate();
        }
    }

    @Test
    void shouldNotTakeAValueOfAQueueInputOnceStopped() throws Exception {
        var input = new DataflowQueue<Object>();
        var fired = new CopyOnWriteArrayList<Object>();
        var op = Dataflow.operator(Map.of("inputs", List.of(input)), Closures.action(x -> fired.add(x)));
        input.bind(1);
        await(() -> fired.size() == 1);

        op.terminate();
        op.join();
        input.bind(2);
        assertThat(input.getVal(5, TimeUnit.SECONDS)).isEqualTo(2);
        assertThat(fired).containsExactly(1);
    }

    @Test
    void shouldDropAtMostTheValueHandedOverWhenTerminatedWhileAWriterBinds() throws Exception {
        for( int round = 0; round < 300; round++ ) {
            var queue = new DataflowQueue<Object>();
            terminateWhileBinding(queue, queue, round);
            var broadcast = new DataflowBroadcast<Object>();
            terminateWhileBinding(broadcast, broadcast.createReadChannel(), round);
        }
    }

    private static void terminateWhileBinding(DataflowWriteChannel<Object> channel, DataflowReadChannel<Object> input, int round) throws Exception {
        final int n = 200;
        var fired = new AtomicInteger();
        var op = Dataflow.operator(Map.of("inputs", List.of(input)), Closures.action(x -> fired.incrementAndGet()));
        var writer = start(() -> {
            for( int i = 0; i < n; i++ )
                channel.bind(i);
        });
        for( int spin = round * 7 % 2000; spin > 0; spin-- )
            Thread.onSpinWait();
        op.terminate();
        op.join();
        writer.join();

        // each value is fired or left in the channel, except one handed over to the stopped operator
        final int firedBeforeStop = fired.get();
        assertThat(firedBeforeStop + input.length()).isBetween(n - 1, n);
        Thread.sleep(1);
        assertThat(fired.get()).isEqualTo(firedBeforeStop);
    }

    @Test
    void shouldYieldThePoolThreadWhileFiringOnBufferedValues() throws Exception {
        // a single thread: the other operator runs only if the first one yields between batches
        final int n = 100_000;
        var group = new DefaultPGroup(new DefaultPool(true, 1));
        try {
            var busy = new DataflowBroadcast<Object>();
            var reader = busy.createReadChannel();
            for( int i = 0; i < n; i++ )
                busy.bind(i);
            busy.bind(STOP);
            var count = new AtomicInteger();
            var seen = new DataflowVariable<Integer>();
            var stopped = new StopCounter();

            group.operator(Map.of("inputs", List.of(reader), "listeners", List.of(stopped)), Closures.action(x -> count.incrementAndGet()));
            group.operator(Map.of("inputs", List.of(constant(1))), Closures.action(x -> seen.bind(count.get())));

            assertThat(seen.getVal(5, TimeUnit.SECONDS)).isLessThan(n);
            await(() -> stopped.stops.get() == 1);
            assertThat(count.get()).isEqualTo(n);
        }
        finally {
            group.shutdown();
        }
    }

    private static DataflowQueue<Object> constant(Object value) {
        var result = new DataflowQueue<Object>();
        result.bind(value);
        return result;
    }
}
