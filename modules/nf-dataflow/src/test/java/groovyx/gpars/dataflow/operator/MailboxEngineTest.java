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
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import groovyx.gpars.Closures;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.api.condition.DisabledIfEnvironmentVariable;

import static groovyx.gpars.TestThreads.poll;
import static groovyx.gpars.TestThreads.await;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Operator engine behaviors beyond the GPars conformance suite.
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class MailboxEngineTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    static class StopCounter implements DataflowEventListener {
        final AtomicInteger stops = new AtomicInteger();

        @Override
        public void afterStop(DataflowProcessor processor) {
            stops.incrementAndGet();
        }
    }

    @Test
    void shouldBeActiveFromStartUntilItStops() throws Exception {
        var input = new DataflowQueue<Object>();
        var op = new DataflowProcessor(List.of(input), List.of(), List.of(), 1, Closures.action(x -> {}));
        assertThat(op.isActive()).isFalse();

        op.start();
        assertThat(op.isActive()).isTrue();

        input.bind(STOP);
        op.join();
        assertThat(op.isActive()).isFalse();
    }

    @Test
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
    void shouldStopTakingInputsWhileAllForksAreBusy() throws Exception {
        var input = new DataflowQueue<Object>();
        var release = new CountDownLatch(1);
        var started = new AtomicInteger();
        var stopped = new StopCounter();

        new DataflowProcessor(List.of(input), List.of(), List.of(stopped), 2, Closures.action(x -> {
            started.incrementAndGet();
            release.await();
        })).start();
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
    void shouldHoldAForkUntilAnAsynchronousBodyCompletes() throws Exception {
        for( int forks : new int[] {1, 2} ) {
            var input = new DataflowQueue<Object>();
            var bodies = new CopyOnWriteArrayList<CompletableFuture<Object>>();
            var afterRun = new AtomicInteger();
            var stopped = new StopCounter() {
                @Override
                public void afterRun(DataflowProcessor processor, List<Object> messages) {
                    afterRun.incrementAndGet();
                }
            };

            new DataflowProcessor(List.of(input), List.of(), List.of(stopped), forks, Closures.function(x -> {
                var result = new CompletableFuture<Object>();
                bodies.add(result);
                return result;
            })).setAsync(true).start();
            for( int i = 0; i < 5; i++ )
                input.bind(i);

            // the bodies return at once, but each holds its fork until its future completes
            await(() -> bodies.size() == forks);
            Thread.sleep(50);
            assertThat(bodies).hasSize(forks);
            assertThat(afterRun.get()).isZero();

            // the poison pill waits for the bodies in flight
            input.bind(STOP);
            for( int i = 0; i < 5; i++ ) {
                final int n = i + 1;
                await(() -> bodies.size() >= n);
                bodies.get(i).complete(null);
            }
            await(() -> stopped.stops.get() == 1);
            assertThat(afterRun.get()).isEqualTo(5);
        }
    }

    @Test
    void shouldReportTheErrorOfAnAsynchronousBody() throws Exception {
        var input = new DataflowQueue<Object>();
        var errors = new CopyOnWriteArrayList<Throwable>();
        var listener = new DataflowEventListener() {
            @Override
            public boolean onException(DataflowProcessor processor, Throwable error) {
                errors.add(error);
                return false;
            }
        };
        var failure = new IllegalStateException("boom");

        var op = new DataflowProcessor(List.of(input), List.of(), List.of(listener), 2, Closures.function(x ->
            CompletableFuture.runAsync(() -> { throw failure; })
        )).setAsync(true).start();
        input.bind(1);

        await(() -> errors.size() == 1);
        assertThat(errors.get(0)).isSameAs(failure);
        assertThat(op.isActive()).isTrue();
        input.bind(STOP);
        op.join();
    }

    @Test
    void shouldForwardAnIncompleteFutureUntouchedUnlessAsync() throws Exception {
        for( int forks : new int[] {1, 2} ) {
            var input = new DataflowQueue<Object>();
            var output = new DataflowQueue<Object>();
            var stopped = new StopCounter();
            var item = new CompletableFuture<Object>();

            // like the forwarding operators, the body binds its input and returns it
            new DataflowProcessor(List.of(input), List.of(output), List.of(stopped), forks, Closures.closure(1, (self, args) -> {
                ((DataflowProcessor) self).bindOutput(args[0]);
                return args[0];
            })).start();
            input.bind(item);
            input.bind(1);
            input.bind(STOP);

            await(() -> stopped.stops.get() == 1);
            var received = new ArrayList<Object>();
            for( int i = 0; i < 3; i++ )
                received.add(output.getVal());
            assertThat(received.subList(0, 2)).containsExactlyInAnyOrder(item, 1);
            assertThat(received.get(2)).isSameAs(STOP);
            assertThat(item).isNotDone();
        }
    }

    @Test
    void shouldFireAfterStopOnceAndReleaseJoinWhenTerminatedConcurrently() throws Exception {
        for( int round = 0; round < 50; round++ ) {
            var input = new DataflowQueue<Object>();
            var stopped = new StopCounter();
            var op = new DataflowProcessor(List.of(input), List.of(), List.of(stopped), 1, Closures.action(x -> {})).start();
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
        var pool = Executors.newFixedThreadPool(1);
        try {
            var started = new CountDownLatch(1);
            var sleeper = new DataflowProcessor(pool, List.of(constant(1)), List.of(), List.of(), 1, Closures.action(x -> {
                started.countDown();
                Thread.sleep(30_000);
            })).start();
            assertThat(started.await(5, TimeUnit.SECONDS)).isTrue();
            sleeper.terminate();
            sleeper.join();

            var interrupted = new DataflowVariable<Boolean>();
            new DataflowProcessor(pool, List.of(constant(2)), List.of(), List.of(), 1, Closures.action(x -> interrupted.bind(Thread.currentThread().isInterrupted()))).start();
            assertThat(poll(interrupted, 5, TimeUnit.SECONDS)).isFalse();
        }
        finally {
            pool.shutdown();
        }
    }

    @Test
    void shouldNotScheduleAnUnstartedOperator() throws Exception {
        var executions = new AtomicInteger();
        var pool = Executors.newFixedThreadPool(2);
        Executor counting = task -> {
            executions.incrementAndGet();
            pool.execute(task);
        };
        try {
            var input = new DataflowQueue<Object>();
            var fired = new CopyOnWriteArrayList<Object>();
            var op = new DataflowProcessor(counting, List.of(input), List.of(), List.of(), 1, Closures.action(x -> fired.add(x)));
            input.bind(1);
            Thread.sleep(100);
            assertThat(executions.get()).isZero();

            op.start();
            await(() -> fired.size() == 1);
            op.terminate();
            op.join();
        }
        finally {
            pool.shutdown();
        }
    }

    @Test
    void shouldClearTheInterruptOfAnExternalTerminateBeforeAfterStop() throws Exception {
        var started = new CountDownLatch(1);
        var interrupted = new DataflowVariable<Boolean>();
        var listener = new DataflowEventListener() {
            @Override
            public void afterStop(DataflowProcessor processor) {
                interrupted.bind(Thread.currentThread().isInterrupted());
            }
        };
        var op = new DataflowProcessor(List.of(constant(1)), List.of(), List.of(listener), 1, Closures.action(x -> {
            started.countDown();
            // wait for the interrupt of the terminate, leaving the flag set
            final long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(5);
            while( !Thread.currentThread().isInterrupted() && System.nanoTime() < deadline )
                Thread.onSpinWait();
        })).start();
        assertThat(started.await(5, TimeUnit.SECONDS)).isTrue();
        op.terminate();
        assertThat(poll(interrupted, 5, TimeUnit.SECONDS)).isFalse();
    }

    @Test
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
    void shouldForwardThePoisonPillWhenTerminatedWhileWaitingForTheBodiesInFlight() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var started = new CountDownLatch(1);
        var poisoned = new CountDownLatch(1);
        var release = new CountDownLatch(1);
        var listener = new DataflowEventListener() {
            @Override
            public Object controlMessageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
                poisoned.countDown();
                return message;
            }
        };
        var op = new DataflowProcessor(List.of(input), List.of(output), List.of(listener), 2, Closures.action(x -> {
            started.countDown();
            release.await();
        })).start();
        try {
            input.bind(1);
            assertThat(started.await(5, TimeUnit.SECONDS)).isTrue();
            input.bind(STOP);
            assertThat(poisoned.await(5, TimeUnit.SECONDS)).isTrue();

            op.terminate();
            op.join();
            assertThat(poll(output, 5, TimeUnit.SECONDS)).isInstanceOf(PoisonPill.class);
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

        var op = new DataflowProcessor(List.of(input), List.of(), List.of(stopped), 2, Closures.action(x -> {
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
        var op = new DataflowProcessor(List.of(a, b), List.of(), List.of(), 1, Closures.action((x, y) -> fired.incrementAndGet())).start();
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
        var op = new DataflowProcessor(List.of(input), List.of(), List.of(), 1, Closures.action(x -> fired.add(x))).start();
        input.bind(1);
        await(() -> fired.size() == 1);

        op.terminate();
        op.join();
        input.bind(2);
        assertThat(poll(input, 5, TimeUnit.SECONDS)).isEqualTo(2);
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
        var op = new DataflowProcessor(List.of(input), List.of(), List.of(), 1, Closures.action(x -> fired.incrementAndGet())).start();
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
        var pool = Executors.newFixedThreadPool(1);
        try {
            var busy = new DataflowBroadcast<Object>();
            var reader = busy.createReadChannel();
            for( int i = 0; i < n; i++ )
                busy.bind(i);
            busy.bind(STOP);
            var count = new AtomicInteger();
            var seen = new DataflowVariable<Integer>();
            var stopped = new StopCounter();

            new DataflowProcessor(pool, List.of(reader), List.of(), List.of(stopped), 1, Closures.action(x -> count.incrementAndGet())).start();
            new DataflowProcessor(pool, List.of(constant(1)), List.of(), List.of(), 1, Closures.action(x -> seen.bind(count.get()))).start();

            assertThat(poll(seen, 5, TimeUnit.SECONDS)).isLessThan(n);
            await(() -> stopped.stops.get() == 1);
            assertThat(count.get()).isEqualTo(n);
        }
        finally {
            pool.shutdown();
        }
    }

    private static DataflowQueue<Object> constant(Object value) {
        var result = new DataflowQueue<Object>();
        result.bind(value);
        return result;
    }
}
