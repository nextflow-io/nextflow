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
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import groovyx.gpars.Closures;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.api.condition.DisabledIfEnvironmentVariable;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

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
    void shouldForwardASinglePoisonPillWhenSeveralInputsHandOverOneInTheSameDrain() throws Exception {
        // a manual pool, so that both pills are handed over before the drain that takes them
        var tasks = new ArrayList<Runnable>();
        Executor manual = tasks::add;
        var a = new DataflowQueue<Object>();
        var b = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var counter = new StopCounter();
        var op = new DataflowProcessor(manual, List.of(a, b), List.of(output), List.of(counter), 1, Closures.action((x, y) -> {}));

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

    /**
     * Records the errors, the runs and the stop of an operator
     */
    static class Recorder extends StopCounter {
        final List<Throwable> errors = new CopyOnWriteArrayList<>();
        final AtomicInteger before = new AtomicInteger();
        final AtomicInteger after = new AtomicInteger();

        @Override
        public boolean onException(DataflowProcessor processor, Throwable e) {
            errors.add(e);
            return true;
        }

        @Override
        public List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
            before.incrementAndGet();
            return messages;
        }

        @Override
        public void afterRun(DataflowProcessor processor, List<Object> messages) {
            after.incrementAndGet();
        }
    }

    /**
     * The failures of a pool that cannot run a task: a rejection, or an error such as no thread could be created
     */
    static Stream<Throwable> poolFailures() {
        // not an OutOfMemoryError, which JUnit rethrows, aborting the run, if it fails a test
        return Stream.of(new RejectedExecutionException("rejected"), new Error("unable to create native thread"));
    }

    static RuntimeException raise(Throwable e) {
        if( e instanceof Error )
            throw (Error) e;
        throw (RuntimeException) e;
    }

    @ParameterizedTest
    @MethodSource("poolFailures")
    void shouldReportAFailedDrainAndStopInsteadOfFailingTheProducer(Throwable failure) throws Exception {
        var tasks = new ArrayList<Runnable>();
        var reject = new AtomicBoolean();
        Executor pool = task -> {
            if( reject.get() )
                throw raise(failure);
            tasks.add(task);
        };
        var input = new DataflowQueue<Object>();
        var recorder = new Recorder();
        var op = new DataflowProcessor(pool, List.of(input), List.of(), List.of(recorder), 1, Closures.action(x -> {}));

        op.start();
        // subscribe the input
        runAll(tasks);
        reject.set(true);
        input.bind(1);

        op.join();
        assertThat(op.isActive()).isFalse();
        assertThat(recorder.errors).containsExactly(failure);
        assertThat(recorder.stops.get()).isEqualTo(1);
        assertThat(recorder.before.get()).isZero();
    }

    @ParameterizedTest
    @MethodSource("poolFailures")
    void shouldReportAFailedRescheduleAndStopInsteadOfStrandingTheWork(Throwable failure) throws Exception {
        var tasks = new ArrayList<Runnable>();
        var reject = new AtomicBoolean();
        Executor pool = task -> {
            if( reject.get() )
                throw raise(failure);
            tasks.add(task);
        };
        var input = new DataflowQueue<Object>();
        var recorder = new Recorder();
        var fired = new AtomicInteger();
        var op = new DataflowProcessor(pool, List.of(input), List.of(), List.of(recorder), 1, Closures.action(x -> fired.incrementAndGet()));

        // more values than a batch, so that the drain reschedules
        for( int i = 0; i < 1000; i++ )
            input.bind(i);
        op.start();
        reject.set(true);
        runAll(tasks);

        op.join();
        assertThat(op.isActive()).isFalse();
        assertThat(recorder.errors).containsExactly(failure);
        assertThat(recorder.stops.get()).isEqualTo(1);
        assertThat(fired.get()).isBetween(1, 999);
    }

    @ParameterizedTest
    @MethodSource("poolFailures")
    void shouldStopWhenTerminatedOnAPoolThatFailsToRunTheDrain(Throwable failure) throws Exception {
        var tasks = new ArrayList<Runnable>();
        var reject = new AtomicBoolean();
        Executor pool = task -> {
            if( reject.get() )
                throw raise(failure);
            tasks.add(task);
        };
        var recorder = new Recorder();
        var op = new DataflowProcessor(pool, List.of(new DataflowQueue<Object>()), List.of(), List.of(recorder), 1, Closures.action(x -> {}));

        op.start();
        runAll(tasks);
        reject.set(true);
        op.terminate();

        op.join();
        assertThat(recorder.errors).containsExactly(failure);
        assertThat(recorder.stops.get()).isEqualTo(1);
    }

    @ParameterizedTest
    @MethodSource("poolFailures")
    void shouldPairBeforeRunAndAfterRunAndStopWhenAForkFails(Throwable failure) throws Exception {
        var tasks = new ArrayList<Runnable>();
        // run the drains, fail the forks
        Executor pool = task -> {
            if( !(task instanceof MailboxEngine) )
                throw raise(failure);
            tasks.add(task);
        };
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var recorder = new Recorder();
        var started = new AtomicInteger();
        var op = new DataflowProcessor(pool, List.of(input), List.of(output), List.of(recorder), 2, Closures.action(x -> started.incrementAndGet()));

        op.start();
        input.bind(1);
        input.bind(STOP);
        runAll(tasks);

        op.join();
        assertThat(recorder.errors).containsExactly(failure);
        assertThat(recorder.before.get()).isEqualTo(1);
        assertThat(recorder.after.get()).isEqualTo(1);
        assertThat(recorder.stops.get()).isEqualTo(1);
        assertThat(started.get()).isZero();
    }

    @Test
    void shouldFireAfterStopAndReportTheErrorWhenThePoisonPillCannotBeForwardedOnTerminate() throws Exception {
        // run the drains on demand and never the forks, so that the pill waits for a body in flight
        var drains = new ArrayList<Runnable>();
        Executor pool = task -> {
            if( task instanceof MailboxEngine )
                drains.add(task);
        };
        var input = new DataflowQueue<Object>();
        var output = new DataflowVariable<Object>();
        output.bind("x");
        var recorder = new Recorder();
        var op = new DataflowProcessor(pool, List.of(input), List.of(output), List.of(recorder), 2, Closures.action(x -> {}));

        op.start();
        input.bind(1);
        input.bind(STOP);
        runAll(drains);
        assertThat(recorder.before.get()).isEqualTo(1);
        assertThat(recorder.stops.get()).isZero();

        op.terminate();
        runAll(drains);

        op.join();
        assertThat(recorder.stops.get()).isEqualTo(1);
        assertThat(recorder.errors).singleElement().isInstanceOf(IllegalStateException.class);
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
    @DisabledIfEnvironmentVariable(named = "NXF_DATAFLOW_POOL_SIZE", matches = "\\s*1\\s*", disabledReason = "needs more than one dataflow thread")
    void shouldPairBeforeRunAndAfterRunWhenAParkedInputSetIsDroppedOnStop() throws Exception {
        var input = new DataflowQueue<Object>();
        var release = new CountDownLatch(1);
        var started = new AtomicInteger();
        var before = new AtomicInteger();
        var after = new AtomicInteger();
        var stopped = new StopCounter();
        var listener = new DataflowEventListener() {
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

        var op = new DataflowProcessor(List.of(input), List.of(), List.of(listener, stopped), 2, Closures.action(x -> {
            started.incrementAndGet();
            release.await();
        }));
        op.start();
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
        for( int round = 0; round < 100; round++ ) {
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
        // spread the terminate over the binds
        for( int spin = round * 20 % 2000; spin > 0; spin-- )
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
