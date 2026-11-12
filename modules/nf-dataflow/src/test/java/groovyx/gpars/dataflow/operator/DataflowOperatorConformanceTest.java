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
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

import groovy.lang.Closure;
import groovyx.gpars.Closures;
import groovyx.gpars.dataflow.Dataflow;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.Closures.function;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Pins the GPars operator semantics that Nextflow relies on: {@code DataflowHelper.newOperator}
 * ({@code Dataflow.operator}) and {@code TaskProcessor} ({@code new DataflowOperator(group, params, code).start()}).
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class DataflowOperatorConformanceTest {

    static final long WAIT = 5;

    static final PoisonPill STOP = PoisonPill.getInstance();

    /**
     * Records the listener callbacks that Nextflow uses, in the order they fire.
     */
    static class Recorder extends DataflowEventAdapter {
        final List<String> events = new CopyOnWriteArrayList<>();
        final DataflowVariable<Object> stopped = new DataflowVariable<>();
        final AtomicInteger stops = new AtomicInteger();
        final boolean terminateOnError;

        Recorder() { this(true); }

        Recorder(boolean terminateOnError) { this.terminateOnError = terminateOnError; }

        @Override
        public Object messageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
            events.add("arrived:" + index + ":" + message);
            return message;
        }

        @Override
        public Object controlMessageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
            events.add("control:" + index + ":" + (message == STOP ? "STOP" : message));
            return message;
        }

        @Override
        public List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
            events.add("beforeRun:" + messages);
            return messages;
        }

        @Override
        public void afterRun(DataflowProcessor processor, List<Object> messages) {
            events.add("afterRun:" + messages);
        }

        @Override
        public boolean onException(DataflowProcessor processor, Throwable e) {
            events.add("error:" + e.getMessage());
            return terminateOnError;
        }

        @Override
        public void afterStop(DataflowProcessor processor) {
            events.add("afterStop");
            stops.incrementAndGet();
            stopped.bind(processor);
        }
    }

    @FunctionalInterface
    interface Body1 { void call(DataflowProcessor op, Object x) throws Throwable; }

    @FunctionalInterface
    interface Body2 { void call(DataflowProcessor op, Object x, Object y) throws Throwable; }

    /**
     * An operator body, which receives the operator (the closure delegate) like a Groovy body would.
     */
    static Closure<Object> body(Body1 body) {
        return Closures.closure(1, (op, args) -> { body.call((DataflowProcessor) op, args[0]); return null; });
    }

    static Closure<Object> body(Body2 body) {
        return Closures.closure(2, (op, args) -> { body.call((DataflowProcessor) op, args[0], args[1]); return null; });
    }

    /**
     * A closure declaring its arity, like Nextflow's InvokeTaskAdapter and ForwardClosure,
     * which records its delegate and arguments.
     */
    static Closure<Object> recording(int arity, DataflowQueue<Object> calls) {
        return Closures.closure(arity, (delegate, args) -> {
            calls.bind(new Call(delegate, Arrays.asList(args)));
            return null;
        });
    }

    record Call(Object delegate, List<Object> args) {}

    static List<Object> read(DataflowReadChannel<Object> ch, int n) throws InterruptedException {
        var result = new ArrayList<Object>(n);
        for( int i = 0; i < n; i++ )
            result.add(ch.getVal(WAIT, TimeUnit.SECONDS));
        return result;
    }

    static boolean isEmpty(DataflowReadChannel<Object> ch) throws InterruptedException {
        return ch.getVal(100, TimeUnit.MILLISECONDS) == null;
    }

    static DataflowVariable<Object> joinAsync(DataflowProcessor op) {
        final var result = new DataflowVariable<Object>();
        start(() -> { op.join(); result.bind(true); });
        return result;
    }

    static List<Object> range(int n) {
        return new ArrayList<>(IntStream.range(0, n).boxed().toList());
    }

    static void send(DataflowQueue<Object> ch, Object... values) {
        for( var value : values )
            ch.bind(value);
    }

    // -- creation and outputs

    @Test
    void shouldCreateAStartedOperatorWithDataflowOperator() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output)), body((self, x) -> self.bindOutput((Integer) x * 2)));
        send(input, 1, 2, 3);
        assertThat(op).isInstanceOf(DataflowOperator.class);
        assertThat(read(output, 3)).containsExactly(2, 4, 6);
    }

    @Test
    void shouldNotFireBeforeStartWhenCreatedWithTheConstructor() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var group = Dataflow.retrieveCurrentDFPGroup();

        var op = new DataflowOperator(group, Map.of("inputs", List.of(input), "outputs", List.of(output)), body((self, x) -> self.bindOutput((Integer) x + 1)));
        send(input, 1, 2);
        assertThat(group).isNotNull();
        assertThat(isEmpty(output)).isTrue();

        var started = op.start();
        assertThat(started).isSameAs(op);
        assertThat(read(output, 2)).containsExactly(2, 3);
    }

    @Test
    void shouldUseTheProcessorAsTheBodyDelegate() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var seen = new DataflowVariable<List<Object>>();

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output)), body((self, x) -> {
            seen.bind(List.of(self, self.getOutput(), self.getOutputs()));
            self.bindOutput(x);
        }));
        input.bind("a");
        var result = seen.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(result.get(0)).isSameAs(op);
        assertThat(result.get(1)).isSameAs(output);
        assertThat(result.get(2)).isEqualTo(List.of(output));
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("a");
    }

    @Test
    void shouldSetTheProcessorAsDelegateOfAClosureSubclassAndPassOneArgumentPerInput() throws Exception {
        var a = new DataflowQueue<Object>();
        var b = new DataflowQueue<Object>();
        var calls = new DataflowQueue<Object>();

        var op = new DataflowOperator(Dataflow.retrieveCurrentDFPGroup(), Map.of("inputs", List.of(a, b)), recording(2, calls)).start();
        a.bind(1);
        b.bind("x");
        var call = (Call) calls.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(call.delegate()).isSameAs(op);
        assertThat(call.args()).containsExactly(1, "x");
    }

    @Test
    void shouldExposeTheOutputsAndBindThemByIndexOrAllAtOnce() throws Exception {
        var input = new DataflowQueue<Object>();
        var o1 = new DataflowQueue<Object>();
        var o2 = new DataflowQueue<Object>();

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(o1, o2)), body((self, x) -> {
            self.bindOutput(x);
            self.bindOutput(1, -(Integer) x);
            self.bindAllOutputs("all");
            self.bindAllOutputValues("v1", "v2");
            self.bindAllOutputsAtomically("atomic");
        }));
        input.bind(7);
        assertThat(op.getOutputs()).containsExactly(o1, o2);
        assertThat(op.getOutput()).isSameAs(o1);
        assertThat(op.getOutputs(1)).isSameAs(o2);
        assertThat(read(o1, 4)).containsExactly(7, "all", "v1", "atomic");
        assertThat(read(o2, 4)).containsExactly(-7, "all", "v2", "atomic");
    }

    @Test
    void shouldRejectInvalidParameters() {
        var group = Dataflow.retrieveCurrentDFPGroup();
        var q = new DataflowQueue<Object>();

        assertThatThrownBy(() -> new DataflowOperator(group, Map.of("inputs", List.of()), body((self, x) -> {})))
            .isInstanceOf(IllegalArgumentException.class);

        assertThatThrownBy(() -> new DataflowOperator(group, Map.of("inputs", List.of(q)), body((self, x, y) -> {})))
            .isInstanceOf(IllegalArgumentException.class);

        assertThatThrownBy(() -> new DataflowOperator(group, Map.of("inputs", List.of(q), "maxForks", 0), body((self, x) -> {})))
            .isInstanceOf(IllegalArgumentException.class);
    }

    // -- firing

    @Test
    void shouldFireOncePerCompleteSetOfInputsWithArgumentsInInputOrder() throws Exception {
        var a = new DataflowQueue<Object>();
        var b = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();

        Dataflow.operator(Map.of("inputs", List.of(a, b), "outputs", List.of(output)), body((self, x, y) -> self.bindOutput(List.of(x, y))));
        send(b, "p", "q");
        send(a, 1, 2, 3);
        assertThat(read(output, 2)).containsExactly(List.of(1, "p"), List.of(2, "q"));
        assertThat(isEmpty(output)).isTrue();

        // the third item stays unmatched until b provides a partner
        b.bind("r");
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(List.of(3, "r"));
    }

    @Test
    void shouldRunTheBodySeriallyAndInOrderWithoutMaxForks() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var running = new AtomicInteger();
        var maxRunning = new AtomicInteger();
        var n = 50;

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output)), body((self, x) -> {
            maxRunning.accumulateAndGet(running.incrementAndGet(), Math::max);
            // widen the window in which overlapping bodies would be observed
            Thread.sleep(1);
            running.decrementAndGet();
            self.bindOutput(x);
        }));
        for( int i = 0; i < n; i++ )
            input.bind(i);
        assertThat(read(output, n)).isEqualTo(range(n));
        assertThat(maxRunning.get()).isEqualTo(1);
    }

    @Test
    void shouldRunUpToMaxForksBodiesConcurrentlyAndNoMore() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var running = new AtomicInteger();
        var maxRunning = new AtomicInteger();
        var allIn = new CountDownLatch(3);
        var n = 30;

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "maxForks", 3), body((self, x) -> {
            maxRunning.accumulateAndGet(running.incrementAndGet(), Math::max);
            // the first three bodies only complete if they run at the same time
            allIn.countDown();
            var together = allIn.await(WAIT, TimeUnit.SECONDS);
            // widen the window in which a fourth overlapping body would be observed
            Thread.sleep(2);
            running.decrementAndGet();
            self.bindOutput(together);
        }));
        for( int i = 0; i < n; i++ )
            input.bind(i);
        assertThat(read(output, n)).containsOnly(true);
        assertThat(maxRunning.get()).isEqualTo(3);
    }

    @Test
    void shouldCallBeforeRunSeriallyInArrivalOrderWithMaxForks() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var order = new CopyOnWriteArrayList<Object>();
        var listener = new DataflowEventAdapter() {
            @Override
            public List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
                order.add(messages.get(0));
                return messages;
            }
        };
        var n = 200;

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "maxForks", 4, "listeners", List.of(listener)), body((self, x) -> self.bindOutput(x)));
        for( int i = 0; i < n; i++ )
            input.bind(i);
        assertThat(read(output, n)).containsExactlyInAnyOrderElementsOf(range(n));
        assertThat(order).isEqualTo(range(n));
    }

    @Test
    void shouldPassTheResultOfMessageArrivedAndBeforeRunToTheBody() throws Exception {
        var input = new DataflowQueue<Object>();
        var calls = new DataflowQueue<Object>();
        var listener = new DataflowEventAdapter() {
            @Override
            public Object messageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
                return (Integer) message * 10;
            }

            @Override
            public List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
                // like TaskProcessorInterceptor: replace the arguments with [params, messages]
                return List.of("params", messages);
            }
        };

        new DataflowOperator(Dataflow.retrieveCurrentDFPGroup(), Map.of("inputs", List.of(input), "listeners", List.of(listener)), recording(1, calls)).start();
        input.bind(1);
        assertThat(((Call) calls.getVal(WAIT, TimeUnit.SECONDS)).args()).containsExactly("params", List.of(10));
    }

    @Test
    void shouldReReadAValueInputOnEveryFiringWithoutConsumingIt() throws Exception {
        var items = new DataflowQueue<Object>();
        var value = new DataflowVariable<Object>();
        var output = new DataflowQueue<Object>();

        Dataflow.operator(Map.of("inputs", List.of(items, value), "outputs", List.of(output)), body((self, x, v) -> self.bindOutput((Integer) x + (Integer) v)));
        send(items, 1, 2, 3);
        assertThat(isEmpty(output)).isTrue();

        value.bind(10);
        assertThat(read(output, 3)).containsExactly(11, 12, 13);

        items.bind(4);
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(14);
        assertThat(value.getVal()).isEqualTo(10);
    }

    @Test
    void shouldFireOnceOnAValueInputWhenTheBodyTerminatesTheOperator() throws Exception {
        // the DataflowHelper.subscribeImpl pattern for value channels
        var value = new DataflowVariable<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder();
        var count = new AtomicInteger();

        var op = Dataflow.operator(Map.of("inputs", List.of(value), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x) -> {
            count.incrementAndGet();
            self.bindOutput(x);
            self.terminate();
        }));
        value.bind("v");
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("v");
        assertThat(joinAsync(op).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
        assertThat(isEmpty(output)).isTrue();
        assertThat(count.get()).isEqualTo(1);
        assertThat(rec.stops.get()).isEqualTo(1);
    }

    @Test
    void shouldReadFromABroadcastReadChannel() throws Exception {
        var bc = new DataflowBroadcast<Object>();
        var o1 = new DataflowQueue<Object>();
        var o2 = new DataflowQueue<Object>();

        Dataflow.operator(Map.of("inputs", List.of(bc.createReadChannel()), "outputs", List.of(o1)), body((self, x) -> self.bindOutput(x)));
        Dataflow.operator(Map.of("inputs", List.of(bc.createReadChannel()), "outputs", List.of(o2)), body((self, x) -> self.bindOutput((Integer) x * 2)));
        bc.leftShift(1).leftShift(2).leftShift(STOP);
        assertThat(read(o1, 3)).containsExactly(1, 2, STOP);
        assertThat(read(o2, 3)).containsExactly(2, 4, STOP);
    }

    @Test
    void shouldCopyItemsWithAChainWithClosure() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output)), new ChainWithClosure<>(function(x -> (Integer) x * 3)));
        send(input, 1, 2, STOP);
        assertThat(read(output, 3)).containsExactly(3, 6, STOP);
    }

    // -- poison pill

    @Test
    void shouldForwardAPoisonPillToEveryOutputAfterThePrecedingItemsAndStop() throws Exception {
        var input = new DataflowQueue<Object>();
        var o1 = new DataflowQueue<Object>();
        var o2 = new DataflowQueue<Object>();
        var rec = new Recorder();

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(o1, o2), "listeners", List.of(rec)), body((self, x) -> self.bindAllOutputs(x)));
        send(input, 1, 2, STOP);
        assertThat(read(o1, 3)).containsExactly(1, 2, STOP);
        assertThat(read(o2, 3)).containsExactly(1, 2, STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(joinAsync(op).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
    }

    @Test
    void shouldStopWhenAPoisonPillArrivesOnAnyInput() throws Exception {
        var a = new DataflowQueue<Object>();
        var b = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder();

        var op = Dataflow.operator(Map.of("inputs", List.of(a, b), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x, y) -> self.bindOutput((Integer) x + (Integer) y)));
        send(a, 1, 2, 3);
        send(b, 10, STOP);
        assertThat(read(output, 2)).containsExactly(11, STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(rec.events).contains("control:1:STOP");
        assertThat(isEmpty(output)).isTrue();
    }

    @Test
    void shouldStopOnAPoisonPillInAQueueInputCombinedWithAValueInput() throws Exception {
        var items = new DataflowQueue<Object>();
        var value = new DataflowVariable<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder();
        value.bind(100);

        var op = Dataflow.operator(Map.of("inputs", List.of(items, value), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x, v) -> self.bindOutput((Integer) x + (Integer) v)));
        send(items, 1, STOP);
        assertThat(read(output, 2)).containsExactly(101, STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
    }

    @Test
    void shouldForwardASinglePoisonPillPerOutputWhenSeveralInputsArePoisoned() throws Exception {
        var a = new DataflowQueue<Object>();
        var b = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder();

        Dataflow.operator(Map.of("inputs", List.of(a, b), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x, y) -> self.bindOutput(x)));
        a.bind(STOP);
        b.bind(STOP);
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isNotNull();
        assertThat(isEmpty(output)).isTrue();
        assertThat(rec.stops.get()).isEqualTo(1);
    }

    @Test
    void shouldWaitForInFlightBodiesBeforeForwardingThePoisonPillWithMaxForks() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var started = new CountDownLatch(4);
        var release = new CountDownLatch(1);
        var rec = new Recorder();

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "maxForks", 4, "listeners", List.of(rec)), body((self, x) -> {
            started.countDown();
            release.await(WAIT, TimeUnit.SECONDS);
            self.bindOutput(x);
        }));
        send(input, 1, 2, 3, 4, STOP);
        assertThat(started.await(WAIT, TimeUnit.SECONDS)).isTrue();
        assertThat(isEmpty(output)).isTrue();
        assertThat(rec.stopped.getVal(50, TimeUnit.MILLISECONDS)).isNull();

        release.countDown();
        var result = read(output, 5);
        assertThat(result.subList(0, 4)).containsExactlyInAnyOrder(1, 2, 3, 4);
        assertThat(result.get(4)).isSameAs(STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isNotNull();
        assertThat(rec.events).filteredOn(it -> it.startsWith("afterRun")).hasSize(4);
    }

    // -- listeners and termination

    @Test
    void shouldFireListenerCallbacksInOrder() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder();

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x) -> {
            rec.events.add("body:" + x);
            self.bindOutput(x);
        }));
        send(input, 1, 2, STOP);
        rec.stopped.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(rec.events).containsExactly(
            "arrived:0:1", "beforeRun:[1]", "body:1", "afterRun:[1]",
            "arrived:0:2", "beforeRun:[2]", "body:2", "afterRun:[2]",
            "control:0:STOP",
            "afterStop"
        );
        assertThat(read(output, 3)).containsExactly(1, 2, STOP);
    }

    @Test
    void shouldAddAListenerAfterCreation() throws Exception {
        var input = new DataflowQueue<Object>();
        var rec = new Recorder();
        var group = Dataflow.retrieveCurrentDFPGroup();

        var op = new DataflowOperator(group, Map.of("inputs", List.of(input)), body((self, x) -> {}));
        op.addDataflowEventListener(rec);
        op.start();
        send(input, 1, STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(rec.events).containsExactly("arrived:0:1", "beforeRun:[1]", "afterRun:[1]", "control:0:STOP", "afterStop");
    }

    @Test
    void shouldKeepRunningWhenOnExceptionReturnsFalse() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder(false);

        Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x) -> {
            if( x.equals(2) ) throw new IllegalStateException("bad " + x);
            self.bindOutput(x);
        }));
        send(input, 1, 2, 3, STOP);
        assertThat(read(output, 3)).containsExactly(1, 3, STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isNotNull();
        // afterRun is called also for the failed firing
        assertThat(rec.events).contains("error:bad 2");
        assertThat(rec.events.indexOf("error:bad 2")).isLessThan(rec.events.indexOf("afterRun:[2]"));
    }

    @Test
    void shouldTerminateWhenOnExceptionReturnsTrue() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder(true);

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x) -> {
            if( x.equals(2) ) throw new IllegalStateException("bad " + x);
            self.bindOutput(x);
        }));
        send(input, 1, 2, 3);
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(1);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(joinAsync(op).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
        assertThat(rec.events).contains("error:bad 2");
        // an error does not forward a poison pill and later items are not processed
        assertThat(isEmpty(output)).isTrue();
        assertThat(rec.stops.get()).isEqualTo(1);
    }

    @Test
    void shouldTerminateOnAnExceptionWhenThereAreNoListeners() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output)), body((self, x) -> {
            if( x.equals(2) ) throw new IllegalStateException("expected test error");
            self.bindOutput(x);
        }));
        send(input, 1, 2, 3);
        assertThat(output.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(1);
        assertThat(joinAsync(op).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
        assertThat(isEmpty(output)).isTrue();
    }

    @Test
    void shouldFireAfterStopAndReleaseJoinOnAnExternalTerminateWhileIdle() throws Exception {
        var input = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var rec = new Recorder();
        var op = Dataflow.operator(Map.of("inputs", List.of(input), "outputs", List.of(output), "listeners", List.of(rec)), body((self, x) -> self.bindOutput(x)));
        input.bind(1);
        output.getVal(WAIT, TimeUnit.SECONDS);
        var joined = joinAsync(op);

        assertThat(joined.getVal(100, TimeUnit.MILLISECONDS)).isNull();

        op.terminate();
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(joined.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);

        // no pill is forwarded and later items are not processed
        input.bind(2);
        assertThat(isEmpty(output)).isTrue();
    }

    @Test
    @Tag("gpars-incidental")
    void shouldInterruptARunningBodyOnAnExternalTerminate() throws Exception {
        var input = new DataflowQueue<Object>();
        var rec = new Recorder();
        var started = new CountDownLatch(1);
        var interrupted = new DataflowVariable<Object>();
        var op = Dataflow.operator(Map.of("inputs", List.of(input), "listeners", List.of(rec)), body((self, x) -> {
            started.countDown();
            try {
                Thread.sleep(30_000);
                interrupted.bind(false);
            }
            catch( InterruptedException e ) {
                interrupted.bind(true);
            }
        }));
        input.bind(1);
        assertThat(started.await(WAIT, TimeUnit.SECONDS)).isTrue();

        op.terminate();
        assertThat(interrupted.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(joinAsync(op).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
    }

    @Test
    void shouldFireAfterStopOnlyOnce() throws Exception {
        var input = new DataflowQueue<Object>();
        var rec = new Recorder();
        var op = Dataflow.operator(Map.of("inputs", List.of(input), "listeners", List.of(rec)), body((self, x) -> {}));

        input.bind(STOP);
        rec.stopped.getVal(WAIT, TimeUnit.SECONDS);
        op.terminate();
        op.terminate();
        // a short wait to catch a duplicate callback: it can only miss a bug, never cause a false failure
        Thread.sleep(50);
        assertThat(rec.stops.get()).isEqualTo(1);
    }

    @Test
    void shouldBlockJoinUntilTheOperatorStops() throws Exception {
        var input = new DataflowQueue<Object>();
        var op = Dataflow.operator(Map.of("inputs", List.of(input)), body((self, x) -> {}));
        var joined = joinAsync(op);

        assertThat(joined.getVal(100, TimeUnit.MILLISECONDS)).isNull();

        input.bind(STOP);
        assertThat(joined.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(true);
    }

    @Test
    void shouldStopAnOperatorWithNoOutputsOnAPoisonPill() throws Exception {
        var input = new DataflowQueue<Object>();
        var rec = new Recorder();
        var seen = new CopyOnWriteArrayList<Object>();

        var op = Dataflow.operator(Map.of("inputs", List.of(input), "listeners", List.of(rec)), body((self, x) -> seen.add(x)));
        send(input, 1, 2, STOP);
        assertThat(rec.stopped.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(op);
        assertThat(seen).containsExactly(1, 2);
        assertThat(op.getOutputs()).isEmpty();
        assertThat(op.getOutput()).isNull();
    }

}
