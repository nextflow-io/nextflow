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
package groovyx.gpars.dataflow;

import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Supplier;

import groovy.lang.Binding;
import groovy.lang.GroovyShell;
import groovy.lang.MissingMethodException;
import groovy.lang.MissingPropertyException;
import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.dataflow.impl.Waiter;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.TestThreads.await;
import static groovyx.gpars.TestThreads.awaitBlocked;
import static groovyx.gpars.TestThreads.poll;
import static groovyx.gpars.TestThreads.start;
import static groovyx.gpars.dataflow.impl.ReadCallback.NONE;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * The lightweight channel internals: non-blocking reads, null values, lost-value safety of interrupted
 * reads, failing callbacks, garbage collection of broadcast nodes and readers, and the dynamic API of
 * value channels.
 */
@Timeout(value = 20, unit = TimeUnit.SECONDS)
class ChannelInternalsTest {

    static final long WAIT = 5;

    static final class Recorder implements ReadCallback {
        final List<Object> items = new CopyOnWriteArrayList<>();

        @Override
        public void onItem(Object item) {
            items.add(item);
        }
    }

    /**
     * A channel whose subscriptions cannot be cancelled, so that a read that gives up waits for the hand-over
     */
    static final class HandingOverChannel implements DataflowReadChannel<Object> {
        volatile ReadCallback callback;
        final AtomicInteger unsubscribes = new AtomicInteger();

        @Override
        public Object takeOrSubscribe(ReadCallback callback) {
            this.callback = callback;
            return NONE;
        }

        @Override
        public boolean unsubscribe(ReadCallback callback) {
            unsubscribes.incrementAndGet();
            return false;
        }

        @Override
        public Object getVal() throws InterruptedException {
            return Waiter.read(this);
        }

        @Override
        public void into(DataflowWriteChannel<Object> target) {
            throw new UnsupportedOperationException();
        }

        @Override
        public boolean isBound() {
            return false;
        }

        @Override
        public int length() {
            return 0;
        }
    }

    @Test
    void shouldTakeABufferedValueOrSubscribeOnceOnAQueue() {
        var q = new DataflowQueue<Object>();
        var first = new Recorder();
        var second = new Recorder();
        q.bind("a");

        assertThat(q.takeOrSubscribe(first)).isEqualTo("a");
        assertThat(q.takeOrSubscribe(first)).isSameAs(NONE);
        assertThat(q.takeOrSubscribe(second)).isSameAs(NONE);

        // each value goes to one subscriber, in subscription order, and is not buffered
        q.bind("b");
        q.bind("c");
        q.bind("d");
        assertThat(first.items).containsExactly("b");
        assertThat(second.items).containsExactly("c");
        assertThat(q.length()).isEqualTo(1);
    }

    @Test
    void shouldNotHandAValueToAnUnsubscribedReader() {
        var q = new DataflowQueue<Object>();
        var reader = new Recorder();
        assertThat(q.takeOrSubscribe(reader)).isSameAs(NONE);

        assertThat(q.unsubscribe(reader)).isTrue();
        assertThat(q.unsubscribe(reader)).isFalse();
        q.bind("a");
        assertThat(reader.items).isEmpty();
        assertThat(q.length()).isEqualTo(1);
    }

    @Test
    void shouldReadAValueWithoutConsumingIt() {
        var v = new DataflowVariable<Object>();
        var early = new Recorder();
        assertThat(v.takeOrSubscribe(early)).isSameAs(NONE);

        v.bind("x");
        assertThat(early.items).containsExactly("x");
        assertThat(v.takeOrSubscribe(new Recorder())).isEqualTo("x");
        assertThat(v.takeOrSubscribe(new Recorder())).isEqualTo("x");
    }

    @Test
    void shouldPassNullValuesThroughEveryChannel() throws Exception {
        var q = new DataflowQueue<Object>();
        q.bind(null);
        q.bind("x");
        assertThat(q.length()).isEqualTo(2);
        assertThat(q.getVal()).isNull();
        assertThat(q.getVal()).isEqualTo("x");

        var b = new DataflowBroadcast<Object>();
        var r = b.createReadChannel();
        b.bind(null);
        assertThat(r.isBound()).isTrue();
        assertThat(r.getVal()).isNull();
        assertThat(r.isBound()).isFalse();
    }

    @Test
    void shouldNotLoseAValueWhenAReadIsInterrupted() throws Exception {
        var q = new DataflowQueue<Object>();
        var error = new DataflowVariable<Object>();
        var thread = start(() -> {
            try {
                q.getVal();
            }
            catch( InterruptedException e ) {
                error.bind(e);
            }
        });
        awaitBlocked(thread);
        thread.interrupt();
        assertThat(poll(error, WAIT, TimeUnit.SECONDS)).isInstanceOf(InterruptedException.class);

        q.bind("a");
        assertThat(poll(q, WAIT, TimeUnit.SECONDS)).isEqualTo("a");
    }

    @Test
    void shouldCallEveryWaiterOfAValueWhenSomeFail() {
        var v = new DataflowVariable<Object>();
        var first = new IllegalStateException("first");
        var second = new IllegalStateException("second");
        var recorder = new Recorder();
        v.takeOrSubscribe(item -> { throw first; });
        v.takeOrSubscribe(item -> { throw second; });
        v.takeOrSubscribe(recorder);

        assertThatThrownBy(() -> v.bind("x")).isSameAs(first).hasSuppressedException(second);
        assertThat(recorder.items).containsExactly("x");
    }

    @Test
    void shouldWakeABlockedReadWhenAnEarlierCallbackFails() throws Exception {
        var v = new DataflowVariable<Object>();
        v.takeOrSubscribe(item -> { throw new IllegalStateException("boom"); });
        var result = new DataflowVariable<Object>();
        var reader = start(() -> result.bind(v.getVal()));
        awaitBlocked(reader);

        assertThatThrownBy(() -> v.bind("x")).hasMessage("boom");
        assertThat(poll(result, WAIT, TimeUnit.SECONDS)).isEqualTo("x");
    }

    @Test
    void shouldWakeEveryBroadcastReaderWhenSomeCallbackFails() {
        var b = new DataflowBroadcast<Object>();
        var failing = b.createReadChannel();
        var other = b.createReadChannel();
        var next = new Recorder();
        var recorder = new Recorder();
        failing.takeOrSubscribe(item -> { throw new IllegalStateException("boom"); });
        failing.takeOrSubscribe(next);
        other.takeOrSubscribe(recorder);

        assertThatThrownBy(() -> b.bind(1)).hasMessage("boom");
        assertThat(recorder.items).containsExactly(1);
        // the failing reader still waits for the next value
        b.bind(2);
        assertThat(next.items).containsExactly(2);
    }

    @Test
    void shouldParkUntilTheValueIsHandedOverWhenAReadIsInterrupted() throws Exception {
        var channel = new HandingOverChannel();
        var result = new DataflowVariable<Object>();
        var reader = start(() -> result.bind(List.of(channel.getVal(), Thread.currentThread().isInterrupted())));
        await(() -> channel.callback != null && reader.getState() == Thread.State.WAITING);

        reader.interrupt();
        await(() -> channel.unsubscribes.get() == 1 && reader.getState() == Thread.State.WAITING);
        channel.callback.onItem("x");
        // the value is not lost and the interrupt is kept
        assertThat(poll(result, WAIT, TimeUnit.SECONDS)).isEqualTo(List.of("x", true));
    }

    @Test
    void shouldHandBufferedValuesToTheReadsOfABroadcastReaderInOrder() throws Exception {
        var b = new DataflowBroadcast<Object>();
        var r = (DataflowStreamReadAdapter<Object>) b.createReadChannel();
        var first = new Recorder();
        var second = new Recorder();

        assertThat(r.takeOrSubscribe(first)).isSameAs(NONE);
        assertThat(r.takeOrSubscribe(second)).isSameAs(NONE);
        b.bind(1);
        b.bind(2);
        b.bind(3);
        assertThat(first.items).containsExactly(1);
        assertThat(second.items).containsExactly(2);
        // the blocking and the asynchronous reads share the cursor
        assertThat(r.length()).isEqualTo(1);
        assertThat(r.getVal()).isEqualTo(3);
    }

    @Test
    void shouldIterateOverTheBufferedValuesWithoutConsumingThem() throws Exception {
        var q = new DataflowQueue<Object>();
        var b = new DataflowBroadcast<Object>();
        var r = (DataflowStreamReadAdapter<Object>) b.createReadChannel();
        for( int i = 0; i < 3; i++ ) {
            q.bind(i);
            b.bind(i);
        }
        r.getVal();

        assertThat(q.iterator()).toIterable().containsExactly(0, 1, 2);
        assertThat(r.iterator()).toIterable().containsExactly(1, 2);
        assertThat(q.length()).isEqualTo(3);
        assertThat(r.length()).isEqualTo(2);
    }

    @Test
    void shouldDeliverEveryValueOnceToCompetingQueueReaders() throws Exception {
        final int n = 20_000;
        var q = new DataflowQueue<Object>();
        var received = new CopyOnWriteArrayList<Object>();
        var readers = new ArrayList<Thread>();
        for( int i = 0; i < 4; i++ ) {
            readers.add(start(() -> {
                Object item;
                while( (item = q.getVal()) != null )
                    received.add(item);
            }));
        }
        var writers = new ArrayList<Thread>();
        for( int w = 0; w < 2; w++ ) {
            final int offset = w * n;
            writers.add(start(() -> { for( int i = 0; i < n; i++ ) q.bind(offset + i); }));
        }
        for( var writer : writers )
            writer.join();
        for( int i = 0; i < readers.size(); i++ )
            q.bind(null);
        for( var reader : readers )
            reader.join(WAIT * 1000);

        assertThat(received).hasSize(2 * n).doesNotHaveDuplicates();
    }

    @Test
    void shouldDeliverEveryValueInOrderToConcurrentBroadcastReaders() throws Exception {
        final int n = 20_000;
        var b = new DataflowBroadcast<Object>();
        var results = new ArrayList<List<Object>>();
        var threads = new ArrayList<Thread>();
        var ready = new CountDownLatch(4);
        for( int i = 0; i < 4; i++ ) {
            final var reader = b.createReadChannel();
            final var result = new ArrayList<Object>();
            results.add(result);
            threads.add(start(() -> {
                ready.countDown();
                for( int k = 0; k < n; k++ )
                    result.add(reader.getVal());
            }));
        }
        ready.await();
        for( int i = 0; i < n; i++ )
            b.bind(i);
        for( var thread : threads )
            thread.join(WAIT * 1000);

        for( var result : results ) {
            assertThat(result).hasSize(n);
            for( int i = 0; i < n; i++ )
                assertThat(result.get(i)).isEqualTo(i);
        }
    }

    @Test
    void shouldDeliverEveryValueOfEachWriterInOrderToReadersCreatedMidStream() throws Exception {
        final int writers = 3;
        final int n = 20_000;
        final int readers = 8;
        var b = new DataflowBroadcast<Object>();
        var bound = new AtomicInteger();
        var created = new CountDownLatch(1);
        var errors = new CopyOnWriteArrayList<String>();
        var threads = new ArrayList<Thread>();
        // each writer binds its end marker once every reader is created, so that every reader sees it
        for( int w = 0; w < writers; w++ ) {
            final int writer = w;
            threads.add(start(() -> {
                for( int i = 0; i < n; i++ ) {
                    b.bind(new int[] { writer, i });
                    bound.incrementAndGet();
                }
                created.await();
                b.bind(new int[] { writer, -1 });
            }));
        }
        for( int k = 0; k < readers; k++ ) {
            final int progress = k * writers * n / readers;
            await(() -> bound.get() >= progress);
            final var reader = (DataflowStreamReadAdapter<Object>) b.createReadChannel();
            final boolean async = k % 2 == 0;
            threads.add(start(() -> {
                // the last index read of each writer
                final int[] last = new int[writers];
                Arrays.fill(last, Integer.MIN_VALUE);
                int ended = 0;
                while( ended < writers ) {
                    final int length = reader.length(5);
                    if( length < 0 || length > 5 )
                        errors.add("length " + length);
                    Object item;
                    if( async ) {
                        var future = new CompletableFuture<Object>();
                        item = reader.takeOrSubscribe(future::complete);
                        if( item == NONE )
                            item = future.get(WAIT, TimeUnit.SECONDS);
                    }
                    else
                        item = poll(reader, WAIT, TimeUnit.SECONDS);
                    if( item == null ) {
                        errors.add("timeout");
                        return;
                    }
                    final int writer = ((int[]) item)[0];
                    final int index = ((int[]) item)[1];
                    if( index == -1 ) {
                        ended++;
                        if( last[writer] != Integer.MIN_VALUE && last[writer] != n - 1 )
                            errors.add("writer " + writer + " ended at " + last[writer]);
                    }
                    else if( last[writer] != Integer.MIN_VALUE && index != last[writer] + 1 )
                        errors.add("writer " + writer + " skipped from " + last[writer] + " to " + index);
                    else
                        last[writer] = index;
                }
            }));
        }
        created.countDown();
        for( var thread : threads )
            thread.join(2 * WAIT * 1000);

        assertThat(threads).noneMatch(Thread::isAlive);
        assertThat(errors).isEmpty();
    }

    @Test
    void shouldNotRetainTheValuesOfABroadcastWithoutReaders() {
        var b = new DataflowBroadcast<Object>();
        var item = bindWeak(b);
        b.bind("last");

        assertCollected(item);
    }

    @Test
    void shouldReleaseTheValuesReadByEveryReader() throws Exception {
        var b = new DataflowBroadcast<Object>();
        var r1 = b.createReadChannel();
        var r2 = b.createReadChannel();
        var item = bindWeak(b);
        b.bind("last");

        // a reader holds the last node it read
        r1.getVal();
        r1.getVal();
        assertRetained(item);
        r2.getVal();
        r2.getVal();
        assertCollected(item);
    }

    @Test
    void shouldCollectAnUnreferencedReaderAndTheValuesItDidNotRead() {
        var b = new DataflowBroadcast<Object>();
        var reader = new WeakReference<>(b.createReadChannel());
        var item = bindWeak(b);
        b.bind("last");

        assertCollected(reader::get);
        assertCollected(item);
    }

    @Test
    void shouldCollectAReaderThatWasWaitingOnceTheNextValueIsBound() {
        var b = new DataflowBroadcast<Object>();
        var reader = b.createReadChannel();
        assertThat(reader.takeOrSubscribe(new Recorder())).isSameAs(NONE);
        var ref = new WeakReference<>(reader);
        reader = null;

        b.bind("x");
        assertCollected(ref::get);
    }

    @Test
    void shouldFailOnAnUnknownPropertyOrMethodOfAValueChannel() throws Exception {
        var v = new DataflowVariable<Object>();
        var binding = new Binding();
        binding.setVariable("v", v);
        var shell = new GroovyShell(binding);

        assertThat(shell.evaluate("v.bound")).isEqualTo(false);
        v.bind("hello");
        assertThat(shell.evaluate("v.val")).isEqualTo("hello");
        assertThat(shell.evaluate("v.getVal()")).isEqualTo("hello");
        assertThatThrownBy(() -> shell.evaluate("v.toUpperCase()"))
            .isInstanceOf(MissingMethodException.class)
            .hasMessageContaining("toUpperCase");
        assertThatThrownBy(() -> shell.evaluate("v.bytes"))
            .isInstanceOf(MissingPropertyException.class)
            .hasMessageContaining("bytes");
    }

    @Test
    void shouldCopyAValueWithInto() throws Exception {
        var v = new DataflowVariable<Object>();
        var target = new DataflowQueue<Object>();
        v.into(target);
        v.bind("x");
        assertThat(poll(target, WAIT, TimeUnit.SECONDS)).isEqualTo("x");
    }

    private static Supplier<Object> bindWeak(DataflowWriteChannel<Object> channel) {
        var item = new Object();
        channel.bind(item);
        var ref = new WeakReference<>(item);
        return ref::get;
    }

    private static void assertCollected(Supplier<?> ref) {
        for( int i = 0; i < 50 && ref.get() != null; i++ ) {
            System.gc();
            try {
                Thread.sleep(10);
            }
            catch( InterruptedException e ) {
                throw new AssertionError(e);
            }
        }
        assertThat(ref.get()).isNull();
    }

    private static void assertRetained(Supplier<?> ref) {
        System.gc();
        assertThat(ref.get()).isNotNull();
    }
}
