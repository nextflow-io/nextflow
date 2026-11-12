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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import groovyx.gpars.dataflow.operator.PoisonPill;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.Closures.action;
import static groovyx.gpars.TestThreads.awaitBlocked;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Pins the GPars {@link DataflowQueue} semantics that Nextflow relies on
 * (queue channels, process control channels).
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class DataflowQueueConformanceTest {

    static final long WAIT = 5;

    static List<Object> read(DataflowReadChannel<Object> ch, int n) throws InterruptedException {
        var result = new ArrayList<Object>(n);
        for( int i = 0; i < n; i++ )
            result.add(ch.getVal(WAIT, TimeUnit.SECONDS));
        return result;
    }

    static List<Object> range(int n) {
        return new ArrayList<>(IntStream.range(0, n).boxed().toList());
    }

    @Test
    void shouldDeliverItemsInFifoOrder() throws Exception {
        var q = new DataflowQueue<Object>();
        var n = 10_000;

        for( int i = 0; i < n; i++ )
            q.bind(i);
        assertThat(read(q, n)).isEqualTo(range(n));
    }

    @Test
    void shouldPreserveFifoOrderWithAConcurrentWriter() throws Exception {
        var q = new DataflowQueue<Object>();
        var n = 10_000;
        start(() -> { for( int i = 0; i < n; i++ ) q.bind(i); });

        assertThat(read(q, n)).isEqualTo(range(n));
    }

    @Test
    void shouldChainLeftShift() throws Exception {
        var q = new DataflowQueue<Object>();

        var result = q.leftShift(1).leftShift(2).leftShift(3);
        assertThat(result).isSameAs(q);
        assertThat(q.getVal()).isEqualTo(1);
        assertThat(q.getVal()).isEqualTo(2);
        assertThat(q.getVal()).isEqualTo(3);
    }

    @Test
    void shouldReportTheNumberOfBufferedItemsWithLength() throws Exception {
        var q = new DataflowQueue<Object>();
        assertThat(q.length()).isZero();
        assertThat(q.isBound()).isFalse();

        q.leftShift("a").leftShift("b").leftShift("c");
        assertThat(q.length()).isEqualTo(3);
        assertThat(q.isBound()).isTrue();

        q.getVal();
        assertThat(q.length()).isEqualTo(2);

        q.getVal();
        q.getVal();
        assertThat(q.length()).isZero();
        assertThat(q.isBound()).isFalse();
    }

    @Test
    void shouldReturnNullWhenGetValTimesOutOnAnEmptyQueue() throws Exception {
        var q = new DataflowQueue<Object>();

        var t0 = System.nanoTime();
        var result = q.getVal(100, TimeUnit.MILLISECONDS);
        var elapsed = TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - t0);
        assertThat(result).isNull();
        assertThat(elapsed).isGreaterThanOrEqualTo(90);

        q.bind("late");
        assertThat(q.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("late");
        assertThat(q.length()).isZero();
    }

    @Test
    void shouldBlockGetValUntilAnItemIsBound() throws Exception {
        var q = new DataflowQueue<Object>();
        var reader = new DataflowVariable<Object>();
        var thread = start(() -> reader.bind(q.getVal()));

        awaitBlocked(thread);
        assertThat(reader.isBound()).isFalse();

        q.bind("x");
        assertThat(reader.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("x");
    }

    @Test
    void shouldConsumeOneItemPerWhenBoundCallback() throws Exception {
        var q = new DataflowQueue<Object>();
        var received = new DataflowVariable<Object>();

        q.whenBound(action(x -> received.bind(x)));
        q.leftShift("a").leftShift("b");
        assertThat(received.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("a");
        assertThat(q.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("b");
        assertThat(q.length()).isZero();
    }

    @Test
    void shouldConsumeAnAlreadyBufferedItemWithWhenBound() throws Exception {
        var q = new DataflowQueue<Object>();
        q.leftShift("a").leftShift("b");
        var received = new DataflowVariable<Object>();

        q.whenBound(action(x -> received.bind(x)));
        assertThat(received.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("a");
        assertThat(q.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("b");
    }

    @Test
    void shouldHandEachItemToExactlyOnePendingWhenBoundCallback() throws Exception {
        var q = new DataflowQueue<Object>();
        var received = new CopyOnWriteArrayList<Object>();
        var done = new CountDownLatch(3);
        for( int i = 0; i < 3; i++ )
            q.whenBound(action(x -> { received.add(x); done.countDown(); }));

        q.leftShift(1).leftShift(2).leftShift(3).leftShift(4);
        assertThat(done.await(WAIT, TimeUnit.SECONDS)).isTrue();
        assertThat(received).containsExactlyInAnyOrder(1, 2, 3);
        assertThat(q.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo(4);
    }

    @Test
    void shouldPassAPoisonPillAsAnOrdinaryItem() throws Exception {
        var q = new DataflowQueue<Object>();

        q.leftShift(1).leftShift(PoisonPill.getInstance()).leftShift(2);
        assertThat(q.length()).isEqualTo(3);
        assertThat(q.getVal()).isEqualTo(1);
        assertThat(q.getVal()).isSameAs(PoisonPill.getInstance());
        assertThat(q.getVal()).isEqualTo(2);
    }

    @Test
    void shouldCopyEveryItemAndThePoisonPillWithInto() throws Exception {
        var q = new DataflowQueue<Object>();
        var target = new DataflowBroadcast<Object>();
        var reader = target.createReadChannel();

        q.into(target);
        q.leftShift(1).leftShift(2).leftShift(3).leftShift(PoisonPill.getInstance());
        assertThat(read(reader, 4)).containsExactly(1, 2, 3, PoisonPill.getInstance());
    }

}
