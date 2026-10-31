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

import static groovyx.gpars.TestThreads.poll;
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
            result.add(poll(ch, WAIT, TimeUnit.SECONDS));
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
    void shouldReportTheNumberOfBufferedItemsWithLength() throws Exception {
        var q = new DataflowQueue<Object>();
        assertThat(q.length()).isZero();
        assertThat(q.isBound()).isFalse();

        q.bind("a");
        q.bind("b");
        q.bind("c");
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
    void shouldUnsubscribeAnExpiredReadWithoutConsumingALaterValue() throws Exception {
        var q = new DataflowQueue<Object>();

        assertThat(poll(q, 50, TimeUnit.MILLISECONDS)).isNull();

        q.bind("late");
        assertThat(q.length()).isEqualTo(1);
        assertThat(poll(q, WAIT, TimeUnit.SECONDS)).isEqualTo("late");
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
        assertThat(poll(reader, WAIT, TimeUnit.SECONDS)).isEqualTo("x");
    }

    @Test
    void shouldConsumeOneItemPerWhenBoundCallback() throws Exception {
        var q = new DataflowQueue<Object>();
        var received = new DataflowVariable<Object>();

        q.whenBound(action(x -> received.bind(x)));
        q.bind("a");
        q.bind("b");
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isEqualTo("a");
        assertThat(poll(q, WAIT, TimeUnit.SECONDS)).isEqualTo("b");
        assertThat(q.length()).isZero();
    }

    @Test
    void shouldConsumeAnAlreadyBufferedItemWithWhenBound() throws Exception {
        var q = new DataflowQueue<Object>();
        q.bind("a");
        q.bind("b");
        var received = new DataflowVariable<Object>();

        q.whenBound(action(x -> received.bind(x)));
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isEqualTo("a");
        assertThat(poll(q, WAIT, TimeUnit.SECONDS)).isEqualTo("b");
    }

    @Test
    void shouldHandEachItemToExactlyOnePendingWhenBoundCallback() throws Exception {
        var q = new DataflowQueue<Object>();
        var received = new CopyOnWriteArrayList<Object>();
        var done = new CountDownLatch(3);
        for( int i = 0; i < 3; i++ )
            q.whenBound(action(x -> { received.add(x); done.countDown(); }));

        q.bind(1);
        q.bind(2);
        q.bind(3);
        q.bind(4);
        assertThat(done.await(WAIT, TimeUnit.SECONDS)).isTrue();
        assertThat(received).containsExactlyInAnyOrder(1, 2, 3);
        assertThat(poll(q, WAIT, TimeUnit.SECONDS)).isEqualTo(4);
    }

    @Test
    void shouldPassAPoisonPillAsAnOrdinaryItem() throws Exception {
        var q = new DataflowQueue<Object>();

        q.bind(1);
        q.bind(PoisonPill.getInstance());
        q.bind(2);
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
        q.bind(1);
        q.bind(2);
        q.bind(3);
        q.bind(PoisonPill.getInstance());
        assertThat(read(reader, 4)).containsExactly(1, 2, 3, PoisonPill.getInstance());
    }

}
