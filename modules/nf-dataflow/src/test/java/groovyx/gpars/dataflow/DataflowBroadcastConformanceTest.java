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
 * Pins the GPars {@link DataflowBroadcast} semantics that Nextflow relies on
 * (queue channels with several consumers, topics).
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class DataflowBroadcastConformanceTest {

    static final long WAIT = 5;

    static final PoisonPill STOP = PoisonPill.getInstance();

    static List<Object> read(DataflowReadChannel<Object> ch, int n) throws InterruptedException {
        var result = new ArrayList<Object>(n);
        for( int i = 0; i < n; i++ )
            result.add(poll(ch, WAIT, TimeUnit.SECONDS));
        return result;
    }

    @Test
    void shouldCreateIndependentReadChannels() {
        var ch = new DataflowBroadcast<Object>();

        var r1 = ch.createReadChannel();
        var r2 = ch.createReadChannel();
        assertThat(r1).isInstanceOf(DataflowReadChannel.class);
        assertThat(r2).isInstanceOf(DataflowReadChannel.class);
        assertThat(r1).isNotSameAs(r2);
    }

    @Test
    void shouldDeliverEveryItemInOrderToEveryReader() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var readers = List.of(ch.createReadChannel(), ch.createReadChannel(), ch.createReadChannel());
        var n = 5_000;

        start(() -> { for( int i = 0; i < n; i++ ) ch.bind(i); });
        var expected = IntStream.range(0, n).boxed().toList();
        for( var reader : readers )
            assertThat(read(reader, n)).isEqualTo(expected);
    }

    @Test
    void shouldNotDeliverItemsBoundBeforeTheReaderWasCreated() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var early = ch.createReadChannel();
        ch.bind(1);
        ch.bind(2);
        var late = ch.createReadChannel();
        ch.bind(3);
        assertThat(poll(late, WAIT, TimeUnit.SECONDS)).isEqualTo(3);
        assertThat(poll(late, 50, TimeUnit.MILLISECONDS)).isNull();
        assertThat(read(early, 3)).containsExactly(1, 2, 3);
    }

    @Test
    void shouldNotLetOneReaderConsumeItemsOfAnother() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var r1 = ch.createReadChannel();
        var r2 = ch.createReadChannel();
        ch.bind("a");
        ch.bind("b");
        assertThat(read(r1, 2)).containsExactly("a", "b");
        assertThat(read(r2, 2)).containsExactly("a", "b");
    }

    @Test
    void shouldDeliverThePoisonPillToEveryReader() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var r1 = ch.createReadChannel();
        var r2 = ch.createReadChannel();

        ch.bind(1);
        ch.bind(STOP);
        assertThat(poll(r1, WAIT, TimeUnit.SECONDS)).isEqualTo(1);
        assertThat(poll(r1, WAIT, TimeUnit.SECONDS)).isSameAs(STOP);
        assertThat(poll(r2, WAIT, TimeUnit.SECONDS)).isEqualTo(1);
        assertThat(poll(r2, WAIT, TimeUnit.SECONDS)).isSameAs(STOP);
    }

    @Test
    void shouldBlockAReaderUntilAnItemIsBoundAndTimeOutWithNull() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var r = ch.createReadChannel();

        assertThat(poll(r, 100, TimeUnit.MILLISECONDS)).isNull();

        var received = new DataflowVariable<Object>();
        var reader = start(() -> received.bind(r.getVal()));
        awaitBlocked(reader);
        assertThat(received.isBound()).isFalse();

        ch.bind("x");
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isEqualTo("x");
    }

    @Test
    void shouldNotLoseAValueAfterAnExpiredRead() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var r = ch.createReadChannel();

        assertThat(poll(r, 50, TimeUnit.MILLISECONDS)).isNull();

        ch.bind("a");
        ch.bind("b");
        assertThat(r.length()).isEqualTo(2);
        assertThat(read(r, 2)).containsExactly("a", "b");
    }

    @Test
    void shouldConsumeOneItemPerWhenBoundCallbackOnAReadChannel() throws Exception {
        var ch = new DataflowBroadcast<Object>();
        var r = ch.createReadChannel();
        var received = new DataflowVariable<Object>();

        r.whenBound(action(x -> received.bind(x)));
        ch.bind("a");
        ch.bind("b");
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isEqualTo("a");
        assertThat(poll(r, WAIT, TimeUnit.SECONDS)).isEqualTo("b");
    }

}
