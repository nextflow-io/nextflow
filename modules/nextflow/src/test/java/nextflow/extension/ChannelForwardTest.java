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

package nextflow.extension;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import groovyx.gpars.dataflow.operator.PoisonPill;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * {@link CH#forward}, which bridges queues to broadcasts and feeds the recursion channels,
 * and the copy closures of the copying operators
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class ChannelForwardTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    static List<Object> read(DataflowReadChannel<?> channel, int n) throws InterruptedException {
        final List<Object> result = new ArrayList<>();
        for( int i = 0; i < n; i++ )
            result.add(channel.getVal());
        return result;
    }

    @Test
    void shouldForwardEveryItemAndThePoisonPillOfAQueue() throws Exception {
        var queue = new DataflowQueue<Object>();
        var target = new DataflowBroadcast<Object>();
        var reader = target.createReadChannel();

        CH.forward(queue, target);
        for( var item : List.of(1, 2, 3, STOP) )
            queue.bind(item);
        assertThat(read(reader, 4)).containsExactly(1, 2, 3, STOP);
    }

    @Test
    void shouldForwardEveryItemOfABroadcastReaderWithoutConsumingTheOtherReaders() throws Exception {
        var broadcast = new DataflowBroadcast<Object>();
        var source = broadcast.createReadChannel();
        var other = broadcast.createReadChannel();
        var target = new DataflowQueue<Object>();

        CH.forward(source, target);
        for( var item : List.of(1, 2, STOP) )
            broadcast.bind(item);
        assertThat(read(target, 3)).containsExactly(1, 2, STOP);
        assertThat(read(other, 3)).containsExactly(1, 2, STOP);
    }

    @Test
    void shouldForwardAValueOnceBound() throws Exception {
        var value = new DataflowVariable<Object>();
        var target = new DataflowQueue<Object>();

        CH.forward(value, target);
        value.bind("x");
        assertThat(target.getVal()).isEqualTo("x");
    }

    @Test
    void shouldForwardABoundValueInTheCallerThread() {
        var value = new DataflowVariable<Object>();
        var target = new DataflowQueue<Object>();
        value.bind("x");

        CH.forward(value, target);
        assertThat(target.length()).isEqualTo(1);
    }

    @Test
    void shouldForwardTheStopOfAValueAndItsEndOfChannel() throws Exception {
        // a value bound to STOP, as by CH.close0 on an unbound value
        for( boolean bound : List.of(true, false) ) {
            var value = new DataflowVariable<Object>();
            var target = new DataflowVariable<Object>();
            if( bound )
                value.bind(STOP);
            CH.forward(value, target);
            if( !bound )
                CH.close0(value);
            assertThat(target.getVal()).isSameAs(STOP);
        }
    }

    @Test
    void shouldBindTheResultOfTheWrappedClosureToEveryOutput() throws Exception {
        var input = new DataflowQueue<Object>();
        var o1 = new DataflowQueue<Object>();
        var o2 = new DataflowQueue<Object>();
        var triple = new Closure<Object>(null) {
            @Override
            public Object call(Object... args) {
                return (Integer) args[0] * 3;
            }

            @Override
            public int getMaximumNumberOfParameters() {
                return 1;
            }
        };

        var copy = new DataflowQueue<Object>();

        new DataflowProcessor(List.of(input), List.of(o1, o2), null, 1, new ChainWithClosure<>(triple)).start();
        new DataflowProcessor(List.of(o2), List.of(copy), null, 1, new ChainWithClosure<>(new CopyChannelsClosure<>())).start();
        for( var item : List.of(1, 2, STOP) )
            input.bind(item);
        assertThat(read(o1, 3)).containsExactly(3, 6, STOP);
        assertThat(read(copy, 3)).containsExactly(3, 6, STOP);
    }
}
