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
package groovyx.gpars.scheduler;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import groovyx.gpars.Closures;
import groovyx.gpars.agent.Agent;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.operator.DataflowEventAdapter;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import groovyx.gpars.dataflow.operator.PoisonPill;
import groovyx.gpars.group.DefaultPGroup;
import groovyx.gpars.group.PGroup;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.junit.jupiter.api.Assumptions.assumeFalse;

@Timeout(value = 10, unit = TimeUnit.SECONDS)
class SchedulerTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    PGroup loop;

    @BeforeEach
    void setUp() {
        loop = new DefaultPGroup(Scheduler.newPool(1));
    }

    @AfterEach
    void tearDown() {
        loop.shutdown();
    }

    @Test
    void shouldParseTheEngine() {
        assertThat(Scheduler.parseEngine(null)).isFalse();
        assertThat(Scheduler.parseEngine("")).isFalse();
        assertThat(Scheduler.parseEngine("mailbox")).isFalse();
        assertThat(Scheduler.parseEngine("legacy")).isTrue();
        assertThatThrownBy(() -> Scheduler.parseEngine("actors")).isInstanceOf(IllegalArgumentException.class);
    }

    @Test
    void shouldParseThePoolSize() {
        assertThat(Scheduler.parsePoolSize(null)).isZero();
        assertThat(Scheduler.parsePoolSize("")).isZero();
        assertThat(Scheduler.parsePoolSize("1")).isEqualTo(1);
        assertThat(Scheduler.parsePoolSize(" 8 ")).isEqualTo(8);
        assertThatThrownBy(() -> Scheduler.parsePoolSize("0")).isInstanceOf(IllegalArgumentException.class);
        assertThatThrownBy(() -> Scheduler.parsePoolSize("many")).isInstanceOf(IllegalArgumentException.class);
    }

    @Test
    void shouldFlagOnlyTheThreadOfASingleThreadPoolAsTheLoopThread() throws Exception {
        assertThat(Scheduler.isLoopThread()).isFalse();
        assertThat(onPool(loop.getThreadPool(), Scheduler::isLoopThread)).isTrue();

        var pool = Scheduler.newPool(2);
        try {
            assertThat(onPool(pool, Scheduler::isLoopThread)).isFalse();
        }
        finally {
            pool.shutdown();
        }
    }

    @Test
    void shouldFailInsteadOfBlockingTheLoopThreadOnAnUnboundValue() throws Exception {
        var unbound = new DataflowVariable<Object>();
        var bound = new DataflowVariable<Object>();
        bound.bind("x");

        assertThat(onPool(loop.getThreadPool(), () -> (Object) bound.getVal())).isEqualTo("x");
        assertThat(onPool(loop.getThreadPool(), () -> {
            try {
                return unbound.getVal();
            }
            catch( IllegalStateException e ) {
                return e;
            }
        })).isInstanceOf(IllegalStateException.class);
    }

    @Test
    void shouldFailInsteadOfBlockingTheLoopThreadOnAnEmptyQueueOrStreamWithoutTakingTheNextValue() throws Exception {
        var queue = new DataflowQueue<Object>();
        var broadcast = new DataflowBroadcast<Object>();
        var reader = broadcast.createReadChannel();

        for( var channel : List.<DataflowReadChannel<Object>>of(queue, reader) ) {
            assertThat(onPool(loop.getThreadPool(), () -> {
                try {
                    return channel.getVal();
                }
                catch( IllegalStateException e ) {
                    return e;
                }
            })).isInstanceOf(IllegalStateException.class);
        }
        queue.bind("x");
        broadcast.bind("y");
        assertThat(queue.getVal(5, TimeUnit.SECONDS)).isEqualTo("x");
        assertThat(reader.getVal(5, TimeUnit.SECONDS)).isEqualTo("y");
    }

    @Test
    void shouldRunOperatorsForksAndAgentsOnTheLoopThread() throws Exception {
        // the legacy engine blocks the actor thread while all forks are busy, so it deadlocks on a single thread
        assumeFalse(Scheduler.isLegacyEngine());
        var input = new DataflowQueue<Object>();
        var middle = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var threads = new CopyOnWriteArrayList<Thread>();
        var agent = loop.agent(0);

        loop.operator(Map.of("inputs", List.of(input), "outputs", List.of(middle)), Closures.closure(1, (self, args) -> {
            threads.add(Thread.currentThread());
            ((DataflowProcessor) self).bindOutput((Integer) args[0] * 10);
            return null;
        }));
        // with maxForks, the bodies are queued on the loop and run one at a time
        loop.operator(Map.of("inputs", List.of(middle), "outputs", List.of(output), "maxForks", 4), Closures.closure(1, (self, args) -> {
            threads.add(Thread.currentThread());
            agent.send(Closures.action(it -> threads.add(Thread.currentThread())));
            ((DataflowProcessor) self).bindOutput((Integer) args[0] + 1);
            return null;
        }));
        for( int i = 0; i < 100; i++ )
            input.bind(i);
        input.bind(STOP);

        var result = new ArrayList<Object>();
        for( int i = 0; i < 101; i++ )
            result.add(output.getVal(5, TimeUnit.SECONDS));
        assertThat(result.get(100)).isSameAs(STOP);
        assertThat(result.subList(0, 100)).containsExactlyInAnyOrderElementsOf(IntStream.range(0, 100).map(i -> i * 10 + 1).boxed().toList());
        agent.await();
        assertThat(threads).hasSize(300).allMatch(t -> t == threads.get(0));
    }

    @Test
    void shouldReportABlockingReadInAnOperatorOnTheLoopThreadAsAnError() throws Exception {
        var input = new DataflowQueue<Object>();
        var unbound = new DataflowVariable<Object>();
        var error = new DataflowVariable<Throwable>();
        var listener = new DataflowEventAdapter() {
            @Override
            public boolean onException(DataflowProcessor processor, Throwable e) {
                error.bind(e);
                return true;
            }
        };

        loop.operator(Map.of("inputs", List.of(input), "listeners", List.of(listener)), Closures.action(x -> unbound.getVal()));
        input.bind(1);
        assertThat(error.getVal(5, TimeUnit.SECONDS))
            .isInstanceOf(IllegalStateException.class)
            .hasMessageContaining(Scheduler.POOL_SIZE_ENV);
    }

    @Test
    void shouldNotLeakAnInterruptSetByAnAgentClosureToTheNextMessageOrPoolTask() throws Exception {
        var agent = loop.agent(0);
        var interrupted = new CopyOnWriteArrayList<Boolean>();
        // hold the loop thread, so that the agent handles both messages in one batch
        var hold = new CountDownLatch(1);
        loop.getThreadPool().execute(() -> {
            try {
                hold.await();
            }
            catch( InterruptedException e ) {
                throw new RuntimeException(e);
            }
        });
        agent.send(Closures.action(it -> Thread.currentThread().interrupt()));
        agent.send(Closures.action(it -> interrupted.add(Thread.currentThread().isInterrupted())));
        hold.countDown();
        agent.await();

        assertThat(interrupted).containsExactly(false);
        assertThat(agent.getErrors()).singleElement().isInstanceOf(InterruptedException.class);
        assertThat(onPool(loop.getThreadPool(), () -> Thread.currentThread().isInterrupted())).isFalse();
    }

    interface Call<T> { T call() throws Exception; }

    static <T> T onPool(Pool pool, Call<T> call) throws Exception {
        var result = new DataflowVariable<Object>();
        pool.execute(() -> {
            try {
                result.bind(call.call());
            }
            catch( Throwable e ) {
                result.bind(e);
            }
        });
        @SuppressWarnings("unchecked")
        var value = (T) result.getVal(5, TimeUnit.SECONDS);
        return value;
    }
}
