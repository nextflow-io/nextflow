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
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import groovyx.gpars.Closures;
import groovyx.gpars.agent.TestAgents;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.operator.DataflowEventListener;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import groovyx.gpars.dataflow.operator.PoisonPill;
import groovyx.gpars.dataflow.operator.TestOperators;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.TestThreads.poll;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

@Timeout(value = 10, unit = TimeUnit.SECONDS)
class SchedulerTest {

    static final PoisonPill STOP = PoisonPill.getInstance();

    ExecutorService loop;

    @BeforeEach
    void setUp() {
        loop = Scheduler.newPool(1);
    }

    @AfterEach
    void tearDown() {
        loop.shutdown();
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
        assertThat(onPool(loop, Scheduler::isLoopThread)).isTrue();

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

        assertThat(onPool(loop, () -> (Object) bound.getVal())).isEqualTo("x");
        assertThat(onPool(loop, () -> {
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
            assertThat(onPool(loop, () -> {
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
        assertThat(poll(queue, 5, TimeUnit.SECONDS)).isEqualTo("x");
        assertThat(poll(reader, 5, TimeUnit.SECONDS)).isEqualTo("y");
    }

    @Test
    void shouldRunOperatorsForksAndAgentsOnTheLoopThread() throws Exception {
        var input = new DataflowQueue<Object>();
        var middle = new DataflowQueue<Object>();
        var output = new DataflowQueue<Object>();
        var threads = new CopyOnWriteArrayList<Thread>();
        var agent = TestAgents.create(loop, 0);

        TestOperators.start(loop, List.of(input), List.of(middle), List.of(), 1, Closures.closure(1, (self, args) -> {
            threads.add(Thread.currentThread());
            ((DataflowProcessor) self).bindOutput((Integer) args[0] * 10);
            return null;
        }));
        // with maxForks, the bodies are queued on the loop and run one at a time
        TestOperators.start(loop, List.of(middle), List.of(output), List.of(), 4, Closures.closure(1, (self, args) -> {
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
            result.add(poll(output, 5, TimeUnit.SECONDS));
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
        var listener = new DataflowEventListener() {
            @Override
            public boolean onException(DataflowProcessor processor, Throwable e) {
                error.bind(e);
                return true;
            }
        };

        TestOperators.start(loop, List.of(input), List.of(), List.of(listener), 1, Closures.action(x -> unbound.getVal()));
        input.bind(1);
        assertThat(poll(error, 5, TimeUnit.SECONDS))
            .isInstanceOf(IllegalStateException.class)
            .hasMessageContaining(Scheduler.POOL_SIZE_ENV);
    }

    @Test
    void shouldNotLeakAnInterruptSetByAnAgentClosureToTheNextMessageOrPoolTask() throws Exception {
        var agent = TestAgents.create(loop, 0);
        var interrupted = new CopyOnWriteArrayList<Boolean>();
        // hold the loop thread, so that the agent handles both messages in one batch
        var hold = new CountDownLatch(1);
        loop.execute(() -> {
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
        assertThat(onPool(loop, () -> Thread.currentThread().isInterrupted())).isFalse();
    }

    @Test
    void shouldHandleTheAgentMessagesInPlaceWhenTheLoopThreadWaitsForThem() throws Exception {
        var agent = TestAgents.create(loop, new ArrayList<Object>());
        // the messages are queued behind the task that waits for them
        var result = onPool(loop, () -> {
            for( int i = 0; i < 3; i++ ) {
                final int x = i;
                agent.send(Closures.action(it -> ((List<Object>) it).add(x)));
            }
            agent.await();
            var awaited = List.copyOf((List<?>) agent.getVal());
            agent.send(Closures.action(it -> ((List<Object>) it).add(3)));
            return List.of(awaited, agent.getVal(), Thread.currentThread());
        });

        assertThat(result.get(0)).isEqualTo(List.of(0, 1, 2));
        assertThat(result.get(1)).isEqualTo(List.of(0, 1, 2, 3));
        assertThat(onPool(loop, Thread::currentThread)).isSameAs(result.get(2));
    }

    @Test
    void shouldThrowWhenTheLoopThreadIsInterruptedBeforeWaitingForAnAgent() throws Exception {
        var agent = TestAgents.create(loop, 0);
        var result = onPool(loop, () -> {
            Thread.currentThread().interrupt();
            try {
                agent.await();
                return null;
            }
            catch( InterruptedException e ) {
                return List.of(e, Thread.currentThread().isInterrupted());
            }
        });

        assertThat((List<?>) result).first().isInstanceOf(InterruptedException.class);
        assertThat((List<?>) result).last().isEqualTo(false);
        assertThat(agent.hasErrors()).isFalse();
        assertThat(agent.getVal()).isEqualTo(0);
    }

    @Test
    void shouldNameTheDataflowThreads() throws Exception {
        var pool = Scheduler.newDefaultPool();
        try {
            var thread = onPool(pool, Thread::currentThread);
            assertThat(thread.getName()).matches("Dataflow Thread \\d+");
            assertThat(thread.isDaemon()).isTrue();
            assertThat(onPool(pool, Scheduler::isLoopThread)).isFalse();
        }
        finally {
            pool.shutdown();
        }
        var thread = onPool(loop, Thread::currentThread);
        assertThat(thread.getName()).matches("Dataflow Thread \\d+");
        assertThat(thread.isDaemon()).isTrue();
        // a custom handler, not the thread group
        assertThat(thread.getUncaughtExceptionHandler()).isNotSameAs(thread.getThreadGroup());
        assertThat(Scheduler.newThreadName()).isNotEqualTo(Scheduler.newThreadName());
    }

    interface Call<T> { T call() throws Exception; }

    static <T> T onPool(Executor pool, Call<T> call) throws Exception {
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
        var value = (T) poll(result, 5, TimeUnit.SECONDS);
        return value;
    }
}
