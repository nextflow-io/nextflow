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
package groovyx.gpars.agent;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

import groovy.lang.Closure;
import groovyx.gpars.Closures;
import groovyx.gpars.dataflow.DataflowVariable;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.TestThreads.poll;
import static groovyx.gpars.Closures.action;
import static groovyx.gpars.TestThreads.awaitBlocked;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Pins the GPars {@link Agent} semantics that Nextflow relies on: the {@code TaskProcessor} state
 * agent ({@code ChannelEx.update}), and the cache, trace and Tower report writers.
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class AgentConformanceTest {

    @FunctionalInterface
    interface Update<T> { T apply(T state); }

    /**
     * The {@code ChannelEx.update} pattern: the agent is the closure delegate.
     */
    @SuppressWarnings("unchecked")
    static <T> Closure<Object> updateValue(Update<T> fn) {
        return Closures.closure(1, (delegate, args) -> {
            final var agent = (Agent<T>) delegate;
            agent.updateValue(fn.apply((T) args[0]));
            return null;
        });
    }

    @SuppressWarnings("unchecked")
    static List<Object> list(Object state) {
        return (List<Object>) state;
    }

    @Test
    void shouldRunSentClosuresInSendOrderOnTheAgentState() throws Exception {
        var agent = new Agent<List<Object>>(new ArrayList<>());
        var n = 1_000;

        for( int i = 0; i < n; i++ ) {
            final int x = i;
            agent.send(action(it -> list(it).add(x)));
        }
        assertThat(agent.getVal()).isEqualTo(IntStream.range(0, n).boxed().toList());
    }

    @Test
    void shouldRunSentClosuresSerially() throws Exception {
        var agent = new Agent<Integer>(0);
        var running = new AtomicInteger();
        var maxRunning = new AtomicInteger();
        var calls = new AtomicInteger();
        var senders = new ArrayList<Thread>();
        for( int i = 0; i < 4; i++ ) {
            senders.add(start(() -> {
                for( int j = 0; j < 100; j++ ) {
                    agent.send(action(it -> {
                        maxRunning.accumulateAndGet(running.incrementAndGet(), Math::max);
                        calls.incrementAndGet();
                        running.decrementAndGet();
                    }));
                }
            }));
        }

        for( var sender : senders )
            sender.join();
        agent.await();
        assertThat(calls.get()).isEqualTo(400);
        assertThat(maxRunning.get()).isEqualTo(1);
    }

    @Test
    void shouldWaitForAllPreviouslySentMessagesWithAwait() throws Exception {
        var agent = new Agent<List<Object>>(new ArrayList<>());
        var caller = Thread.currentThread();
        var release = new CountDownLatch(1);
        // the first message only completes once the caller is blocked in await()
        start(() -> { awaitBlocked(caller); release.countDown(); });

        agent.send(action(it -> { release.await(); list(it).add("slow"); }));
        agent.send(action(it -> list(it).add("fast")));
        agent.await();
        assertThat(agent.getVal()).containsExactly("slow", "fast");
    }

    @Test
    void shouldReplaceTheStateAndNotifyListenersWithUpdateValue() throws Exception {
        var agent = new Agent<Integer>(1);
        var changes = new CopyOnWriteArrayList<List<Object>>();
        agent.addListener(action((old, value) -> changes.add(List.of(old, value))));

        agent.send(updateValue((Integer it) -> it + 1));
        agent.send(updateValue((Integer it) -> it * 10));
        assertThat(agent.getVal()).isEqualTo(20);
        assertThat(changes).containsExactly(List.of(1, 2), List.of(2, 20));
    }

    @Test
    void shouldPassTheAgentToAThreeArgumentListener() throws Exception {
        var agent = new Agent<String>("a");
        var seen = new DataflowVariable<List<Object>>();
        agent.addListener(action((self, old, value) -> seen.bind(List.of(self, old, value))));

        agent.send(updateValue((String it) -> "b"));
        var result = poll(seen, 5, TimeUnit.SECONDS);
        assertThat(result.get(0)).isSameAs(agent);
        assertThat(result.get(1)).isEqualTo("a");
        assertThat(result.get(2)).isEqualTo("b");
    }

    @Test
    void shouldReplaceTheStateWhenSentAPlainValue() throws Exception {
        var agent = new Agent<String>("a");
        var changes = new CopyOnWriteArrayList<List<Object>>();
        agent.addListener(action((old, value) -> changes.add(List.of(old, value))));

        agent.send("b");
        agent.send("c");
        assertThat(agent.getVal()).isEqualTo("c");
        assertThat(changes).containsExactly(List.of("a", "b"), List.of("b", "c"));
    }

    @Test
    @Tag("gpars-incidental")
    void shouldNotNotifyListenersWhenTheClosureMutatesTheStateInPlace() throws Exception {
        var agent = new Agent<List<Object>>(new ArrayList<>());
        var calls = new AtomicInteger();
        agent.addListener(action((old, value) -> calls.incrementAndGet()));

        agent.send(action(it -> list(it).add(1)));
        agent.await();
        assertThat(agent.getVal()).containsExactly(1);
        assertThat(calls.get()).isZero();
    }

    @Test
    void shouldKeepProcessingMessagesAfterAFailingClosure() throws Exception {
        var agent = new Agent<List<Object>>(new ArrayList<>());

        agent.send(action(it -> list(it).add(1)));
        agent.send(action(it -> { throw new IllegalStateException("agent boom"); }));
        agent.send(action(it -> list(it).add(2)));
        assertThat(agent.getVal()).containsExactly(1, 2);
        assertThat(agent.hasErrors()).isTrue();
        assertThat(agent.getErrors()).extracting(Throwable::getMessage).containsExactly("agent boom");
        assertThat(agent.hasErrors()).isFalse();
    }

    @Test
    void shouldWorkWithANullInitialState() throws Exception {
        // the CacheDB writer pattern
        var agent = new Agent<Object>();
        var done = new DataflowVariable<Boolean>();

        agent.send(action(it -> done.bind(it == null)));
        agent.await();
        assertThat(poll(done, 5, TimeUnit.SECONDS)).isTrue();
    }

    @Test
    @Tag("gpars-incidental")
    void shouldFailWhenAwaitingFromInsideTheAgent() throws Exception {
        var agent = new Agent<Integer>(0);
        var error = new DataflowVariable<Object>();

        agent.send(action(it -> {
            try {
                agent.await();
                error.bind(null);
            }
            catch( Throwable e ) {
                error.bind(e);
            }
        }));
        assertThat(poll(error, 5, TimeUnit.SECONDS)).isInstanceOf(IllegalStateException.class);
    }

}
