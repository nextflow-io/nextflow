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
import java.util.concurrent.Executor;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import groovyx.gpars.Closures;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Agent behaviors beyond the GPars conformance suite.
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class AgentTest {

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

    @SuppressWarnings("unchecked")
    static List<Object> list(Object state) {
        return (List<Object>) state;
    }

    @ParameterizedTest
    @MethodSource("poolFailures")
    void shouldReportAFailedDrainAndHandleTheMessageInsteadOfFailingTheSender(Throwable failure) throws Exception {
        Executor pool = task -> { throw raise(failure); };
        var agent = new Agent<List<Object>>(pool, new ArrayList<>());

        agent.send(Closures.action(it -> list(it).add(1)));
        agent.send(Closures.action(it -> list(it).add(2)));

        assertThat(agent.getVal()).containsExactly(1, 2);
        assertThat(agent.getErrors()).isNotEmpty().allSatisfy(e -> assertThat(e).isSameAs(failure));
    }

    @ParameterizedTest
    @MethodSource("poolFailures")
    void shouldReportAFailedRescheduleAndHandleTheRemainingMessages(Throwable failure) throws Exception {
        var tasks = new ArrayList<Runnable>();
        var reject = new AtomicBoolean();
        Executor pool = task -> {
            if( reject.get() )
                throw raise(failure);
            tasks.add(task);
        };
        var agent = new Agent<List<Object>>(pool, new ArrayList<>());
        final int n = 1000;

        // more messages than a batch, so that the drain reschedules
        for( int i = 0; i < n; i++ ) {
            final int x = i;
            agent.send(Closures.action(it -> list(it).add(x)));
        }
        reject.set(true);
        tasks.remove(0).run();

        assertThat(tasks).isEmpty();
        assertThat(agent.getErrors()).isNotEmpty().allSatisfy(e -> assertThat(e).isSameAs(failure));
        assertThat(agent.getVal()).isEqualTo(IntStream.range(0, n).boxed().toList());
    }

    @Test
    void shouldKeepTheInterruptOfASenderThatRunsTheDrain() throws Exception {
        // like a saturated caller-runs pool
        Executor callerRuns = Runnable::run;
        var agent = new Agent<List<Object>>(callerRuns, new ArrayList<>());
        try {
            Thread.currentThread().interrupt();
            agent.send(Closures.action(it -> list(it).add(Thread.currentThread().isInterrupted())));
            agent.send(Closures.action(it -> list(it).add(2)));
            assertThat(Thread.currentThread().isInterrupted()).isTrue();

            // an interrupt raised by a closure is still an error, and the sender's interrupt is kept
            agent.send(Closures.action(it -> Thread.currentThread().interrupt()));
            agent.send(Closures.action(it -> list(it).add(Thread.currentThread().isInterrupted())));
            assertThat(Thread.currentThread().isInterrupted()).isTrue();
        }
        finally {
            Thread.interrupted();
        }

        assertThat(agent.getVal()).containsExactly(false, 2, false);
        assertThat(agent.getErrors()).singleElement().isInstanceOf(InterruptedException.class);
    }
}
