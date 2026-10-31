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
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.Closures.action;
import static groovyx.gpars.Closures.function;
import static groovyx.gpars.TestThreads.awaitBlocked;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Pins the GPars {@link DataflowVariable} semantics that Nextflow relies on
 * (value channels, reply variables, {@code LazyDataflowVariable} caches).
 */
@Timeout(value = 10, unit = TimeUnit.SECONDS)
class DataflowVariableConformanceTest {

    static final long WAIT = 5;

    @Test
    void shouldBindValueOnceAndReadItManyTimes() throws Throwable {
        var v = new DataflowVariable<Object>();
        assertThat(v.isBound()).isFalse();

        v.bind("x");
        assertThat(v.isBound()).isTrue();
        assertThat(v.getVal()).isEqualTo("x");
        assertThat(v.getVal()).isEqualTo("x");
        assertThat(v.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("x");
        assertThat(v.get()).isEqualTo("x");
        assertThat(v.isError()).isFalse();
    }

    @Test
    void shouldBindWithLeftShift() throws Exception {
        var v = new DataflowVariable<Object>();

        var result = v.leftShift(42);
        assertThat(result).isSameAs(v);
        assertThat(v.getVal()).isEqualTo(42);
    }

    @Test
    void shouldBindNullValue() throws Exception {
        var v = new DataflowVariable<Object>();

        v.bind(null);
        assertThat(v.isBound()).isTrue();
        assertThat(v.getVal()).isNull();
    }

    @Test
    void shouldFailWhenRebindingDifferentValue() throws Exception {
        var v = new DataflowVariable<Object>();
        v.bind(1);

        assertThatThrownBy(() -> v.bind(2)).isInstanceOf(IllegalStateException.class);
        assertThat(v.getVal()).isEqualTo(1);
    }

    @Test
    @Tag("gpars-incidental")
    void shouldTolerateRebindingEqualValue() throws Exception {
        var v = new DataflowVariable<Object>();
        v.bind("a");

        v.bind("a");
        v.bind(new String("a"));
        assertThat(v.getVal()).isEqualTo("a");
    }

    @Test
    void shouldReturnNullWhenGetValTimesOutOnUnboundVariable() throws Exception {
        var v = new DataflowVariable<Object>();

        var t0 = System.nanoTime();
        var result = v.getVal(100, TimeUnit.MILLISECONDS);
        var elapsed = TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - t0);
        assertThat(result).isNull();
        assertThat(elapsed).isGreaterThanOrEqualTo(90);
        assertThat(v.isBound()).isFalse();
    }

    @Test
    void shouldBlockGetValUntilAnotherThreadBinds() throws Exception {
        var v = new DataflowVariable<Object>();
        var reader = new DataflowVariable<Object>();
        var thread = start(() -> reader.bind(v.getVal()));

        awaitBlocked(thread);
        assertThat(reader.isBound()).isFalse();

        v.bind("done");
        assertThat(reader.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("done");
    }

    @Test
    void shouldThrowInterruptedExceptionWhenBlockedGetValIsInterrupted() throws Exception {
        var v = new DataflowVariable<Object>();
        var error = new DataflowVariable<Object>();
        var thread = start(() -> {
            try {
                v.getVal();
                error.bind("no error");
            }
            catch( InterruptedException e ) {
                error.bind(e);
            }
        });

        awaitBlocked(thread);
        thread.interrupt();
        assertThat(error.getVal(WAIT, TimeUnit.SECONDS)).isInstanceOf(InterruptedException.class);
    }

    @Test
    void shouldFireWhenBoundRegisteredBeforeTheBind() throws Exception {
        var v = new DataflowVariable<Object>();
        var received = new DataflowVariable<Object>();
        v.whenBound(action(received::bind));

        assertThat(received.getVal(50, TimeUnit.MILLISECONDS)).isNull();

        v.bind("hello");
        assertThat(received.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("hello");
    }

    @Test
    void shouldFireWhenBoundRegisteredAfterTheBind() throws Exception {
        var v = new DataflowVariable<Object>();
        v.bind("hello");
        var received = new DataflowVariable<Object>();

        v.whenBound(action(received::bind));
        assertThat(received.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("hello");
    }

    @Test
    void shouldFireEveryWhenBoundCallbackExactlyOnce() throws Exception {
        var v = new DataflowVariable<Object>();
        var calls = new CopyOnWriteArrayList<String>();
        var done = new CountDownLatch(4);
        for( int i = 0; i < 2; i++ ) {
            final var name = "before-" + i;
            v.whenBound(action(x -> { calls.add(name + ":" + x); done.countDown(); }));
        }

        v.bind("x");
        for( int i = 0; i < 2; i++ ) {
            final var name = "after-" + i;
            v.whenBound(action(x -> { calls.add(name + ":" + x); done.countDown(); }));
        }
        assertThat(done.await(WAIT, TimeUnit.SECONDS)).isTrue();
        // a short wait to catch a duplicate callback: it can only miss a bug, never cause a false failure
        Thread.sleep(50);
        assertThat(calls).containsExactlyInAnyOrder("after-0:x", "after-1:x", "before-0:x", "before-1:x");
    }

    @Test
    void shouldBindAnError() throws Exception {
        var v = new DataflowVariable<Object>();
        var error = new RuntimeException("boom");

        v.bindError(error);
        assertThat(v.isBound()).isTrue();
        assertThat(v.isError()).isTrue();
        assertThat(v.getError()).isSameAs(error);
        // getVal does not throw: it returns null and the error is only visible through get()/isError()
        assertThat(v.getVal()).isNull();
        assertThat(v.getVal(WAIT, TimeUnit.SECONDS)).isNull();

        assertThatThrownBy(v::get).isSameAs(error);
    }

    @Test
    void shouldPassTheBoundErrorToWhenBoundCallbacks() throws Exception {
        var v = new DataflowVariable<Object>();
        var error = new RuntimeException("boom");
        var received = new DataflowVariable<Object>();
        v.whenBound(action(received::bind));

        v.bindError(error);
        assertThat(received.getVal(WAIT, TimeUnit.SECONDS)).isSameAs(error);
    }

    @Test
    void shouldFailWhenBindingAfterAnErrorOrAnErrorAfterAValue() throws Exception {
        var failed = new DataflowVariable<Object>();
        failed.bindError(new RuntimeException("boom"));
        var bound = new DataflowVariable<Object>();
        bound.bind(1);

        assertThatThrownBy(() -> failed.bind(1)).isInstanceOf(IllegalStateException.class);

        assertThatThrownBy(() -> bound.bindError(new RuntimeException("boom"))).isInstanceOf(IllegalStateException.class);
        assertThat(bound.getVal()).isEqualTo(1);
        assertThat(bound.isError()).isFalse();
    }

    @Test
    @Tag("gpars-incidental")
    void shouldTreatABoundThrowableValueAsAnError() throws Exception {
        var v = new DataflowVariable<Object>();
        var error = new IllegalArgumentException("x");

        v.bind(error);
        assertThat(v.isError()).isTrue();
        assertThat(v.getError()).isSameAs(error);
        assertThat(v.getVal()).isSameAs(error);
    }

    @Test
    void shouldChainThenOnTheBoundValue() throws Exception {
        // the ToListOp pattern: the result of then() is used as a DataflowVariable
        var before = new DataflowVariable<Object>();
        var early = before.then(function(x -> "early:" + x));
        before.bind(1);

        var after = new DataflowVariable<Object>();
        after.bind(2);
        var late = after.then(function(x -> "late:" + x));

        assertThat(early).isInstanceOf(DataflowVariable.class);
        assertThat(late).isInstanceOf(DataflowVariable.class);
        assertThat(((DataflowVariable<?>) early).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("early:1");
        assertThat(((DataflowVariable<?>) late).getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("late:2");
    }

    @Test
    void shouldPropagateErrorsThroughThen() throws Exception {
        var calls = new AtomicInteger();
        var failed = new DataflowVariable<Object>();
        var error = new RuntimeException("boom");
        var skipped = (DataflowVariable<?>) failed.then(function(x -> calls.incrementAndGet()));
        failed.bindError(error);

        var bound = new DataflowVariable<Object>();
        var thrown = (DataflowVariable<?>) bound.then(function(x -> { throw new IllegalStateException("then boom"); }));
        bound.bind(1);

        skipped.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(skipped.isError()).isTrue();
        assertThat(skipped.getError()).isSameAs(error);
        assertThat(calls.get()).isZero();

        thrown.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(thrown.isError()).isTrue();
        assertThat(thrown.getError()).isInstanceOf(IllegalStateException.class).hasMessage("then boom");
    }

    @Test
    void shouldComputeALazyVariableOnceOnFirstRead() throws Exception {
        var count = new AtomicInteger();
        var v = new LazyDataflowVariable<Object>(function(() -> { count.incrementAndGet(); return "value"; }));

        assertThat(count.get()).isZero();

        var a = v.getVal(WAIT, TimeUnit.SECONDS);
        var b = v.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(a).isEqualTo("value");
        assertThat(b).isEqualTo("value");
        assertThat(count.get()).isEqualTo(1);
        assertThat(v.isBound()).isTrue();
    }

    @Test
    void shouldComputeALazyVariableOnceUnderConcurrentReads() throws Exception {
        // the conda/spack/container cache pattern: concurrent tasks wait for the same environment
        var count = new AtomicInteger();
        var release = new CountDownLatch(1);
        var v = new LazyDataflowVariable<Object>(function(() -> {
            count.incrementAndGet();
            release.await(WAIT, TimeUnit.SECONDS);
            return "value";
        }));
        var results = new CopyOnWriteArrayList<Object>();
        var readers = new ArrayList<Thread>();
        for( int i = 0; i < 8; i++ )
            readers.add(start(() -> results.add(v.getVal())));

        for( var reader : readers )
            awaitBlocked(reader);
        release.countDown();
        for( var reader : readers )
            reader.join(WAIT * 1000);

        assertThat(results).hasSize(8).containsOnly("value");
        assertThat(count.get()).isEqualTo(1);
    }

    @Test
    void shouldComputeALazyVariableOnceWhenTheFirstReadsRace() throws Exception {
        for( int round = 0; round < 200; round++ ) {
            var count = new AtomicInteger();
            var v = new LazyDataflowVariable<Object>(function(() -> { count.incrementAndGet(); return "value"; }));
            var go = new CountDownLatch(1);
            var readers = new ArrayList<Thread>();
            for( int i = 0; i < 6; i++ ) {
                final int kind = i % 3;
                readers.add(start(() -> {
                    go.await();
                    if( kind == 0 )
                        v.getVal();
                    else if( kind == 1 )
                        v.touch();
                    else
                        v.isError();
                }));
            }
            go.countDown();
            for( var reader : readers )
                reader.join(WAIT * 1000);

            assertThat(v.getVal(WAIT, TimeUnit.SECONDS)).isEqualTo("value");
            assertThat(count.get()).isEqualTo(1);
        }
    }

    @Test
    void shouldReportTheErrorOfAFailingLazyVariable() throws Exception {
        var v = new LazyDataflowVariable<Object>(function(() -> { throw new IllegalStateException("lazy boom"); }));

        v.getVal(WAIT, TimeUnit.SECONDS);
        assertThat(v.isError()).isTrue();
        assertThat(v.getError()).isInstanceOf(IllegalStateException.class).hasMessage("lazy boom");
    }

}
