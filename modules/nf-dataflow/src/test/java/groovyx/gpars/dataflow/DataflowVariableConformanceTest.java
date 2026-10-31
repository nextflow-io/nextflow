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

import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static groovyx.gpars.TestThreads.poll;
import static groovyx.gpars.Closures.action;
import static groovyx.gpars.TestThreads.awaitBlocked;
import static groovyx.gpars.TestThreads.start;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Pins the GPars {@link DataflowVariable} semantics that Nextflow relies on
 * (value channels and reply variables).
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
        assertThat(poll(v, WAIT, TimeUnit.SECONDS)).isEqualTo("x");
        assertThat(v.get()).isEqualTo("x");
        assertThat(v.isError()).isFalse();
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
    void shouldUnsubscribeAnExpiredReadWithoutConsumingALaterValue() throws Exception {
        var v = new DataflowVariable<Object>();

        assertThat(poll(v, 50, TimeUnit.MILLISECONDS)).isNull();
        assertThat(v.isBound()).isFalse();

        v.bind("late");
        assertThat(poll(v, WAIT, TimeUnit.SECONDS)).isEqualTo("late");
    }

    @Test
    void shouldBlockGetValUntilAnotherThreadBinds() throws Exception {
        var v = new DataflowVariable<Object>();
        var reader = new DataflowVariable<Object>();
        var thread = start(() -> reader.bind(v.getVal()));

        awaitBlocked(thread);
        assertThat(reader.isBound()).isFalse();

        v.bind("done");
        assertThat(poll(reader, WAIT, TimeUnit.SECONDS)).isEqualTo("done");
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
        assertThat(poll(error, WAIT, TimeUnit.SECONDS)).isInstanceOf(InterruptedException.class);
    }

    @Test
    void shouldFireWhenBoundRegisteredBeforeTheBind() throws Exception {
        var v = new DataflowVariable<Object>();
        var received = new DataflowVariable<Object>();
        v.whenBound(action(received::bind));

        assertThat(poll(received, 50, TimeUnit.MILLISECONDS)).isNull();

        v.bind("hello");
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isEqualTo("hello");
    }

    @Test
    void shouldFireWhenBoundRegisteredAfterTheBind() throws Exception {
        var v = new DataflowVariable<Object>();
        v.bind("hello");
        var received = new DataflowVariable<Object>();

        v.whenBound(action(received::bind));
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isEqualTo("hello");
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
        assertThat(poll(v, WAIT, TimeUnit.SECONDS)).isNull();

        assertThatThrownBy(v::get).isSameAs(error);
    }

    @Test
    void shouldPassTheBoundErrorToWhenBoundCallbacks() throws Exception {
        var v = new DataflowVariable<Object>();
        var error = new RuntimeException("boom");
        var received = new DataflowVariable<Object>();
        v.whenBound(action(received::bind));

        v.bindError(error);
        assertThat(poll(received, WAIT, TimeUnit.SECONDS)).isSameAs(error);
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

}
