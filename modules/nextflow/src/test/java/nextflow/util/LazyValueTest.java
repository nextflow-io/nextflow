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

package nextflow.util;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

@Timeout(value = 30, unit = TimeUnit.SECONDS)
class LazyValueTest {

    @Test
    void shouldComputeOnceInTheThreadOfTheFirstReader() throws Exception {
        var count = new AtomicInteger();
        var computing = new CountDownLatch(1);
        var release = new CountDownLatch(1);
        var value = new LazyValue<String>(() -> {
            count.incrementAndGet();
            computing.countDown();
            release.await();
            return Thread.currentThread().getName();
        });

        ExecutorService readers = Executors.newFixedThreadPool(8);
        try {
            List<Future<String>> results = new ArrayList<>();
            results.add(readers.submit(value::get));
            computing.await();
            for( int i = 1; i < 8; i++ )
                results.add(readers.submit(value::get));
            release.countDown();

            var first = results.get(0).get(10, TimeUnit.SECONDS);
            for( var result : results )
                assertThat(result.get(10, TimeUnit.SECONDS)).isEqualTo(first);
            assertThat(first).startsWith("pool-");
            assertThat(count).hasValue(1);
        }
        finally {
            readers.shutdownNow();
        }
    }

    @Test
    void shouldGiveTheErrorToEveryReader() {
        var count = new AtomicInteger();
        var value = new LazyValue<String>(() -> {
            count.incrementAndGet();
            throw new IllegalArgumentException("lazy boom");
        });

        for( int i = 0; i < 3; i++ ) {
            assertThatThrownBy(value::get)
                .isInstanceOf(IllegalStateException.class)
                .cause().isInstanceOf(IllegalArgumentException.class).hasMessage("lazy boom");
        }
        assertThat(count).hasValue(1);
    }

    @Test
    void shouldRecomputeAfterAnInterruptedComputationAndRestoreTheInterrupt() throws Exception {
        var count = new AtomicInteger();
        var value = new LazyValue<String>(() -> {
            if( count.incrementAndGet() == 1 )
                throw new InterruptedException("lazy interrupt");
            return "ok";
        });

        assertThatThrownBy(value::get)
            .isInstanceOf(IllegalStateException.class)
            .cause().isInstanceOf(InterruptedException.class);
        assertThat(Thread.interrupted()).isTrue();

        assertThat(value.get()).isEqualTo("ok");
        assertThat(value.get()).isEqualTo("ok");
        assertThat(count).hasValue(2);
    }

    @Test
    void shouldRetryInAReaderWaitingForAnInterruptedComputation() throws Exception {
        var count = new AtomicInteger();
        var computing = new CountDownLatch(1);
        var waiting = new CountDownLatch(1);
        var value = new LazyValue<String>(() -> {
            if( count.incrementAndGet() == 1 ) {
                computing.countDown();
                // interrupted, e.g. while a conda build runs
                new CountDownLatch(1).await();
            }
            return "ok";
        });

        var outcome = new AtomicReference<Object>();
        var first = new Thread(() -> {
            try {
                outcome.set(value.get());
            }
            catch( IllegalStateException e ) {
                outcome.set(List.of(e.getCause(), Thread.currentThread().isInterrupted()));
            }
        });
        ExecutorService readers = Executors.newSingleThreadExecutor();
        try {
            first.start();
            computing.await();
            var second = readers.submit(() -> {
                waiting.countDown();
                return value.get();
            });
            waiting.await();
            // let the second reader block on the first computation
            Thread.sleep(100);
            first.interrupt();
            first.join(10_000);

            assertThat((List<?>) outcome.get()).first().isInstanceOf(InterruptedException.class);
            assertThat((List<?>) outcome.get()).last().isEqualTo(true);
            assertThat(second.get(10, TimeUnit.SECONDS)).isEqualTo("ok");
            assertThat(count).hasValue(2);
        }
        finally {
            readers.shutdownNow();
        }
    }
}
