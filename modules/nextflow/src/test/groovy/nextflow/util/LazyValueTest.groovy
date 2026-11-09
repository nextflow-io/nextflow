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

package nextflow.util

import java.util.concurrent.Callable
import java.util.concurrent.CountDownLatch
import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicReference

import spock.lang.Specification
import spock.lang.Timeout

@Timeout(value = 30, unit = TimeUnit.SECONDS)
class LazyValueTest extends Specification {

    /**
     * A closure is both a Runnable and a Callable, so wrap it to submit it as a Callable
     */
    static <T> Callable<T> callable(Closure<T> action) {
        return new Callable<T>() {
            @Override
            T call() throws Exception { action.call() }
        }
    }

    def 'should compute once in the thread of the first reader' () {
        given:
        def count = new AtomicInteger()
        def computing = new CountDownLatch(1)
        def release = new CountDownLatch(1)
        def value = new LazyValue<String>({
            count.incrementAndGet()
            computing.countDown()
            release.await()
            return Thread.currentThread().getName()
        })
        def readers = Executors.newFixedThreadPool(8)

        when:
        List<Future<String>> results = []
        results.add(readers.submit(callable { value.get() }))
        computing.await()
        for( int i = 1; i < 8; i++ )
            results.add(readers.submit(callable { value.get() }))
        release.countDown()
        def first = results[0].get(10, TimeUnit.SECONDS)

        then:
        results.every { it.get(10, TimeUnit.SECONDS) == first }
        first.startsWith('pool-')
        count.get() == 1

        cleanup:
        readers.shutdownNow()
    }

    def 'should give the error to every reader' () {
        given:
        def count = new AtomicInteger()
        def value = new LazyValue<String>({
            count.incrementAndGet()
            throw new IllegalArgumentException('lazy boom')
        })

        when:
        def errors = (1..3).collect {
            try { value.get(); return null }
            catch( IllegalStateException e ) { return e }
        }

        then:
        errors.size() == 3
        errors.every { it instanceof IllegalStateException }
        errors.every { it.cause instanceof IllegalArgumentException }
        errors.every { it.cause.message == 'lazy boom' }
        count.get() == 1
    }

    def 'should cache an interrupted exception thrown by the computation' () {
        given:
        // e.g. a FileMutex lock timeout, which is not a real interrupt
        def count = new AtomicInteger()
        def value = new LazyValue<String>({
            count.incrementAndGet()
            throw new InterruptedException('lock timeout')
        })

        when:
        def errors = (1..3).collect {
            try { value.get(); return null }
            catch( IllegalStateException e ) { return e }
        }

        then:
        errors.every { it.cause instanceof InterruptedException }
        !Thread.currentThread().isInterrupted()
        count.get() == 1
    }

    def 'should fail a waiting reader that is interrupted without affecting the computation' () {
        given:
        def count = new AtomicInteger()
        def computing = new CountDownLatch(1)
        def release = new CountDownLatch(1)
        def value = new LazyValue<String>({
            count.incrementAndGet()
            computing.countDown()
            release.await()
            return 'ok'
        })
        def outcome = new AtomicReference<Object>()
        def readers = Executors.newSingleThreadExecutor()
        def waiter = new Thread({
            try {
                outcome.set(value.get())
            }
            catch( IllegalStateException e ) {
                outcome.set([e.cause, Thread.currentThread().isInterrupted()])
            }
        } as Runnable)

        when:
        def first = readers.submit(callable { value.get() })
        computing.await()
        waiter.start()
        // let the waiter block on the computation
        Thread.sleep(100)
        waiter.interrupt()
        waiter.join(10_000)
        release.countDown()

        then:
        (outcome.get() as List).first() instanceof InterruptedException
        (outcome.get() as List).last() == true
        first.get(10, TimeUnit.SECONDS) == 'ok'
        value.get() == 'ok'
        count.get() == 1

        cleanup:
        readers.shutdownNow()
    }
}
