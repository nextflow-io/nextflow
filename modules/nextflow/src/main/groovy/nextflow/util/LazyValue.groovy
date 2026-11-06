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
import java.util.concurrent.ExecutionException
import java.util.concurrent.FutureTask

import groovy.transform.CompileStatic
/**
 * A value computed once, in the thread of its first reader. Concurrent readers wait
 * for it, and every reader gets the error of a failed computation. A computation
 * interrupted with an {@link InterruptedException} is not cached: the next reader
 * computes the value again.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class LazyValue<T> {

    private final Callable<T> initializer

    private volatile Computation<T> computation

    LazyValue(Callable<T> initializer) {
        this.initializer = initializer
        this.computation = new Computation<T>(initializer)
    }

    /**
     * @return the value, computing it on the first call
     * @throws IllegalStateException wrapping the error of the computation
     */
    T get() {
        while( true ) {
            final current = computation
            // a no-op once the computation has started
            current.task.run()
            try {
                return current.task.get()
            }
            catch( ExecutionException e ) {
                if( e.cause !instanceof InterruptedException )
                    throw new IllegalStateException(e.cause)
                synchronized (this) {
                    if( computation.is(current) )
                        computation = new Computation<T>(initializer)
                }
                // the interrupted reader fails, the readers waiting for it retry
                if( current.runner.is(Thread.currentThread()) ) {
                    Thread.currentThread().interrupt()
                    throw new IllegalStateException(e.cause)
                }
            }
        }
    }

    /**
     * One attempt at the computation, which records the reader that runs it
     */
    @CompileStatic
    private static class Computation<T> implements Callable<T> {

        private final Callable<T> initializer

        final FutureTask<T> task = new FutureTask<T>(this)

        volatile Thread runner

        Computation(Callable<T> initializer) {
            this.initializer = initializer
        }

        @Override
        T call() {
            runner = Thread.currentThread()
            return initializer.call()
        }
    }
}
