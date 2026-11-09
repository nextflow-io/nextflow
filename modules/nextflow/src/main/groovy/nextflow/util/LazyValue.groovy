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
 * for it, and every reader gets the error of a failed computation.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class LazyValue<T> {

    private final FutureTask<T> task

    LazyValue(Callable<T> initializer) {
        this.task = new FutureTask<T>(initializer)
    }

    /**
     * @return the value, computing it on the first call
     * @throws IllegalStateException wrapping the error of the computation, or the
     *      {@link InterruptedException} of a reader interrupted while waiting for it
     */
    T get() {
        // a no-op once the computation has started
        task.run()
        try {
            return task.get()
        }
        catch( ExecutionException e ) {
            throw new IllegalStateException(e.cause)
        }
        catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
            throw new IllegalStateException(e)
        }
    }
}
