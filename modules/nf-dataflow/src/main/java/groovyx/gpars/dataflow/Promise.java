// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2012  The original author or authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Modified by Seqera Labs for Nextflow (2026): trimmed to get, join, whenBound, then and the error accessors;
// no longer a SelectableChannel.

package groovyx.gpars.dataflow;

import groovy.lang.Closure;

import java.util.concurrent.TimeUnit;

/**
 * A unifying future-like interface to dataflow variables, asynchronous functions and active objects.
 * Represents the read end of DataflowVariables.
 *
 * @author Vaclav Pech
 */
public interface Promise<T> {

    /**
     * Retrieves the value of the variable, blocking until a value is available
     *
     * @return The value stored in the variable
     * @throws Throwable If the stored value is an exception instance it gets re-thrown
     */
    @SuppressWarnings({"ProhibitedExceptionDeclared"})
    T get() throws Throwable;

    /**
     * Retrieves the value of the variable, blocking up to given timeout, if the value has not been assigned yet.
     *
     * @param timeout The timeout value
     * @param units   Units for the timeout
     * @return The value stored in the variable
     * @throws Throwable If the stored value is an exception instance it gets re-thrown
     */
    @SuppressWarnings({"ProhibitedExceptionDeclared"})
    T get(final long timeout, final TimeUnit units) throws Throwable;

    /**
     * Blocks until a value gets bound
     *
     * @throws InterruptedException If the current thread gets interrupted while waiting for the value
     */
    void join() throws InterruptedException;

    /**
     * Blocks until a value gets bound, but no longer than the given timeout
     *
     * @param timeout The timeout value
     * @param units   Units for the timeout
     * @throws InterruptedException If the current thread gets interrupted while waiting for the value
     */
    void join(final long timeout, final TimeUnit units) throws InterruptedException;

    /**
     * Schedules the closure to run on the current parallel group once the value is bound
     *
     * @param closure closure to execute when data becomes available. The closure should take at most one argument.
     */
    <V> void whenBound(final Closure<V> closure);

    /**
     * Schedules the closure to run once the value is bound, and binds its result to the returned promise.
     * A bound error is propagated to the returned promise without running the closure.
     *
     * @param closure closure to execute when data becomes available. The closure should take at most one argument.
     * @return A promise for the results of the supplied closure
     */
    <V> Promise<V> then(final Closure<V> closure);

    /**
     * Schedules the closure to run once the value is bound, or the error handler once an error is bound,
     * and binds the result to the returned promise.
     *
     * @param closure      closure to execute when data becomes available. The closure should take at most one argument.
     * @param errorHandler closure to execute when an error (instance of Throwable) gets bound. The closure should take at most one argument.
     * @return A promise for the results of the supplied closure
     */
    <V> Promise<V> then(final Closure<V> closure, final Closure<V> errorHandler);

    /**
     * Check if value has been set already for this expression
     *
     * @return true if bound already
     */
    boolean isBound();

    /**
     * Checks if the promise is bound to an error
     *
     * @return True, if an error has been bound
     */
    boolean isError();

    /**
     * Returns the error bound to the promise
     *
     * @return The error
     * @throws IllegalStateException If not bound or not bound to an error
     */
    Throwable getError();

    /**
     * May be used by lazy implementations to warm up
     */
    void touch();
}
