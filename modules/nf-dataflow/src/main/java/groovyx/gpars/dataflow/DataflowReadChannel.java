// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2013  The original author or authors
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
// Modified by Seqera Labs for Nextflow (2026): trimmed to the methods used by Nextflow and its plugins
// (removed SelectableChannel, rightShift, then, chainWith, into, filter, split, tap, merge, choice, separate,
// the Pool/PGroup/Map variants, getValAsync, the timed getVal and the event manager); added the non-blocking reads of the operator engine.

package groovyx.gpars.dataflow;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.impl.ReadCallback;

/**
 * A common interface for all dataflow variables, streams or queues
 *
 * @author Vaclav Pech
 *         Date: 21st Sep 2010
 */
public interface DataflowReadChannel<T> {

    /**
     * Reads the current value of the channel. Blocks, if the value has not been assigned yet.
     *
     * @return The actual value
     * @throws InterruptedException If the current thread gets interrupted while waiting for the channel to be bound
     */
    T getVal() throws InterruptedException;

    /**
     * Schedules the closure to run on the dataflow pool once a value is available.
     * A queue or a stream reader hands each value to a single callback.
     *
     * @param closure closure to execute when data becomes available. The closure should take at most one argument.
     */
    default <V> void whenBound(final Closure<V> closure) {
        final DataCallback callback = new DataCallback(closure);
        final Object item = takeOrSubscribe(callback);
        if (item != ReadCallback.NONE)
            callback.onItem(item);
    }

    /**
     * Check if value has been set already for this expression
     *
     * @return true if bound already
     */
    boolean isBound();

    /**
     * Reports the current number of elements in the channel
     *
     * @return The current snapshot of the number of elements in the channel
     */
    int length();

    /**
     * Takes the next value without blocking. When there is none, subscribes the callback, which is then
     * called once with the next value. Engine-internal: used by the operator engine and the blocking reads.
     *
     * @return the value, or {@link ReadCallback#NONE} when the callback was subscribed
     */
    Object takeOrSubscribe(final ReadCallback callback);

    /**
     * Cancels a subscription of {@link #takeOrSubscribe}. Engine-internal.
     *
     * @return true if cancelled, false if the callback has been or is being called
     */
    boolean unsubscribe(final ReadCallback callback);
}
