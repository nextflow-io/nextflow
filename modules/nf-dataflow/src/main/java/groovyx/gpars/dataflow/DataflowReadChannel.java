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
// (removed SelectableChannel, rightShift, then, chainWith, filter, split, tap, merge, choice, separate,
// the Pool/PGroup/Map variants and the event manager); added the non-blocking reads of the operator engine.

package groovyx.gpars.dataflow;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import groovy.lang.Closure;
import groovyx.gpars.actor.impl.MessageStream;
import groovyx.gpars.dataflow.impl.ReadCallback;

/**
 * A common interface for all dataflow variables, streams or queues
 *
 * @author Vaclav Pech
 *         Date: 21st Sep 2010
 */
public interface DataflowReadChannel<T> {

    /**
     * Asynchronously retrieves the value from the channel. Sends the actual value of the channel as a message
     * back the the supplied actor once the value has been bound.
     *
     * @param callback The actor to notify when a value is bound
     */
    default void getValAsync(final MessageStream callback) {
        getValAsync(null, callback);
    }

    /**
     * Asynchronously retrieves the value from the channel. Sends a message back the the supplied MessageStream
     * with a map holding the supplied attachment under the 'attachment' key and the actual value of the channel under
     * the 'result' key once the value has been bound.
     *
     * @param attachment arbitrary non-null attachment if reader needs better identification of result
     * @param callback   The actor to notify when a value is bound
     */
    default void getValAsync(final Object attachment, final MessageStream callback) {
        if (callback == null)
            throw new NullPointerException();
        final ReadCallback reader = item -> {
            if (attachment == null) {
                callback.send(item);
                return;
            }
            final Map<String, Object> message = new HashMap<String, Object>(2);
            message.put("attachment", attachment);
            message.put("result", item);
            callback.send(message);
        };
        final Object item = takeOrSubscribe(reader);
        if (item != ReadCallback.NONE)
            reader.onItem(item);
    }

    /**
     * Reads the current value of the channel. Blocks, if the value has not been assigned yet.
     *
     * @return The actual value
     * @throws InterruptedException If the current thread gets interrupted while waiting for the channel to be bound
     */
    T getVal() throws InterruptedException;

    /**
     * Reads the current value of the channel. Blocks up to given timeout, if the value has not been assigned yet.
     *
     * @param timeout The timeout value
     * @param units   Units for the timeout
     * @return The actual value, or null when the timeout expires
     * @throws InterruptedException If the current thread gets interrupted while waiting for the channel to be bound
     */
    T getVal(final long timeout, final TimeUnit units) throws InterruptedException;

    /**
     * Schedules the closure to run on the current parallel group once a value is available.
     * A queue or a stream reader hands each value to a single callback.
     *
     * @param closure closure to execute when data becomes available. The closure should take at most one argument.
     */
    default <V> void whenBound(final Closure<V> closure) {
        getValAsync(new DataCallback(closure, Dataflow.retrieveCurrentDFPGroup()));
    }

    /**
     * Copies all the values of the channel into the target channel
     *
     * @param target The channel to copy the values into
     */
    void into(final DataflowWriteChannel<T> target);

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
