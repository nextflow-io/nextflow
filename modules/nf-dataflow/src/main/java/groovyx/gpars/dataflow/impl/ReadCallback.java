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
package groovyx.gpars.dataflow.impl;

import groovyx.gpars.dataflow.DataflowReadChannel;

/**
 * Receives the item of a non-blocking channel read, see {@link DataflowReadChannel#takeOrSubscribe}.
 * <p>
 * The channel calls it once, without holding a lock, on the thread that binds the item, or on a thread
 * reading the same channel, e.g. a broadcast reader can hand an item to a waiting read while serving a
 * later one. So it must return quickly: hand the item over, then wake or schedule its consumer.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@FunctionalInterface
public interface ReadCallback {

    /**
     * Returned by a read that found no item
     */
    Object NONE = new Object() {
        @Override
        public String toString() {
            return "NONE";
        }
    };

    void onItem(Object item);

    /**
     * Calls the callback, so that a failing callback does not keep the next ones from being called
     *
     * @param failure the error of a previous callback, or null
     * @return the first error, with the next ones suppressed, or null
     */
    static Throwable deliver(final ReadCallback callback, final Object item, final Throwable failure) {
        try {
            callback.onItem(item);
            return failure;
        }
        catch (Throwable e) {
            if (failure == null)
                return e;
            if (failure != e)
                failure.addSuppressed(e);
            return failure;
        }
    }

    /**
     * Rethrows the error collected by {@link #deliver}, if any
     */
    static void rethrow(final Throwable failure) {
        if (failure == null)
            return;
        if (failure instanceof RuntimeException)
            throw (RuntimeException) failure;
        if (failure instanceof Error)
            throw (Error) failure;
        throw new RuntimeException(failure);
    }
}
