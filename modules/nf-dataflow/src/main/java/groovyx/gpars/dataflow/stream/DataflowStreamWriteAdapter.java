// GPars - Groovy Parallel Systems
//
// Copyright © 2008-11  The original author or authors
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

package groovyx.gpars.dataflow.stream;

import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;

/**
 * Adapts a DataflowStream to accommodate for the DataflowWriteChannel interface.
 * It also synchronizes all changes to the underlying stream allowing multiple threads accessing the stream concurrently.
 *
 * @param <T> The type of messages to pass through the stream
 * @author Vaclav Pech
 */
@SuppressWarnings({"SynchronizedMethod"})
public class DataflowStreamWriteAdapter<T> implements DataflowWriteChannel<T> {

    private StreamCore<T> head;

    /**
     * Creates a new adapter
     *
     * @param stream The stream to wrap
     */
    public DataflowStreamWriteAdapter(final StreamCore<T> stream) {
        this.head = stream;
    }

    @Override
    public final DataflowWriteChannel<T> leftShift(final T value) {
        updateHead().leftShift(value);
        return this;
    }

    @Override
    public final DataflowWriteChannel<T> leftShift(final DataflowReadChannel<T> ref) {
        updateHead().leftShift(ref);
        return this;
    }

    @Override
    public final void bind(final T value) {
        updateHead().leftShift(value);
    }

    /**
     * Moves head
     *
     * @return The old head
     */
    private synchronized StreamCore<T> updateHead() {
        final StreamCore<T> oldHead = head;
        head = (StreamCore<T>) head.getRest();
        return oldHead;
    }

    @Override
    public synchronized String toString() {
        return head.toString();
    }

    protected final synchronized StreamCore<T> getHead() {
        return head;
    }
}

