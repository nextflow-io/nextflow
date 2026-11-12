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
package groovyx.gpars.dataflow.stream;

import java.util.ArrayList;
import java.util.List;

import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import groovyx.gpars.dataflow.impl.ReadCallback;

/**
 * The write end of a stream that every reader reads in full: a singly linked list of nodes, one per
 * value. The writer holds the tail and each {@link DataflowStreamReadAdapter} holds a cursor, so a node
 * is garbage once every reader has read it, and an unreferenced reader is garbage too.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2011 the original author or authors),
 * which wrapped a stream of dataflow variables.
 *
 * @param <T> The type of messages to pass through the stream
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings("unchecked")
public class DataflowStreamWriteAdapter<T> implements DataflowWriteChannel<T> {

    /**
     * A stream node: a value and its position in the stream
     */
    static final class Node {

        final Object item;

        final int seq;

        volatile Node next;

        Node(final Object item, final int seq) {
            this.item = item;
            this.seq = seq;
        }
    }

    // the last node, written under this
    private volatile Node tail = new Node(null, 0);

    // the readers waiting for the next value, guarded by this: a reader or a list of them
    private Object waiters;

    @Override
    public final void bind(final T value) {
        final Object readers;
        synchronized (this) {
            final Node last = tail;
            final Node node = new Node(value, last.seq + 1);
            // the tail first, so that a reader cannot reach a node past the tail
            tail = node;
            last.next = node;
            readers = waiters;
            waiters = null;
        }
        if (readers instanceof DataflowStreamReadAdapter)
            ((DataflowStreamReadAdapter<T>) readers).wake();
        else if (readers != null) {
            // a failing callback of a reader must not keep the next readers from being woken
            Throwable failure = null;
            for (final DataflowStreamReadAdapter<T> reader : (List<DataflowStreamReadAdapter<T>>) readers)
                failure = ReadCallback.deliver(it -> reader.wake(), null, failure);
            ReadCallback.rethrow(failure);
        }
    }

    @Override
    public final DataflowWriteChannel<T> leftShift(final T value) {
        bind(value);
        return this;
    }

    /**
     * Binds the next value of the read channel once available. Unlike GPars, it does not reserve a slot:
     * the value takes its place among the values bound meanwhile when it becomes available, so
     * {@code ch << unboundVar; ch << 1} yields {@code 1} first.
     */
    @Override
    public final DataflowWriteChannel<T> leftShift(final DataflowReadChannel<T> ref) {
        final Object item = ref.takeOrSubscribe(it -> bind((T) it));
        if (item != ReadCallback.NONE)
            bind((T) item);
        return this;
    }

    final Node tail() {
        return tail;
    }

    /**
     * Registers a reader to be woken by the next bind, unless a value follows its cursor already
     *
     * @return false if a value follows the cursor
     */
    final synchronized boolean await(final DataflowStreamReadAdapter<T> reader, final Node cursor) {
        if (cursor.next != null)
            return false;
        if (waiters == null)
            waiters = reader;
        else if (waiters instanceof DataflowStreamReadAdapter) {
            final List<Object> list = new ArrayList<>(4);
            list.add(waiters);
            list.add(reader);
            waiters = list;
        }
        else
            ((List<Object>) waiters).add(reader);
        return true;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + "(written=" + tail.seq + ')';
    }
}
