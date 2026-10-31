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

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.dataflow.Dataflow;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.dataflow.impl.Waiter;
import groovyx.gpars.dataflow.operator.ChainWithClosure;
import groovyx.gpars.dataflow.operator.CopyChannelsClosure;

/**
 * A reader of a stream: it reads every value bound to the stream after its creation, independently
 * of the other readers.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2011 the original author or authors).
 * Unlike GPars, the blocking and the asynchronous reads share a single cursor.
 *
 * @param <T> The type of messages to pass through the stream
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings({"unchecked", "rawtypes"})
public class DataflowStreamReadAdapter<T> implements DataflowReadChannel<T> {

    private final DataflowStreamWriteAdapter<T> stream;

    // the last node read, written under this
    private volatile DataflowStreamWriteAdapter.Node head;

    // the reads waiting for a value, in arrival order, guarded by this
    private ArrayDeque<ReadCallback> pending;

    // whether the reader is registered with the stream for the next bind, guarded by this
    private boolean waiting;

    public DataflowStreamReadAdapter(final DataflowStreamWriteAdapter<T> stream) {
        this.stream = stream;
        this.head = stream.tail();
    }

    @Override
    public final Object takeOrSubscribe(final ReadCallback callback) {
        synchronized (this) {
            if (pending == null || pending.isEmpty()) {
                final DataflowStreamWriteAdapter.Node next = head.next;
                if (next != null) {
                    head = next;
                    return next.item;
                }
            }
            if (pending == null)
                pending = new ArrayDeque<>(2);
            pending.add(callback);
            if (waiting || (waiting = stream.await(this, head)))
                return ReadCallback.NONE;
        }
        // a value was bound meanwhile
        wake();
        return ReadCallback.NONE;
    }

    @Override
    public final synchronized boolean unsubscribe(final ReadCallback callback) {
        return pending != null && pending.remove(callback);
    }

    /**
     * Hands the available values to the waiting reads, then registers with the stream again if some are left
     */
    final void wake() {
        Throwable failure = null;
        for (;;) {
            final ReadCallback callback;
            final Object item;
            synchronized (this) {
                waiting = false;
                if (pending == null || pending.isEmpty())
                    break;
                final DataflowStreamWriteAdapter.Node next = head.next;
                if (next == null) {
                    if (waiting = stream.await(this, head))
                        break;
                    continue;
                }
                head = next;
                item = next.item;
                callback = pending.poll();
            }
            failure = ReadCallback.deliver(callback, item, failure);
        }
        ReadCallback.rethrow(failure);
    }

    @Override
    public final T getVal() throws InterruptedException {
        return (T) Waiter.read(this, -1);
    }

    /**
     * Takes the next value, blocking until there is one or the timeout expires. Returns null on timeout.
     */
    @Override
    public final T getVal(final long timeout, final TimeUnit units) throws InterruptedException {
        final Object result = Waiter.read(this, Math.max(0, units.toNanos(timeout)));
        return result != ReadCallback.NONE ? (T) result : null;
    }

    /**
     * Copies all the values into the target channel, with an operator on the current parallel group
     */
    @Override
    public void into(final DataflowWriteChannel<T> target) {
        Dataflow.retrieveCurrentDFPGroup().operator(this, target, new ChainWithClosure(new CopyChannelsClosure()));
    }

    @Override
    public final boolean isBound() {
        return length() > 0;
    }

    /**
     * The number of values bound to the stream and not read yet. Safe to call from any thread.
     */
    @Override
    public final int length() {
        // the head first, so that it cannot be past the tail
        final int read = head.seq;
        return stream.tail().seq - read;
    }

    /**
     * The number of values not read yet, counting at most {@code max} of them. Safe to call from any thread.
     */
    public final int length(final int max) {
        return Math.min(length(), max);
    }

    /**
     * Iterates over a snapshot of the values not read yet, without consuming them
     */
    public Iterator<T> iterator() {
        final List<T> result = new ArrayList<>();
        // the head first, so that it cannot be past the tail
        DataflowStreamWriteAdapter.Node node = head;
        final DataflowStreamWriteAdapter.Node last = stream.tail();
        while (node != last) {
            node = node.next;
            result.add((T) node.item);
        }
        return Collections.unmodifiableList(result).iterator();
    }

    @Override
    public String toString() {
        return "DataflowStreamReadAdapter(length=" + length() + ')';
    }
}
