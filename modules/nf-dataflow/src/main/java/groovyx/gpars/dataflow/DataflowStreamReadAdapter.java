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
package groovyx.gpars.dataflow;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.dataflow.impl.Waiter;
import groovyx.gpars.dataflow.operator.ChainWithClosure;
import groovyx.gpars.dataflow.operator.CopyChannelsClosure;
import groovyx.gpars.dataflow.operator.DataflowOperator;
import groovyx.gpars.dataflow.operator.DataflowProcessor;

/**
 * A reader of a stream: it reads every value bound to the stream after its creation, independently
 * of the other readers.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2011 the original author or authors).
 * Unlike GPars, the blocking and the asynchronous reads share a single cursor, and there is no timed read.
 *
 * @param <T> The type of messages to pass through the stream
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings({"unchecked", "rawtypes"})
public class DataflowStreamReadAdapter<T> implements DataflowReadChannel<T> {

    private final DataflowBroadcast<T> stream;

    // the last node read, written under this
    private volatile DataflowBroadcast.Node headNode;

    // the reads waiting for a value, in arrival order, guarded by this
    private ArrayDeque<ReadCallback> pending;

    // whether the reader is registered with the stream for the next bind, guarded by this
    private boolean waiting;

    public DataflowStreamReadAdapter(final DataflowBroadcast<T> stream) {
        this.stream = stream;
        this.headNode = stream.tail();
    }

    @Override
    public final Object takeOrSubscribe(final ReadCallback callback) {
        synchronized (this) {
            if (pending == null || pending.isEmpty()) {
                final DataflowBroadcast.Node next = headNode.next;
                if (next != null) {
                    headNode = next;
                    return next.item;
                }
            }
            if (pending == null)
                pending = new ArrayDeque<>(2);
            pending.add(callback);
            if (waiting || (waiting = stream.await(this, headNode)))
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
                final DataflowBroadcast.Node next = headNode.next;
                if (next == null) {
                    if (waiting = stream.await(this, headNode))
                        break;
                    continue;
                }
                headNode = next;
                item = next.item;
                callback = pending.poll();
            }
            failure = ReadCallback.deliver(callback, item, failure);
        }
        ReadCallback.rethrow(failure);
    }

    @Override
    public final T getVal() throws InterruptedException {
        return (T) Waiter.read(this);
    }

    /**
     * Copies all the values into the target channel, with an operator
     */
    @Override
    public void into(final DataflowWriteChannel<T> target) {
        final Map<String, Object> params = Map.of(DataflowProcessor.INPUTS, List.of(this), DataflowProcessor.OUTPUTS, List.of(target));
        new DataflowOperator(params, new ChainWithClosure(new CopyChannelsClosure())).start();
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
        final int read = headNode.seq;
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
        DataflowBroadcast.Node node = headNode;
        final DataflowBroadcast.Node last = stream.tail();
        // the tail is published before it is linked, so the last link may not be set yet
        while (node != last && node.next != null) {
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
