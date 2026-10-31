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
import java.util.concurrent.TimeUnit;

import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.dataflow.impl.Waiter;
import groovyx.gpars.dataflow.operator.ChainWithClosure;
import groovyx.gpars.dataflow.operator.CopyChannelsClosure;

/**
 * An unbounded FIFO channel: each value is read once, by a single reader.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2013 the original author or authors)
 * with its public API and semantics. The buffered values are held in an array deque instead of a
 * dataflow variable each, and a waiting read is handed the next value directly.
 *
 * @param <T> The type of messages to pass through the queue
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings({"unchecked", "rawtypes"})
public class DataflowQueue<T> implements DataflowChannel<T> {

    // stands for a null value, which the deque does not accept
    private static final Object NULL = new Object();

    // guarded by this
    private final ArrayDeque<Object> items = new ArrayDeque<>();

    // the reads waiting for a value, in arrival order, guarded by this
    private ArrayDeque<ReadCallback> waiters;

    @Override
    public final void bind(final T value) {
        final ReadCallback waiter;
        synchronized (this) {
            waiter = waiters != null ? waiters.poll() : null;
            if (waiter == null)
                items.add(value != null ? value : NULL);
        }
        if (waiter != null)
            waiter.onItem(value);
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

    @Override
    public final Object takeOrSubscribe(final ReadCallback callback) {
        synchronized (this) {
            final Object item = items.poll();
            if (item != null)
                return item != NULL ? item : null;
            if (waiters == null)
                waiters = new ArrayDeque<>(2);
            waiters.add(callback);
            return ReadCallback.NONE;
        }
    }

    @Override
    public final synchronized boolean unsubscribe(final ReadCallback callback) {
        return waiters != null && waiters.remove(callback);
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
    public final synchronized boolean isBound() {
        return !items.isEmpty();
    }

    @Override
    public final synchronized int length() {
        return items.size();
    }

    /**
     * Iterates over a snapshot of the buffered values, without consuming them
     */
    public final Iterator<T> iterator() {
        return Collections.unmodifiableList(snapshot()).iterator();
    }

    private synchronized List<T> snapshot() {
        final List<T> result = new ArrayList<>(items.size());
        for (final Object item : items)
            result.add(item != NULL ? (T) item : null);
        return result;
    }

    @Override
    public String toString() {
        return "DataflowQueue(queue=" + snapshot() + ')';
    }
}
