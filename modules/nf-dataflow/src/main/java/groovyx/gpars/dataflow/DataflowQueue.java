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
 * An unbounded FIFO channel: each value is read once, by a single reader.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2013 the original author or authors)
 * with its public API and semantics, minus the timed read and leftShift. The buffered values are held in
 * an array deque instead of a dataflow variable each, and a waiting read is handed the next value directly.
 *
 * @param <T> The type of messages to pass through the queue
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings({"unchecked", "rawtypes"})
public class DataflowQueue<T> implements DataflowReadChannel<T>, DataflowWriteChannel<T> {

    // stands for a null value, which the deque does not accept
    private static final Object NULL = new Object();

    // guarded by this
    private final ArrayDeque<Object> queued = new ArrayDeque<>();

    // the reads waiting for a value, in arrival order, guarded by this
    private ArrayDeque<ReadCallback> waiters;

    @Override
    public final void bind(final T value) {
        final ReadCallback waiter;
        synchronized (this) {
            waiter = waiters != null ? waiters.poll() : null;
            if (waiter == null)
                queued.add(value != null ? value : NULL);
        }
        if (waiter != null)
            waiter.onItem(value);
    }

    @Override
    public final Object takeOrSubscribe(final ReadCallback callback) {
        synchronized (this) {
            final Object item = queued.poll();
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
    public final synchronized boolean isBound() {
        return !queued.isEmpty();
    }

    @Override
    public final synchronized int length() {
        return queued.size();
    }

    /**
     * Iterates over a snapshot of the buffered values, without consuming them
     */
    public final Iterator<T> iterator() {
        return Collections.unmodifiableList(snapshot()).iterator();
    }

    private synchronized List<T> snapshot() {
        final List<T> result = new ArrayList<>(queued.size());
        for (final Object item : queued)
            result.add(item != NULL ? (T) item : null);
        return result;
    }

    @Override
    public String toString() {
        return "DataflowQueue(queue=" + snapshot() + ')';
    }
}
