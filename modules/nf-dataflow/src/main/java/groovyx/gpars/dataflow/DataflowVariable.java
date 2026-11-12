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

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicIntegerFieldUpdater;

import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.dataflow.impl.Waiter;

/**
 * A write-once dataflow variable, the value channel of Nextflow: it is bound once, to a value or an error,
 * and read any number of times.
 * <p>
 * Rewrites the GPars 1.2.1 classes {@code DataflowVariable} and {@code DataflowExpression}
 * (Copyright 2008-2013 the original author or authors) with their public API and semantics, minus remote
 * variables, bind error listeners, bindUnique, getValAsync, the timed reads, leftShift, then, chainWith, into
 * and the derived expressions for unknown properties and methods.
 *
 * @param <T> Type of values to bind with the DataflowVariable
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings({"unchecked", "rawtypes"})
public class DataflowVariable<T> implements DataflowReadChannel<T>, DataflowWriteChannel<T> {

    private static final String ALREADY_BOUND = "A DataflowVariable can only be assigned once. Only re-assignments to an equal value are allowed.";

    private static final AtomicIntegerFieldUpdater<DataflowVariable> STATE = AtomicIntegerFieldUpdater.newUpdater(DataflowVariable.class, "bindState");

    private static final int UNBOUND = 0;
    private static final int BINDING = 1;
    private static final int BOUND = 2;

    // fields are named so that Groovy property access (e.g. `.value`) does not resolve to them
    private volatile int bindState;

    // published by the write of bindState
    private T boundValue;

    private Throwable error;

    // the reads waiting for the value, guarded by this
    private List<ReadCallback> waiters;

    public DataflowVariable() {
    }

    @Override
    public final boolean isBound() {
        return bindState == BOUND;
    }

    @Override
    public final int length() {
        return isBound() ? 1 : 0;
    }

    /**
     * Returns the value, or the error if one was bound, without consuming it
     */
    @Override
    public Object takeOrSubscribe(final ReadCallback callback) {
        if (bindState != BOUND) {
            synchronized (this) {
                if (bindState != BOUND) {
                    if (waiters == null)
                        waiters = new ArrayList<>(2);
                    waiters.add(callback);
                    return ReadCallback.NONE;
                }
            }
        }
        return result();
    }

    @Override
    public final synchronized boolean unsubscribe(final ReadCallback callback) {
        return waiters != null && waiters.remove(callback);
    }

    private Object result() {
        return error != null ? error : boundValue;
    }

    /**
     * Reads the value, blocking until it is bound. Returns null if an error was bound.
     */
    @Override
    public T getVal() throws InterruptedException {
        if (bindState != BOUND)
            Waiter.read(this);
        return boundValue;
    }

    /**
     * Reads the value, blocking until it is bound
     *
     * @throws Throwable the error, if one was bound
     */
    @SuppressWarnings({"ProhibitedExceptionDeclared"})
    public final T get() throws Throwable {
        final T result = getVal();
        if (error != null)
            throw error;
        return result;
    }

    /**
     * Binds the value. Binding again an equal value is a no-op, binding a different one fails.
     * Binding a {@link Throwable} also makes it the error of the variable.
     */
    @Override
    public final void bind(final T value) {
        if (!STATE.compareAndSet(this, UNBOUND, BINDING)) {
            try {
                if (Objects.equals(value, getVal()))
                    return;
            }
            catch (InterruptedException ignored) {
                Thread.currentThread().interrupt();
            }
            throw new IllegalStateException(ALREADY_BOUND);
        }
        doBind(value, value instanceof Throwable ? (Throwable) value : null);
    }

    /**
     * Binds an error: the value reads as null, while {@code whenBound} callbacks receive the error
     */
    public final void bindError(final Throwable e) {
        if (!STATE.compareAndSet(this, UNBOUND, BINDING))
            throw new IllegalStateException(ALREADY_BOUND);
        doBind(null, e);
    }

    private void doBind(final T value, final Throwable error) {
        this.boundValue = value;
        this.error = error;
        final List<ReadCallback> callbacks;
        synchronized (this) {
            bindState = BOUND;
            callbacks = waiters;
            waiters = null;
        }
        if (callbacks == null)
            return;
        final Object result = result();
        Throwable failure = null;
        for (final ReadCallback callback : callbacks)
            failure = ReadCallback.deliver(callback, result, failure);
        ReadCallback.rethrow(failure);
    }

    public boolean isError() {
        return isBound() && error != null;
    }

    public Throwable getError() {
        if (isError())
            return error;
        throw new IllegalStateException("No error has been bound to the dataflow variable.");
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + "(value=" + boundValue + ')';
    }
}
