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
package groovyx.gpars.dataflow.expression;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerFieldUpdater;
import java.util.function.BiFunction;

import groovy.lang.Closure;
import groovy.lang.GroovyObject;
import groovy.lang.MetaClass;
import groovy.lang.MetaProperty;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import groovyx.gpars.dataflow.Promise;
import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.dataflow.impl.ThenMessagingRunnable;
import groovyx.gpars.dataflow.impl.Waiter;
import org.codehaus.groovy.runtime.InvokerHelper;

/**
 * A write-once dataflow value: it is bound once, to a value or an error, and read any number of times.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2013 the original author or authors)
 * with its public API and semantics, minus remote binding and bind error listeners.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings({"unchecked", "rawtypes"})
public abstract class DataflowExpression<T> implements DataflowReadChannel<T>, GroovyObject {

    private static final String ALREADY_BOUND = "A DataflowVariable can only be assigned once. Only re-assignments to an equal value are allowed.";

    private static final AtomicIntegerFieldUpdater<DataflowExpression> STATE = AtomicIntegerFieldUpdater.newUpdater(DataflowExpression.class, "state");

    private static final int UNBOUND = 0;
    private static final int BINDING = 1;
    private static final int BOUND = 2;

    private volatile int state;

    // published by the write of state
    protected T value;

    protected Throwable error;

    // the reads waiting for the value, guarded by this
    private List<ReadCallback> waiters;

    // resolved on first use
    private MetaClass metaClass;

    protected DataflowExpression() {
    }

    @Override
    public final boolean isBound() {
        return state == BOUND;
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
        if (state != BOUND) {
            synchronized (this) {
                if (state != BOUND) {
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
        return error != null ? error : value;
    }

    /**
     * Reads the value, blocking until it is bound. Returns null if an error was bound.
     */
    @Override
    public T getVal() throws InterruptedException {
        if (state != BOUND)
            Waiter.read(this, -1);
        return value;
    }

    /**
     * Reads the value, blocking until it is bound or the timeout expires. Returns null on timeout or if an error was bound.
     */
    @Override
    public T getVal(final long timeout, final TimeUnit units) throws InterruptedException {
        if (state != BOUND && Waiter.read(this, Math.max(0, units.toNanos(timeout))) == ReadCallback.NONE)
            return null;
        return value;
    }

    public final void join() throws InterruptedException {
        getVal();
    }

    public final void join(final long timeout, final TimeUnit units) throws InterruptedException {
        getVal(timeout, units);
    }

    /**
     * Binds the value. Binding again an equal value is a no-op, binding a different one fails.
     * Binding a {@link Throwable} also makes it the error of the expression.
     */
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
     * Binds the value, failing if the expression is bound already
     */
    public final void bindUnique(final T value) {
        if (!STATE.compareAndSet(this, UNBOUND, BINDING))
            throw new IllegalStateException("A DataflowVariable can only be assigned once. Use bind() to allow for equal values to be passed into already-bound variables.");
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
        this.value = value;
        this.error = error;
        final List<ReadCallback> callbacks;
        synchronized (this) {
            state = BOUND;
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

    @Override
    public final <V> void whenBound(final Closure<V> closure) {
        DataflowReadChannel.super.whenBound(closure);
    }

    /**
     * Schedules the closure once the value is bound, and binds its result to the returned variable.
     * A bound error is propagated without running the closure.
     */
    public final <V> Promise<V> then(final Closure<V> closure) {
        return then(closure, null);
    }

    public <V> Promise<V> then(final Closure<V> closure, final Closure<V> errorHandler) {
        final DataflowVariable<V> result = new DataflowVariable<V>();
        whenBound(new ThenMessagingRunnable<T, V>(result, closure, errorHandler));
        return result;
    }

    /**
     * Like {@link #then(Closure)}, returning the result as a read channel
     */
    public final <V> DataflowReadChannel<V> chainWith(final Closure<V> closure) {
        return (DataflowReadChannel<V>) then(closure);
    }

    /**
     * Copies the value into the target channel once bound
     */
    @Override
    public void into(final DataflowWriteChannel<T> target) {
        final Object item = takeOrSubscribe(it -> target.bind((T) it));
        if (item != ReadCallback.NONE)
            target.bind((T) item);
    }

    @Override
    public MetaClass getMetaClass() {
        MetaClass result = metaClass;
        if (result == null)
            metaClass = result = InvokerHelper.getMetaClass(getClass());
        return result;
    }

    @Override
    public void setMetaClass(final MetaClass metaClass) {
        this.metaClass = metaClass;
    }

    /**
     * Returns a property of the expression or, for an unknown property, a variable bound to the property of the value
     */
    @Override
    public Object getProperty(final String name) {
        final MetaProperty property = getMetaClass().hasProperty(this, name);
        if (property != null)
            return property.getProperty(this);
        return derive(new Object[0], (receiver, args) -> InvokerHelper.getProperty(receiver, name));
    }

    /**
     * Invokes a method of the expression or, for an unknown method, returns a variable bound to the result of
     * the method on the value. Dataflow arguments are replaced by their values.
     */
    @Override
    public Object invokeMethod(final String name, final Object args) {
        if (!getMetaClass().respondsTo(this, name).isEmpty())
            return getMetaClass().invokeMethod(this, name, args);
        return derive(InvokerHelper.asArray(args), (receiver, values) -> InvokerHelper.invokeMethod(receiver, name, values));
    }

    /**
     * Returns a variable bound to the result of the function once the value and the dataflow arguments are bound,
     * or to the first of their errors. The function runs on the thread that binds the last of them.
     */
    private DataflowVariable<Object> derive(final Object[] args, final BiFunction<Object, Object[], Object> function) {
        final DataflowVariable<Object> result = new DataflowVariable<>();
        final List<DataflowExpression<?>> inputs = new ArrayList<>(1 + args.length);
        inputs.add(this);
        for (final Object arg : args) {
            if (arg instanceof DataflowExpression)
                inputs.add((DataflowExpression<?>) arg);
        }
        final AtomicInteger pending = new AtomicInteger(inputs.size());
        final ReadCallback callback = ignored -> {
            if (pending.decrementAndGet() == 0)
                evaluate(args, function, result);
        };
        for (final DataflowExpression<?> input : inputs) {
            if (input.takeOrSubscribe(callback) != ReadCallback.NONE)
                callback.onItem(null);
        }
        return result;
    }

    private void evaluate(final Object[] args, final BiFunction<Object, Object[], Object> function, final DataflowVariable<Object> result) {
        Throwable failure = error;
        final Object[] values = args.clone();
        for (int i = 0; i < values.length; i++) {
            if (values[i] instanceof DataflowExpression) {
                final DataflowExpression<?> arg = (DataflowExpression<?>) values[i];
                if (failure == null)
                    failure = arg.error;
                values[i] = arg.value;
            }
        }
        Object output = null;
        if (failure == null) {
            try {
                output = function.apply(value, values);
            }
            catch (Throwable e) {
                failure = e;
            }
        }
        if (failure != null)
            result.bindError(failure);
        else
            result.bind(output);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + "(value=" + value + ')';
    }
}
