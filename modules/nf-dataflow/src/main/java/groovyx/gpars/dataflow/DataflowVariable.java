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

import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import groovyx.gpars.dataflow.expression.DataflowExpression;
import groovyx.gpars.dataflow.impl.ReadCallback;

/**
 * A write-once dataflow variable, the value channel of Nextflow.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2013 the original author or authors)
 * with its public API and semantics, minus remote variables.
 *
 * @param <T> Type of values to bind with the DataflowVariable
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@SuppressWarnings("unchecked")
public class DataflowVariable<T> extends DataflowExpression<T> implements DataflowChannel<T>, Promise<T> {

    public DataflowVariable() {
    }

    /**
     * Binds the value, or the value of a read channel once available
     */
    @Override
    public DataflowWriteChannel<T> leftShift(final T value) {
        if (value instanceof DataflowReadChannel)
            return leftShift((DataflowReadChannel<T>) value);
        bind(value);
        return this;
    }

    /**
     * Binds the next value of the read channel once available
     */
    @Override
    public DataflowWriteChannel<T> leftShift(final DataflowReadChannel<T> ref) {
        final Object item = ref.takeOrSubscribe(it -> bind((T) it));
        if (item != ReadCallback.NONE)
            bind((T) item);
        return this;
    }

    /**
     * Reads the value, blocking until it is bound
     *
     * @throws Throwable the error, if one was bound
     */
    @Override
    @SuppressWarnings({"ProhibitedExceptionDeclared"})
    public final T get() throws Throwable {
        final T result = getVal();
        if (error != null)
            throw error;
        return result;
    }

    /**
     * Reads the value, blocking until it is bound or the timeout expires
     *
     * @throws Throwable the error, if one was bound, or a {@link TimeoutException}
     */
    @Override
    @SuppressWarnings({"ProhibitedExceptionDeclared"})
    public final T get(final long timeout, final TimeUnit units) throws Throwable {
        final T result = getVal(timeout, units);
        if (result instanceof Throwable)
            throw (Throwable) result;
        if (result == null) {
            if (!isBound())
                throw new TimeoutException("Timeout expired in DataflowVariable.get().");
            return get();
        }
        return result;
    }

    @Override
    public boolean isError() {
        return isBound() && error != null;
    }

    @Override
    public Throwable getError() {
        if (isError())
            return error;
        throw new IllegalStateException("No error has been bound to the dataflow variable.");
    }

    @Override
    public void touch() {
    }
}
