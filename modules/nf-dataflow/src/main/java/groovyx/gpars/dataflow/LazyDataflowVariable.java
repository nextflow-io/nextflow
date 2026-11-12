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

import java.util.concurrent.atomic.AtomicBoolean;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.group.PGroup;

/**
 * A dataflow variable computed by its initializer, on the group's pool, on the first read.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2012 the original author or authors).
 *
 * @param <T> Type of values to bind with the DataflowVariable
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public final class LazyDataflowVariable<T> extends DataflowVariable<T> {

    private final Closure<T> initializer;

    private final PGroup group;

    private final AtomicBoolean started = new AtomicBoolean();

    public LazyDataflowVariable(final Closure<T> initializer) {
        this(Dataflow.DATA_FLOW_GROUP, initializer);
    }

    public LazyDataflowVariable(final PGroup group, final Closure<T> initializer) {
        if (initializer == null)
            throw new IllegalArgumentException("The LazyDataflowVariable initializer must not be null.");
        if (group == null)
            throw new IllegalArgumentException("The LazyDataflowVariable initialization group must not be null.");
        this.group = group;
        this.initializer = initializer;
    }

    private void initialize() {
        if (started.get() || !started.compareAndSet(false, true))
            return;
        group.getThreadPool().execute(() -> {
            try {
                leftShift(initializer.call());
            }
            catch (Throwable e) {
                bindError(e);
            }
        });
    }

    // every read goes through here: getVal, whenBound, then
    @Override
    public Object takeOrSubscribe(final ReadCallback callback) {
        initialize();
        return super.takeOrSubscribe(callback);
    }

    @Override
    public boolean isError() {
        initialize();
        return super.isError();
    }

    @Override
    public Throwable getError() {
        initialize();
        return super.getError();
    }

    @Override
    public void touch() {
        initialize();
    }
}
