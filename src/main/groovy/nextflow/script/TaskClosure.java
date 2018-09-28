/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.script;

import groovy.lang.Closure;

/**
 * Extends a {@link Closure} class adding the source attribute
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class TaskClosure extends Closure {

    private Closure delegate;

    private String sourceCode;

    public TaskClosure(Closure code, String source) {
        super(code.getOwner(), code.getThisObject());
        this.delegate = code;
        this.sourceCode = source;
    }

    public String getSource() {
        return sourceCode;
    }

    Closure getInnerClosure() {
        return delegate;
    }

    @Override
    public int getMaximumNumberOfParameters() {
        return delegate.getMaximumNumberOfParameters();
    }

    @Override
    public Class[] getParameterTypes() {
        return delegate.getParameterTypes();
    }

    @Override
    public void setDelegate(final Object delegate) {
        super.setDelegate(delegate);
        this.delegate.setDelegate(delegate);
    }

    @Override
    public void setResolveStrategy(final int resolveStrategy) {
        super.setResolveStrategy(resolveStrategy);
        delegate.setResolveStrategy(resolveStrategy);
    }

    @Override
    public Object call(final Object arguments) {
        return delegate.call(arguments);
    }

    @Override
    public Object call(final Object... args) {
        return delegate.call(args);
    }

    @Override
    public Object call() {
        return delegate.call();
    }

    @Override
    public Object clone() {
        TaskClosure result = (TaskClosure)super.clone();
        result.delegate = (Closure)delegate.clone();
        result.sourceCode = sourceCode;
        return result;
    }

}