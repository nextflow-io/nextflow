/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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