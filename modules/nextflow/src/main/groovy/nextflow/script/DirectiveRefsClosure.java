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

package nextflow.script;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import groovy.lang.Closure;

/**
 * Extends a {@link Closure} class carrying the names of the {@code task} directives that
 * it references e.g. {@code memory} for {@code { "-Xmx${task.memory.toGiga()}g" }}.
 *
 * A dynamic directive value defined in the config file is a closure, therefore its source
 * is not available at runtime. The names are instead collected at compile time from the
 * config AST and attached here, so that an executor can tell whether a task depends on a
 * given directive *without* having to evaluate the value.
 *
 * Note this delegates every {@link Closure} method to the closure it wraps, rather than
 * subclassing it, so that the value keeps behaving exactly as the one the user wrote when it
 * is cloned and called by {@code LazyMap#resolveImpl}. The delegation follows the pattern
 * already established by {@link TaskClosure}; the duplication is deliberate, as the two carry
 * unrelated payloads and are attached by different compilation steps.
 *
 * @see nextflow.script.control.DirectiveRefCollector
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class DirectiveRefsClosure extends Closure {

    private Closure delegate;

    private Set<String> directiveRefs;

    public DirectiveRefsClosure(Closure code, List<String> refs) {
        super(code.getOwner(), code.getThisObject());
        this.delegate = code;
        this.directiveRefs = refs != null ? Set.copyOf(refs) : Collections.emptySet();
    }

    public Set<String> getDirectiveRefs() {
        return directiveRefs;
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
        DirectiveRefsClosure result = (DirectiveRefsClosure)super.clone();
        result.delegate = (Closure)delegate.clone();
        result.directiveRefs = directiveRefs;
        return result;
    }

}
