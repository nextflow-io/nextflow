/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
package nextflow.script

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.MissingProcessException

import java.time.Duration
import java.time.Instant

/**
 * Model a testflow DSL definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TestflowDef extends ComponentDef implements ExecutionContext {

    private Closure given
    private Closure when
    private Closure then
    private String name
    private BaseScript owner
    private WorkflowBinding binding
    private Duration runTime

    void given(Closure arg) {
        this.given = arg
    }

    void when(Closure arg) {
        this.when = arg
    }

    void then(Closure arg) {
        this.then = arg
    }

    TestflowDef(BaseScript owner, String name, Closure body) {
        this.owner = owner
        this.name = name
        // the body closure is expected to invokes the `given`, `when`, `then` methods
        // each of them holding a corresponding closure action implementing the
        // testing strategy 
        body.setResolveStrategy(Closure.DELEGATE_ONLY)
        body.setDelegate(this)
        body.call()
    }

    @Override
    String getType() { "testflow" }

    @Override
    String getName() { name }

    Duration getRunTime() { runTime }

    WorkflowBinding getBinding() { binding }

    @Override
    TestflowDef cloneWithName(String name) {
        def result = (TestflowDef) clone()
        result.@name = name
        return result
    }

    @Override
    Object invoke_a(Object[] args) {
        run(args)
        return null
    }

    private void run(Object[] args) {
        binding = new WorkflowBinding(owner)
        ExecutionStack.push(this)
        final runStart = Instant.now()
        try {
            run0(args)
        }
        catch (MissingMethodException e) {
            throw new MissingProcessException(this.binding.scriptMeta, e)
        }
        finally {
            runTime = Duration.between(runStart, Instant.now())
            ExecutionStack.pop()
        }
    }

    private void run0(Object[] args) {
        if (!when)
            throw new IllegalStateException("Missing 'when' statement in testflow: $name")
        if (!then)
            throw new IllegalStateException("Missing 'then' statement in testflow: $name")
        // given part
        if (given) {
            invoke(given, binding)
        }
        // when part
        invoke(when, binding)
        // wrap processes
        // NOTE: this must be invoked here ie. *before* firing the dataflow network (see Session#fireDataflowNetwork)
        // otherwise the DAG is not properly constructed and some output values can be missed
        wrapComponents(binding)
    }


    private void invoke(Closure closure, Object binding) {
        def copy = (Closure) closure.clone()
        copy.setDelegate(binding)
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.call()
    }

    private void wrapComponents(Binding binding) {
        final ctx = binding.variables
        for( String key : ctx.keySet() ) {
            final value = ctx.get(key)
            if( value instanceof ProcessDef ) {
                ctx.put(key, new TestflowDsl(value))
            }
            else if( value instanceof WorkflowDef ) {
                ctx.put(key, new TestflowDsl(value))
            }
        }
    }

    void validateExecution() {
        // then part
        invoke(then, binding)
    }
}
