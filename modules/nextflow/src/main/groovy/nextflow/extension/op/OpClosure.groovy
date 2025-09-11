/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.extension.op

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.prov.OperatorRun
import nextflow.prov.Prov
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * A closure that wraps the execution of an operator target code
 * associating the inputs and outputs to the corresponding operator run.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class OpClosure extends Closure {

    private final Closure target
    private final OpContext context

    OpClosure(Closure code, OpContext context) {
        super(code.getOwner(), code.getThisObject())
        this.target = code
        this.target.setDelegate(code.getDelegate())
        this.target.setResolveStrategy(code.getResolveStrategy())
        this.context = context
    }

    @Override
    Class<?>[] getParameterTypes() {
        return target.getParameterTypes()
    }

    @Override
    int getMaximumNumberOfParameters() {
        return target.getMaximumNumberOfParameters()
    }

    @Override
    Object getDelegate() {
        return target.getDelegate()
    }

    @Override
    int getDirective() {
        return target.getDirective()
    }

    @Override
    void setDelegate(Object delegate) {
        target.setDelegate(delegate)
    }

    @Override
    void setDirective(int directive) {
        target.setDirective(directive)
    }

    @Override
    void setResolveStrategy(int resolveStrategy) {
        target.setResolveStrategy(resolveStrategy)
    }

    @Override
    void setProperty(String propertyName, Object newValue) {
        target.setProperty(propertyName, newValue)
    }

    @Override
    Object call(final Object... args) {
        // when the accumulator flag true, re-use the previous run object
        final OperatorRun run = context.allocateRun()
        // map the inputs
        final List<Object> inputs = Prov.getTracker().receiveInputs(run, Arrays.asList(args))
        final Object result = InvokerHelper.invokeMethod(target, "call", inputs.toArray())
        // return the operation result
        return result
    }

    @Override
    Object call(Object args) {
        return call(InvokerHelper.asArray(args))
    }

    @Override
    Object call() {
        return call(InvokerHelper.EMPTY_ARGS)
    }

}
