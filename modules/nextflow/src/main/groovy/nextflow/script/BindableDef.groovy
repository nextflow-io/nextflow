/*
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
import nextflow.exception.DuplicateProcessInvocation
import static nextflow.Const.SCOPE_SEP

/**
 *  Abstract module component which bind itself in the
 *  current execution context once executed
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class BindableDef extends ComponentDef {

    private Set<String> invocations = new HashSet<>()

    abstract Object run(Object[] args)

    Object invoke_a(Object[] args) {

        // use this instance an workflow template, therefore clone it
        final String prefix = ExecutionStack.workflow()?.name
        final fqName = prefix ? prefix+SCOPE_SEP+name : name
        if( this instanceof ProcessDef && !invocations.add(fqName) ) {
            log.debug "Bindable invocations=$invocations"
            final msg = "Process $name has been already used -- If you need to reuse the same component include it with a different name or include in a different workflow context"
            throw new DuplicateProcessInvocation(msg)
        }

        final comp = (prefix ? this.cloneWithName(fqName) : this.clone()) as BindableDef
        // invoke the process execution
        final result = comp.run(args)

        // register this component invocation in the current context
        // so that it can be accessed in the outer execution scope
        // note the simple name (ie. not the one fully qualified with scope prefix) is used here
        // because the component object is associated in the nested scope
        if( name ) {
            final binding = ExecutionStack.binding()
            binding.setVariable(this.name, comp)
        }
        return result
    }

}
