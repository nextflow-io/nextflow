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

/**
 *  Abstract module component which bind itself in the
 *  current execution context once executed
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class BindableDef extends ComponentDef {

    abstract Object run(Object[] args)

    Object invoke_a(Object[] args) {
        // use this instance an workflow template, therefore clone it
        final comp = (BindableDef)this.clone()
        // invoke the process execution
        final result = comp.run(args)
        // register this component invocation in the current context
        // so that it can be accessed in the outer execution scope
        if( name ) {
            final scope = ExecutionStack.binding()
            scope.setVariable(name, comp)
        }
        return result
    }

}
