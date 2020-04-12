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
import groovy.transform.PackageScope
import nextflow.exception.IllegalInvocationException
import nextflow.extension.OpCall
import nextflow.extension.OperatorEx

/**
 * Models the execution context of a workflow component
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WorkflowBinding extends Binding  {

    static Map<Object,String> lookupTable = new HashMap<>()

    static String lookup(Object value) {
        return lookupTable.get(value)
    }

    static void init() { lookupTable.clear() }

    private BaseScript owner

    private ScriptMeta meta

    WorkflowBinding() { }

    WorkflowBinding(Map vars) {
        super(vars)
    }

    WorkflowBinding(BaseScript owner) {
        super()
        setOwner(owner)
    }

    WorkflowBinding setOwner(BaseScript owner) {
        this.owner = owner
        this.meta = ScriptMeta.get(owner)
        return this
    }

    BaseScript getOwner() {
        return owner
    }

    @Override
    String toString() {
        "${this.getClass().getSimpleName()}[vars=${variables}]"
    }

    @PackageScope void checkScope0(ComponentDef component) {
        if( component instanceof FunctionDef )
            return // OK
        if( component instanceof ChainableDef && !ExecutionStack.withinWorkflow() ) {
            throw new IllegalInvocationException(component)
        }
    }

    @PackageScope ComponentDef getComponent0(String name) {
        meta.getComponent(name)
    }

    @Override
    Object invokeMethod(String name, Object args) {
        if( meta ) {
            final component = getComponent0(name)
            if( component ) {
                checkScope0(component)
                return component.invoke_o(args)
            }

            // check it's an operator name
            if( OperatorEx.OPERATOR_NAMES.contains(name) )
                return OpCall.create(name, args)
        }

        throw new MissingMethodException(name,this.getClass())
    }

    @Override
    void setVariable(String name, Object value) {
        lookupTable.put(value, name)
        super.setVariable(name, value)
    }

    Object getVariable(String name) {
        try {
            super.getVariable(name)
        }
        catch( MissingPropertyException e ) {
            if( !meta )
                 throw e
            
            def component = getComponent0(name)
            if( component )
                return component

            // check it's an operator name
            if( OperatorEx.OPERATOR_NAMES.contains(name) )
                return OpCall.create(name)

            throw e
        }
    }

}
