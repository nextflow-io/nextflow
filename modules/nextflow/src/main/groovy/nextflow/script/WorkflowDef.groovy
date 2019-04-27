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
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.extension.ChannelFactory
/**
 * Models a script workflow component
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WorkflowDef extends BindableDef implements ChainableDef, ExecutionContext {

    private String name

    private TaskBody body

    private List<String> declaredInputs

    private Set<String> variableNames

    private BaseScript owner

    // -- following attributes are mutable and instance dependant
    // -- therefore should not be cloned

    private output

    private WorkflowBinding binding

    WorkflowDef(BaseScript owner, TaskBody body, String name=null, List<String> inputs = Collections.emptyList() ) {
        this.owner = owner
        this.body = body
        this.name = name
        this.declaredInputs = inputs
        this.variableNames = getVarNames0()
    }

    WorkflowDef clone() {
        final copy = (WorkflowDef)super.clone()
        copy.@body = body.clone()
        return copy
    }

    WorkflowDef withName(String name) {
        def result = clone()
        result.@name = name
        return result
    }

    BaseScript getOwner() { owner }

    String getName() { name }

    WorkflowBinding getBinding() { binding }

    Object getOutput() { output }

    @PackageScope TaskBody getBody() { body }

    @PackageScope List<String> getDeclaredInputs() { declaredInputs }

    @PackageScope String getSource() { body.source }

    @PackageScope List<String> getDeclaredVariables() { new ArrayList<String>(variableNames) }

    String getType() { 'workflow' }

    private Set<String> getVarNames0() {
        def variableNames = body.getValNames()
        if( variableNames ) {
            Set<String> declaredNames = []
            declaredNames.addAll( declaredInputs )
            if( declaredNames )
                variableNames = variableNames - declaredNames
        }
        return variableNames
    }


    protected void collectInputs(Binding context, Object[] args) {
        final params = ChannelArrayList.spread(args)
        if( params.size() != declaredInputs.size() )
            throw new IllegalArgumentException("Workflow `$name` declares ${declaredInputs.size()} input channels but ${params.size()} were specified")

        // attach declared inputs with the invocation arguments
        for( int i=0; i< declaredInputs.size(); i++ ) {
            final name = declaredInputs[i]
            context.setProperty( name, params[i] )
        }
    }

    protected Object collectOutputs(Object output) {
        if( output==null )
            return asChannel(null, true)

        if( output instanceof ChannelArrayList )
            return output

        if( output instanceof DataflowWriteChannel )
            return output

        if( !(output instanceof List) ) {
            return asChannel(output, true)
        }

        def result = asChannel(ChannelArrayList.spread(output))
        if( result.size()==0 )
            return null
        if( result.size()==1 )
            return result[0]
        return result
    }

    private List asChannel(List list) {
        final allScalar = ChannelFactory.allScalar(list)
        for( int i=0; i<list.size(); i++ ) {
            def el = list[i]
            if( !ChannelFactory.isChannel(el) ) {
                list[i] = asChannel(el, allScalar)
            }
        }
        return list
    }

    private DataflowWriteChannel asChannel(Object x, boolean var) {
        if( var ) {
            def result = new DataflowVariable()
            result.bind(x)
            return result
        }

        final result = new DataflowQueue()
        if( x != null ) {
            result.bind(x)
            result.bind(Channel.STOP)
        }
        return result
    }

    
    Object run(Object[] args) {
        binding = new WorkflowBinding(owner)
        ExecutionStack.push(this)
        try {
            run0(args)
        }
        finally {
            ExecutionStack.pop()
        }
    }

    private Object run0(Object[] args) {
        collectInputs(binding, args)
        // invoke the workflow execution
        final closure = body.closure
        closure.delegate = binding
        closure.setResolveStrategy(Closure.DELEGATE_FIRST)
        final result = closure.call()
        // collect the workflow outputs
        output = collectOutputs(result)
    }

}
