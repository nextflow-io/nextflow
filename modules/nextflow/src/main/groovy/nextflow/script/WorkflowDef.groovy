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
import groovy.util.logging.Slf4j
import nextflow.exception.MissingValueException
import nextflow.extension.ChannelFactory
/**
 * Models a script workflow component
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowDef extends BindableDef implements ChainableDef, ExecutionContext {

    private String name

    private BodyDef body

    private List<String> declaredInputs

    private List<String> declaredOutputs

    private Set<String> variableNames

    private BaseScript owner

    // -- following attributes are mutable and instance dependant
    // -- therefore should not be cloned

    private ChannelOut output

    private WorkflowBinding binding

    WorkflowDef(BaseScript owner, Closure<BodyDef> rawBody, String name=null) {
        this.owner = owner
        this.name = name
        // invoke the body resolving in/out params
        final copy = (Closure<BodyDef>)rawBody.clone()
        final resolver = new WorkflowParamsResolver()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.setDelegate(resolver)
        this.body = copy.call()
        // now it can access the parameters
        this.declaredInputs = new ArrayList<>(resolver.getGets().keySet())
        this.declaredOutputs = new ArrayList<>(resolver.getEmits().keySet())
        this.variableNames = getVarNames0()
    }

    /* ONLY FOR TESTING PURPOSE */
    protected WorkflowDef() {}

    WorkflowDef clone() {
        final copy = (WorkflowDef)super.clone()
        copy.@body = body.clone()
        return copy
    }

    WorkflowDef cloneWithName(String name) {
        def result = clone()
        result.@name = name
        return result
    }

    BaseScript getOwner() { owner }

    String getName() { name }

    WorkflowBinding getBinding() { binding }

    ChannelOut getOut() { output }

    @PackageScope BodyDef getBody() { body }

    @PackageScope List<String> getDeclaredInputs() { declaredInputs }

    @PackageScope List<String> getDeclaredOutputs() { declaredOutputs }

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
        final params = ChannelOut.spread(args)
        if( params.size() != declaredInputs.size() )
            throw new IllegalArgumentException("Workflow `$name` declares ${declaredInputs.size()} input channels but ${params.size()} were specified")

        // attach declared inputs with the invocation arguments
        for( int i=0; i< declaredInputs.size(); i++ ) {
            final name = declaredInputs[i]
            context.setProperty( name, params[i] )
        }
    }

    protected ChannelOut collectOutputs(List<String> emissions) {
        final channels = new LinkedHashMap<String, ?>(emissions.size())
        for( String name : emissions ) {
            if( !binding.hasVariable(name) )
                throw new MissingValueException("Missing workflow output parameter: $name")
            final obj = binding.getVariable(name)

            if( ChannelFactory.isChannel(obj) ) {
                channels.put(name, obj)
            }

            else if( obj instanceof ChannelOut ) {
                if( obj.size()>1 )
                    throw new IllegalArgumentException("Cannot emit a multi-channel output: $name")
                if( obj.size()==0 )
                    throw new MissingValueException("Cannot emit empty output: $name")
                channels.put(name, obj.get(0))
            }

            else {
                final value = ChannelFactory.create(true)
                value.bind(obj)
                channels.put(name, value)
            }
        }
        return new ChannelOut(channels)
    }


    Object run(Object[] args) {
        binding = new WorkflowBinding(owner)
        ExecutionStack.push(this)
        try {
            return run0(args)
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
        closure.call()
        // collect the workflow outputs
        return output = collectOutputs(declaredOutputs)
    }

}

/**
 * Hold workflow parameters
 */
@CompileStatic
class WorkflowParamsResolver implements GroovyInterceptable {

    Map<String,Object> gets = new LinkedHashMap<>(10)
    Map<String,Object> emits = new LinkedHashMap<>(10)

    @Override
    def invokeMethod(String name, Object args) {
        if( name.startsWith('_get_') )
            gets.put(name.substring(5), args)
        else if( name.startsWith('_emit_') )
            emits.put(name.substring(6), args)
        else
            throw new IllegalArgumentException("Unknown workflow parameter definition: $name")

    }
}
