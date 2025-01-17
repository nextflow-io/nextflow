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
 */

package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.MissingProcessException
import nextflow.exception.MissingValueException
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
/**
 * Models a script workflow component
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowDef extends BindableDef implements ChainableDef, IterableDef, ExecutionContext {

    private String name

    private BodyDef body

    private List<String> declaredInputs

    private List<String> declaredOutputs

    private Set<String> variableNames

    private BaseScript owner

    // -- following attributes are mutable and instance dependent
    // -- therefore should not be cloned

    private ChannelOut output

    private WorkflowBinding binding

    WorkflowDef(BaseScript owner, Closure<BodyDef> rawBody, String name=null) {
        this.owner = owner
        this.name = name
        // invoke the body resolving in/out params
        final copy = (Closure<BodyDef>)rawBody.clone()
        final resolver = new WorkflowParamsDsl()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.setDelegate(resolver)
        this.body = copy.call()
        // now it can access the parameters
        this.declaredInputs = new ArrayList<>(resolver.getTakes().keySet())
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

    ChannelOut getOut() {
        if( output==null )
            throw new ScriptRuntimeException("Access to '${name}.out' is undefined since the workflow '$name' has not been invoked before accessing the output attribute")
        if( output.size()==0 )
            throw new ScriptRuntimeException("Access to '${name}.out' is undefined since the workflow '$name' doesn't declare any output")

        return output
    }

    @PackageScope BodyDef getBody() { body }

    @PackageScope List<String> getDeclaredInputs() { declaredInputs }

    @PackageScope List<String> getDeclaredOutputs() { declaredOutputs }

    @PackageScope Map<String,Map> getDeclaredPublish() { declaredPublish }

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
        if( params.size() != declaredInputs.size() ) {
            final prefix = name ? "Workflow `$name`" : "Main workflow"
            throw new IllegalArgumentException("$prefix declares ${declaredInputs.size()} input channels but ${params.size()} were given")
        }

        // attach declared inputs with the invocation arguments
        for( int i=0; i< declaredInputs.size(); i++ ) {
            final name = declaredInputs[i]
            context.setProperty( name, params[i] )
        }
    }

    protected ChannelOut collectOutputs(List<String> emissions) {
        // make sure feedback channel cardinality matches
        if( feedbackChannels && feedbackChannels.size() != emissions.size() )
            throw new ScriptRuntimeException("Workflow `$name` inputs and outputs do not have the same cardinality - Feedback loop is not supported"  )

        final channels = new LinkedHashMap<String, DataflowWriteChannel>(emissions.size())
        for( int i=0; i<emissions.size(); i++ ) {
            final targetName = emissions[i]
            if( !binding.hasVariable(targetName) )
                throw new MissingValueException("Missing workflow output parameter: $targetName")
            final obj = binding.getVariable(targetName)

            if( CH.isChannel(obj) ) {
                channels.put(targetName, target(i, obj))
            }

            else if( obj instanceof ChannelOut ) {
                if( obj.size()>1 )
                    throw new IllegalArgumentException("Cannot emit a multi-channel output: $targetName")
                if( obj.size()==0 )
                    throw new MissingValueException("Cannot emit empty output: $targetName")
                channels.put(targetName, target(i, obj.get(0)))
            }

            else {
                if( feedbackChannels!=null )
                    throw new ScriptRuntimeException("Workflow `$name` static output is not allowed when using recursion - Check output: $targetName")
                final value = CH.create(true)
                value.bind(obj)
                channels.put(targetName, value)
            }
        }
        return new ChannelOut(channels)
    }

    protected DataflowWriteChannel target(int index, Object output) {
        if( feedbackChannels==null )
            return (DataflowWriteChannel)output
        // otherwise the output should be forwarded into the feedback channel
        final feedback = feedbackChannels[index]
        CH.getReadChannel(output).into( feedback )
        return feedback
    }

    Object run(Object[] args) {
        binding = new WorkflowBinding(owner)
        ExecutionStack.push(this)
        try {
            return run0(args)
        }
        catch (MissingMethodException e) {
            throw new MissingProcessException(this.binding.scriptMeta, e)
        }
        finally {
            ExecutionStack.pop()
        }
    }

    private Object run0(Object[] args) {
        collectInputs(binding, args)
        // invoke the workflow execution
        final closure = body.closure
        closure.setDelegate(binding)
        closure.setResolveStrategy(Closure.DELEGATE_FIRST)
        closure.call()
        // collect the workflow outputs
        output = collectOutputs(declaredOutputs)
        return output
    }

}

/**
 * Implements the DSL for defining workflow takes and emits
 */
@Slf4j
@CompileStatic
class WorkflowParamsDsl {

    static final private String TAKE_PREFIX = '_take_'
    static final private String EMIT_PREFIX = '_emit_'


    Map<String,Object> takes = new LinkedHashMap<>(10)
    Map<String,Object> emits = new LinkedHashMap<>(10)

    @Override
    def invokeMethod(String name, Object args) {
        if( name.startsWith(TAKE_PREFIX) )
            takes.put(name.substring(TAKE_PREFIX.size()), args)

        else if( name.startsWith(EMIT_PREFIX) )
            emits.put(name.substring(EMIT_PREFIX.size()), args)

        else
            throw new MissingMethodException(name, WorkflowDef, args)
    }
}
