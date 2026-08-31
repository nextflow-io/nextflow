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

package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.dataflow.ChannelImpl
import nextflow.dataflow.ValueImpl
import nextflow.exception.IllegalInvocationException
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.OpCall
import nextflow.util.RecordMap
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Models the execution context of a workflow component
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
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

    ScriptMeta getScriptMeta() { meta }

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

    /**
     * Invokes custom methods in the task execution context
     *
     * @see  BaseScript#invokeMethod(java.lang.String, java.lang.Object)
     *
     * @param name the name of the method to call
     * @param args the arguments to use for the method call
     * @return The result of the method invocation
     */
    @Override
    Object invokeMethod(String name, Object args) {
        if( meta ) {
            log.trace "Trying to invoke component: $name - args=${renderArgs(args)}"
            final component = getComponent0(name)
            if( component ) {
                checkScope0(component)
                return invoke0(component, args)
            }

            // check it's an operator name
            if( !owner?.isTypingEnabled() && NF.hasOperator(name) )
                return OpCall.create(name, args)
        }

        throw new MissingMethodException(name,this.getClass())
    }

    private static String renderArgs(Object args) {
        try {
            return InvokerHelper.toString(args)
        }
        catch( Throwable t ) {
            return "<unable to render args: ${t.class.simpleName}>"
        }
    }

    private Object invoke0(ComponentDef component, Object args) {
        final componentTyped = isTypedComponent(component)
        final args1 = DataflowTypeHelper.normalizeArray(args, componentTyped)
        final result = component.invoke_a(args1)
        return DataflowTypeHelper.normalize(result, owner?.isTypingEnabled())
    }

    private static boolean isTypedComponent(ComponentDef component) {
        if( component instanceof WorkflowDef )
            return component.getOwner().isTypingEnabled()
        if( component instanceof PipelineDef )
            return component.getOwner().isTypingEnabled()
        return false
    }

    @Override
    void setProperty(String name, Object value) {
        // when a script variable name matches a BaseScript attribute name
        // set it directly via `setVariable` otherwise groovy will try to
        // set into the base script attribute
        if( name in BaseScriptConsts.PRIVATE_NAMES )
            setVariable(name, value)
        else
            super.setProperty(name, value)
    }

    @Override
    Object getProperty(String name) {
        if( name in BaseScriptConsts.PRIVATE_NAMES )
            return getVariable(name)
        else
            return super.getProperty(name)
    }

    @Override
    void setVariable(String name, Object value) {
        lookupTable.put(value, name)
        setVariable0(name, value)
    }

    protected void setVariable0(String name, Object value) {
        super.setVariable(name, value)
    }

    Object getVariable(String name) {
        try {
            super.getVariable(name)
        }
        catch( MissingPropertyException e ) {
            if( owner?.isTypingEnabled() )
                throw e

            if( !meta )
                throw e

            def component = getComponent0(name)
            if( component )
                return component

            // check it's an operator name
            if( NF.hasOperator(name) )
                return OpCall.create(name)

            throw e
        }
    }

    void _publish_(String name, Object source) {
        if( source instanceof ChannelOut ) {
            if( source.size() > 1 )
                throw new ScriptRuntimeException("Cannot assign a multi-channel to a workflow output: $name")
            source = source[0]
        }

        // a record of channels -- e.g. the outputs of an included pipeline --
        // publishes each of its fields as a separate output
        if( isRecordOfChannels(source) ) {
            for( final entry : (Map<String,Object>)source )
                publish0(entry.key, entry.value)
            return
        }

        publish0(name, source)
    }

    private static boolean isRecordOfChannels(Object source) {
        if( source !instanceof RecordMap || !source )
            return false
        return ((Map)source).values().every { value ->
            value instanceof ChannelImpl || value instanceof ValueImpl || value instanceof DataflowWriteChannel
        }
    }

    private void publish0(String name, Object source) {
        if( owner.session.outputs.containsKey(name) )
            throw new ScriptRuntimeException("Workflow output '${name}' was assigned more than once -- an output record is published under the names of its fields, so two pipelines that share an output name cannot both be published")
        owner.session.outputs[name] =
            source instanceof ChannelImpl ? source.getSource() :
            source instanceof ValueImpl ? source.getSource() :
            source instanceof DataflowWriteChannel ? source :
            CH.value(source)
    }

}
