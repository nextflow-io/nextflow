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
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.dataflow.ChannelImpl
import nextflow.dataflow.ValueImpl
import nextflow.exception.ScriptRuntimeException
import nextflow.script.dsl.Types
import nextflow.script.types.Channel
import nextflow.script.types.Value
import nextflow.util.TypeHelper
import nextflow.util.RecordMap

/**
 * Models a pipeline -- the {@code params} / {@code workflow} / {@code output}
 * trio of a script -- included into another script as a named workflow:
 *
 *   include { workflow as RNASEQ } from './pipelines/rnaseq.nf'
 *
 * The {@code params} block acts as the {@code take:} section, so the pipeline
 * is called with a single record, and the {@code output} block acts as the
 * {@code emit:} section, so the outputs published by the pipeline are emitted
 * to the calling workflow instead of being published to the output directory.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PipelineDef extends BindableDef {

    private BaseScript script

    private String name

    PipelineDef(BaseScript script, String name) {
        this.script = script
        this.name = name
    }

    BaseScript getOwner() { script }

    String getName() { name }

    String getType() { 'workflow' }

    PipelineDef cloneWithName(String name) {
        final result = (PipelineDef)clone()
        result.@name = name
        return result
    }

    @Override
    Object run(Object[] args) {
        final params = new ScriptBinding.ParamsMap(resolveParams(args))
        final published = runEntryWorkflow(params)
        return collectOutputs(script.getOutputNames(), published)
    }

    /**
     * Resolve the record argument of a pipeline call against the
     * declared params of the included pipeline.
     *
     * @param args
     */
    protected Map<String,Object> resolveParams(Object[] args) {
        final given = recordArgument(args)
        final declarations = script.getParamDeclarations()

        for( final name : given.keySet() ) {
            if( !declarations.containsKey(name) )
                throw new ScriptRuntimeException("Pipeline `${this.name}` does not declare a parameter named `${name}`")
        }

        final params = new LinkedHashMap<String,Object>(declarations.size())
        for( final decl : declarations.values() ) {
            final name = decl.name
            final value = given.containsKey(name)
                ? resolveArgument(decl, given.get(name))
                : ParamsHelper.resolveDefault(decl)

            if( value == null && !decl.optional )
                throw new ScriptRuntimeException("Parameter `${name}` of pipeline `${this.name}` is required but no value was provided")

            params.put(name, value)
        }
        return params
    }

    /**
     * The record argument of a pipeline call (e.g. an included params block).
     *
     * @param args
     */
    private Map<String,Object> recordArgument(Object[] args) {
        if( args.length == 0 )
            return Collections.<String,Object>emptyMap()
        if( args.length == 1 && args[0] instanceof RecordMap )
            return (RecordMap)args[0]
        throw new ScriptRuntimeException("Pipeline `${name}` should be called with a record, e.g. `${name}( record(input: params.input) )`")
    }

    private Object resolveArgument(Param decl, Object value) {
        if( value == null )
            return null
        // a channel is passed through as-is, so that the pipeline
        // participates in the calling workflow's dataflow graph
        if( isDataflow(value) ) {
            if( !isDataflowType(decl) )
                throw new ScriptRuntimeException("Parameter `${decl.name}` of pipeline `${this.name}` with type ${Types.getName(decl.type)} cannot be assigned to a dataflow value -- declare the param as a Channel or Value to accept it")
            checkDataflowType(decl, DataflowTypeHelper.normalizeV2(value))
            return DataflowTypeHelper.normalize(value, script.isTypingEnabled())
        }
        final type = TypeHelper.getRawType(decl.type)
        if( type == Channel || type == Value )
            throw new ScriptRuntimeException("Parameter `${decl.name}` of pipeline `${this.name}` with type ${Types.getName(decl.type)} cannot be assigned to ${value} [${Types.getName(value.getClass())}]")

        return ParamsHelper.resolveParam(decl, value, false)
    }

    private static boolean isDataflowType(Param decl) {
        final type = TypeHelper.getRawType(decl.type)
        return type == Object || Channel.isAssignableFrom(type) || Value.isAssignableFrom(type)
    }

    private void checkDataflowType(Param decl, Object value) {
        final type = TypeHelper.getRawType(decl.type)
        if( type == Channel && value !instanceof ChannelImpl || type == Value && value !instanceof ValueImpl ) {
            final actual = value instanceof ChannelImpl ? 'a Channel' : value instanceof ValueImpl ? 'a Value' : 'multiple channels'
            throw new ScriptRuntimeException("Parameter `${decl.name}` of pipeline `${this.name}` with type ${Types.getName(decl.type)} cannot be assigned to ${actual}")
        }
    }

    private static boolean isDataflow(Object value) {
        return value instanceof ChannelImpl
            || value instanceof ValueImpl
            || value instanceof DataflowWriteChannel
            || value instanceof ChannelOut
    }

    /**
     * Execute the entry workflow of the included pipeline and return
     * the outputs that it published.
     *
     * The published outputs are returned by the entry workflow rather
     * than being published, because an included pipeline emits its outputs
     * to the calling workflow -- only the calling pipeline decides what
     * is published.
     */
    protected Map<String,DataflowWriteChannel> runEntryWorkflow(ScriptBinding.ParamsMap params) {
        // the entry workflow is invoked with the pipeline name, so that
        // processes are scoped by it (e.g. `RNASEQ:STAR_ALIGN`), and with
        // the resolved params, so that `params` refers to this call
        final workflow = script.getEntryFlow().cloneWithName(name).withParams(params)
        workflow.run(BaseScriptConsts.EMPTY_ARGS)
        return workflow.getOutput().asMap()
    }

    /**
     * Collect the published outputs of the included pipeline into a record,
     * following the output block of the pipeline.
     *
     * @param declarations
     * @param published
     */
    protected Object collectOutputs(Collection<String> declarations, Map<String,DataflowWriteChannel> published) {
        for( final name : declarations ) {
            if( !published.containsKey(name) )
                throw new ScriptRuntimeException("Output '${name}' of pipeline `${this.name}` was declared in the output block but not assigned in the workflow")
        }
        for( final name : published.keySet() ) {
            if( name !in declarations )
                throw new ScriptRuntimeException("Output '${name}' of pipeline `${this.name}` was assigned in the workflow but not declared in the output block")
        }

        // the outputs are normalized for the calling script, which may be
        // typed or legacy independently of the included pipeline
        final typingEnabled = ExecutionStack.owner().isTypingEnabled()
        final result = new LinkedHashMap<String,Object>(declarations.size())
        for( final name : declarations )
            result.put(name, DataflowTypeHelper.normalize(published.get(name), typingEnabled))
        return new RecordMap(result)
    }

}
