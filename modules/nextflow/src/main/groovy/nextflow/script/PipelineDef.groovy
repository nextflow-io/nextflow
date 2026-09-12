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
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.dataflow.ChannelImpl
import nextflow.dataflow.ValueImpl
import nextflow.exception.ScriptRuntimeException
import nextflow.util.RecordMap

/**
 * Models a pipeline -- the {@code params} / {@code workflow} / {@code output}
 * trio of a script -- included into another script as a named workflow:
 *
 *   include { workflow as RNASEQ } from './pipelines/rnaseq.nf'
 *
 * The {@code params} block acts as the {@code take:} section, so the pipeline
 * is called with named arguments, and the {@code output} block acts as the
 * {@code emit:} section, so the outputs published by the pipeline are emitted
 * to the calling workflow instead of being published to the output directory.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class PipelineDef extends BindableDef implements ChainableDef {

    private BaseScript pipeline

    private String name

    private OutputTypeDef outputType

    PipelineDef(BaseScript pipeline, String name='workflow') {
        this.pipeline = pipeline
        this.name = name
    }

    BaseScript getOwner() { pipeline }

    /**
     * Declare the output type included alongside this pipeline, i.e. the
     * `output` name in the same include statement.
     *
     * The output directives of a pipeline can refer to its params, so they
     * are resolved when the pipeline is executed and handed to the output
     * type, rather than being resolved by the calling pipeline.
     *
     * @param outputType
     */
    void setOutputType(OutputTypeDef outputType) {
        this.outputType = outputType
    }

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
        final declarations = pipeline.getOutputDeclarations(params)
        outputType?.setDeclarations(declarations)
        return collectOutputs(declarations.keySet(), published)
    }

    /**
     * Resolve the named arguments of a pipeline call against the
     * declared params of the included pipeline.
     *
     * @param args
     */
    protected Map<String,Object> resolveParams(Object[] args) {
        final given = namedArguments(args)
        final declarations = pipeline.getParamDeclarations()

        for( final name : given.keySet() ) {
            if( !declarations.containsKey(name) )
                throw new ScriptRuntimeException("Pipeline `${this.name}` does not declare a parameter named `${name}`")
        }

        final params = new LinkedHashMap<String,Object>(declarations.size())
        for( final decl : declarations.values() ) {
            final name = decl.name
            // a null argument is treated as absent, so that the default
            // value of the corresponding param is applied
            final value =
                given.get(name) != null ? resolveArgument(decl, given.get(name)) :
                decl.defaultValue != null ? ParamsHelper.resolveFromCode(decl, decl.defaultValue) :
                ParamsHelper.emptyRecord(decl)

            if( value == null && !decl.optional )
                throw new ScriptRuntimeException("Parameter `${name}` of pipeline `${this.name}` is required but no value was provided")

            params.put(name, value)
        }
        return params
    }

    /**
     * The named arguments of a pipeline call, either as named arguments or
     * as a single record (e.g. an included params block).
     *
     * @param args
     */
    private Map<String,Object> namedArguments(Object[] args) {
        if( args.length == 0 )
            return Collections.<String,Object>emptyMap()
        if( args.length == 1 && args[0] instanceof Map )
            return (Map<String,Object>)args[0]
        throw new ScriptRuntimeException("Pipeline `${name}` should be called with named arguments, e.g. `${name}( input: params.input )`")
    }

    private Object resolveArgument(Param decl, Object value) {
        if( value == null )
            return null
        // a channel is passed through as-is, so that the pipeline
        // participates in the calling workflow's dataflow graph
        if( isDataflow(value) )
            return DataflowTypeHelper.normalize(value, pipeline.isTypingEnabled())

        final result = ParamsHelper.resolveFromCode(decl, value)
        ParamsHelper.checkAssignable(decl, result)
        return result
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
     * The published outputs are captured from the session rather than
     * being published, because an included pipeline emits its outputs
     * to the calling workflow -- only the calling pipeline decides what
     * is published.
     */
    protected Map<String,DataflowWriteChannel> runEntryWorkflow(ScriptBinding.ParamsMap params) {
        final outputs = pipeline.session.outputs
        final saved = new LinkedHashMap<String,DataflowWriteChannel>(outputs)
        outputs.clear()
        try {
            // the entry workflow is invoked with the pipeline name, so that
            // processes are scoped by it (e.g. `RNASEQ:STAR_ALIGN`), and with
            // the resolved params, so that `params` refers to this call
            pipeline.getEntryFlow().cloneWithName(name).withParams(params).run(BaseScriptConsts.EMPTY_ARGS)
            return new LinkedHashMap<String,DataflowWriteChannel>(outputs)
        }
        finally {
            outputs.clear()
            outputs.putAll(saved)
        }
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
