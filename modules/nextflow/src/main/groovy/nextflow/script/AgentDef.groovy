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

import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.agent.AgentRunner
import nextflow.agent.AgentRunnerProvider
import nextflow.agent.AgentRunnerRequest
import nextflow.agent.RecordSchema
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.MapOp
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.util.TypeHelper

/**
 * Runtime model for an agent definition. Holds the captured directives,
 * inputs, outputs and prompt. {@link #run} executes the agent as a dataflow
 * operator, rendering the prompt per input record and delegating the LLM work
 * to an {@link nextflow.agent.AgentRunner} resolved from the active plugins
 * (e.g. nf-agent).
 */
@Slf4j
@CompileStatic
class AgentDef extends BindableDef implements ChainableDef {

    static final String TYPE = 'agent'

    private BaseScript owner
    private String name
    private String simpleName
    private Map<String,Object> directives
    private List<AgentInput> inputs
    private List<AgentOutput> outputs
    private PromptDef prompt

    AgentDef(BaseScript owner, String name, Map<String,Object> directives, List<AgentInput> inputs, List<AgentOutput> outputs, PromptDef prompt) {
        this.owner = owner
        this.name = name
        this.simpleName = name
        this.directives = directives
        this.inputs = inputs
        this.outputs = outputs
        this.prompt = prompt
    }

    @Override String getType() { TYPE }
    @Override String getName() { name }
    String getSimpleName() { simpleName }
    BaseScript getOwner() { owner }

    String getModel() { directives.get('model') as String }
    String getInstruction() { directives.get('instruction') as String }
    List getTools() { (directives.get('tools') ?: []) as List }
    Integer getMaxIterations() { directives.get('maxIterations') as Integer }
    List<AgentInput> getInputs() { inputs }
    List<AgentOutput> getOutputs() { outputs }
    PromptDef getPrompt() { prompt }

    @Override
    ComponentDef cloneWithName(String name) {
        def copy = (AgentDef) this.clone()
        copy.@name = name
        copy.@simpleName = name.contains(':') ? name.tokenize(':').last() : name
        return copy
    }

    private DataflowReadChannel createSourceChannel(Object value) {
        if( value instanceof DataflowReadChannel || value instanceof DataflowBroadcast )
            return CH.getReadChannel(value)
        final result = CH.value()
        result.bind(value)
        return result
    }

    @Override
    @CompileDynamic
    Object run(Object[] args0) {
        final args = ChannelOut.spread(args0)
        if( inputs.size() != 1 )
            throw new ScriptRuntimeException("Agent `${name}` must declare exactly one input (got ${inputs.size()}) - multiple/zero inputs are not yet supported")
        if( args.size() != 1 )
            throw new ScriptRuntimeException("Agent `${name}` expects 1 input channel but received ${args.size()}")

        if( !outputs )
            throw new ScriptRuntimeException("Agent `${name}` must declare exactly one output - zero outputs are not yet supported")

        final inputName = inputs[0].name
        final outputName = outputs[0].name
        final outputClass = outputs[0].type as Class
        final outputSchema = RecordSchema.of(outputClass)
        final source = createSourceChannel(args[0])
        final AgentRunner runner = AgentRunnerProvider.get()
        final agentModel = this.model
        final agentInstruction = this.instruction
        final agentTools = this.tools
        final agentMaxIter = (this.maxIterations != null ? this.maxIterations : 20) as int
        final promptDef = this.prompt

        final mapper = { Object item ->
            final cl = (Closure) promptDef.closure.clone()
            cl.setDelegate([(inputName): item])
            cl.setResolveStrategy(Closure.DELEGATE_FIRST)
            final promptText = cl.call()?.toString()
            final inputJson = toJson(item)
            final req = new AgentRunnerRequest(agentModel, agentInstruction, promptText, agentMaxIter, agentTools, outputSchema, inputJson)
            final json = runner.run(req)
            final map = new JsonSlurper().parseText(stripFences(json)) as Map
            return TypeHelper.asRecordType(map, outputClass)
        }

        final out = new MapOp(source, mapper).apply()
        final channels = new LinkedHashMap<String,Object>()
        channels.put(outputName, out)
        return new ChannelOut(channels)
    }

    /**
     * Serialize the input record (a Map at runtime) to JSON, rendering any
     * {@link Path} value as its absolute string so the model receives a portable
     * representation.
     */
    protected static String toJson(Object item) {
        return JsonOutput.toJson(normalizeForJson(item))
    }

    private static Object normalizeForJson(Object value) {
        if( value instanceof Path )
            return value.toAbsolutePath().toString()
        if( value instanceof Map )
            return value.collectEntries { k, v -> [(k): normalizeForJson(v)] }
        if( value instanceof Collection )
            return value.collect { normalizeForJson(it) }
        return value
    }

    /**
     * Strip a leading ```json (or ```) fence and a trailing ``` fence from the
     * given text, returning the inner content. If no fences are present the
     * original text is returned unchanged.
     */
    protected static String stripFences(String text) {
        if( text == null )
            return null
        def s = text.trim()
        if( !s.startsWith('```') )
            return text
        // drop the opening fence line (``` or ```json)
        final nl = s.indexOf('\n')
        if( nl < 0 )
            return text
        s = s.substring(nl + 1)
        // drop the trailing fence
        final close = s.lastIndexOf('```')
        if( close >= 0 )
            s = s.substring(0, close)
        return s.trim()
    }
}
