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
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.agent.AgentRunner
import nextflow.agent.AgentRunnerProvider
import nextflow.agent.AgentRunnerRequest
import nextflow.agent.ModuleToolBridge
import nextflow.agent.RecordSchema
import nextflow.agent.ToolDescriptor
import nextflow.agent.ToolDispatcher
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.types.Record
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
    List getTools() {
        final value = directives.get('tools')
        if( value == null )
            return []
        return value instanceof List ? (List) value : [value]
    }
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
        // structured output is opt-in: derive a JSON schema and bind the runner's
        // JSON result to a record only when the output type is a record type;
        // otherwise the runner's text is emitted verbatim
        final boolean structured = outputClass != null && Record.isAssignableFrom(outputClass)
        final outputSchema = structured ? RecordSchema.of(outputClass) : null
        final source = createSourceChannel(args[0])
        final AgentRunner runner = AgentRunnerProvider.get()
        final agentModel = this.model
        final agentInstruction = this.instruction
        final agentTools = this.tools
        final agentMaxIter = (this.maxIterations != null ? this.maxIterations : 20) as int
        final promptDef = this.prompt

        // resolve declared `tools` to in-scope processes and pre-wire them into the
        // dataflow network (before ignition); the bridge is poisoned on completion
        final ModuleToolBridge bridge = createToolBridge()
        final List<ToolDescriptor> toolSpecs = bridge != null ? bridge.descriptors() : null

        final mapper = { Object item ->
            final cl = (Closure) promptDef.closure.clone()
            cl.setDelegate([(inputName): item])
            cl.setResolveStrategy(Closure.DELEGATE_FIRST)
            final promptText = cl.call()?.toString()
            final inputJson = toJson(item)
            final req = new AgentRunnerRequest(agentModel, agentInstruction, promptText, agentMaxIter, agentTools, outputSchema, inputJson, toolSpecs, (bridge as ToolDispatcher))
            final result = runner.run(req)
            if( !structured )
                return result
            final map = new JsonSlurper().parseText(stripFences(result)) as Map
            return TypeHelper.asRecordType(map, outputClass)
        }

        final out = applyAgentOperator(source, mapper, bridge)
        final channels = new LinkedHashMap<String,Object>()
        channels.put(outputName, out)
        return new ChannelOut(channels)
    }

    /**
     * Resolve the agent's declared {@code tools} entries (Phase 2: each is a String naming
     * an in-scope process) to their {@link ProcessDef}s and pre-wire them as agent tools.
     *
     * @return a {@link ModuleToolBridge} wiring the resolved processes, or {@code null} when
     *         no tools are declared
     */
    @CompileDynamic
    private ModuleToolBridge createToolBridge() {
        final declared = this.tools
        if( !declared )
            return null
        final meta = ScriptMeta.get(owner)
        final resolved = new LinkedHashMap<String, ProcessDef>()
        for( final entry : declared ) {
            final toolName = entry?.toString()
            final proc = meta?.getProcess(toolName)
            if( proc == null )
                throw new ScriptRuntimeException("Agent `${name}` tool `${toolName}` is not a process in scope")
            resolved.put(toolName, proc)
        }
        return new ModuleToolBridge(resolved)
    }

    /**
     * Build the agent operator over the input {@code source}, applying the {@code mapper}
     * per item. Mirrors {@link nextflow.extension.MapOp} but, when a {@code bridge} is
     * present, attaches a listener whose {@code afterStop} poisons every pre-wired tool
     * input so the session can terminate once the agent's input is exhausted.
     */
    @CompileDynamic
    private DataflowWriteChannel applyAgentOperator(DataflowReadChannel source, Closure mapper, ModuleToolBridge bridge) {
        final DataflowWriteChannel target = CH.createBy(source)
        final boolean stopOnFirst = source instanceof DataflowExpression

        final listener = new DataflowEventAdapter() {
            @Override
            void afterStop(final DataflowProcessor processor) {
                bridge?.close()
            }
            @Override
            boolean onException(final DataflowProcessor processor, final Throwable t) {
                // poison tools so the network can unwind even on failure
                bridge?.close()
                log.error("@unknown", t)
                Global.session?.abort(t)
                return true
            }
        }

        DataflowHelper.newOperator(source, target, listener) { it ->
            final result = mapper.call(it)
            final proc = (DataflowProcessor) getDelegate()
            if( result != Channel.VOID )
                proc.bindOutput(result)
            if( result == Channel.STOP || stopOnFirst )
                proc.terminate()
        }

        return target
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
