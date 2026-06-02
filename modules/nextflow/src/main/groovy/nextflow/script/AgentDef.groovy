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

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.agent.AgentRunner
import nextflow.agent.AgentRunnerProvider
import nextflow.agent.AgentRunnerRequest
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.MapOp
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput

/**
 * Runtime model for an agent definition. Holds the captured directives,
 * inputs, outputs and prompt. Execution is delegated to the agent runner
 * (nf-agent plugin) and is not yet implemented — {@link #run} throws.
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

        final inputName = inputs[0].name
        final outputName = outputs ? outputs[0].name : 'out'
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
            final req = new AgentRunnerRequest(agentModel, agentInstruction, promptText, agentMaxIter, agentTools)
            return runner.run(req)
        }

        final out = new MapOp(source, mapper).apply()
        final channels = new LinkedHashMap<String,Object>()
        channels.put(outputName, out)
        return new ChannelOut(channels)
    }
}
