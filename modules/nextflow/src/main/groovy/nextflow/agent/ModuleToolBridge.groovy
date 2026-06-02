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
package nextflow.agent

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.extension.CH
import nextflow.script.ChannelOut
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessOutput

/**
 * Bridges in-scope Nextflow processes to the agent as LLM tools that execute as
 * real dataflow nodes.
 *
 * The bridge is built in the workflow body (inside {@link nextflow.script.AgentDef#run},
 * <b>before</b> the dataflow network is ignited). For every tool it:
 * <ul>
 *   <li>derives a portable {@link ToolDescriptor} (input/output JSON schema) via {@link ProcessToolSchema};</li>
 *   <li>creates one persistent {@link CH#queue} per process input param and pre-wires the
 *       cloned {@link ProcessDef} over those queues, capturing the single output read channel.</li>
 * </ul>
 *
 * At tool-call time (post-ignition, on the agent operator thread) {@link #call} binds the
 * marshalled args onto the input queues and blocks on the output channel's {@code .val} —
 * the task runs on the real executor in a real work dir. Calls are serialized
 * (synchronized) so that one {@code bind}/{@code .val} pair completes before the next,
 * keeping inputs and outputs correctly correlated.
 *
 * When the agent's input source is exhausted {@link #close} poisons every input queue so
 * the pre-wired operators terminate and {@code session.await()} can return.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleToolBridge implements ToolDispatcher {

    /**
     * Holds the pre-wired dataflow plumbing for a single tool.
     */
    private static class Tool {
        String name
        List<DataflowWriteChannel> toolIns
        List<String> inputParamNames
        DataflowReadChannel toolOut
        String outputParamName
    }

    private final Map<String,Tool> tools = new LinkedHashMap<String,Tool>()

    private final List<ToolDescriptor> descriptors = new ArrayList<ToolDescriptor>()

    private boolean closed

    /**
     * Build the bridge by pre-wiring each in-scope process tool into the dataflow
     * network. Must be invoked in the workflow body before the network is ignited.
     *
     * @param resolved a name -> ProcessDef map of the in-scope tools the agent declared
     */
    ModuleToolBridge(Map<String,ProcessDef> resolved) {
        for( final entry : resolved.entrySet() ) {
            wire(entry.key, entry.value)
        }
    }

    private void wire(String name, ProcessDef proc) {
        final config = proc.getProcessConfig()
        if( !(config instanceof ProcessConfigV2) )
            throw new IllegalArgumentException("Agent tool `${name}` must be a typed process to be used as an agent tool")
        final cfg = (ProcessConfigV2) config

        // -- ordered input param names
        final inputParams = cfg.getInputs().getParams()
        final inputParamNames = new ArrayList<String>(inputParams.size())
        for( final p : inputParams )
            inputParamNames.add(((ProcessInput) p).getName())

        // -- single output param (Phase 2 assumption)
        final outputParams = cfg.getOutputs().getParams()
        if( outputParams.size() != 1 )
            throw new IllegalArgumentException("Agent tool `${name}` must declare exactly one output to be used as an agent tool (got ${outputParams.size()})")
        final outputParamName = outputKey((ProcessOutput) outputParams[0])

        // -- portable schema descriptor
        final descriptor = new ToolDescriptor(
            name,
            descriptionOf(proc, name),
            ProcessToolSchema.inputSchema(proc),
            ProcessToolSchema.outputSchema(proc) )

        // -- pre-wire the tool: one queue per input param, clone & invoke once
        final toolIns = new ArrayList<DataflowWriteChannel>(inputParamNames.size())
        for( int i=0; i<inputParamNames.size(); i++ )
            toolIns.add(CH.queue())
        final cloned = proc.cloneWithName(name)
        final ChannelOut out = (ChannelOut) cloned.run(toolIns as Object[])
        final DataflowReadChannel toolOut = CH.getReadChannel(out[0])

        final tool = new Tool(
            name: name,
            toolIns: toolIns,
            inputParamNames: inputParamNames,
            toolOut: toolOut,
            outputParamName: outputParamName )
        tools.put(name, tool)
        descriptors.add(descriptor)
    }

    /**
     * The descriptors of the wired tools, for the agent runner request.
     */
    List<ToolDescriptor> descriptors() {
        return descriptors
    }

    /**
     * Execute a tool call. Parse the JSON args, bind each input value onto the matching
     * pre-wired queue, then block on the output channel. Calls are serialized so that
     * inputs and outputs stay correlated.
     *
     * @param toolName the name of the tool to invoke
     * @param argsJson the LLM-supplied arguments as a JSON object string
     * @return the tool result serialized as a JSON object string; on any failure a
     *         {@code {"error": "..."}} JSON object so the LLM can recover
     */
    @Override
    synchronized String call(String toolName, String argsJson) {
        try {
            final tool = tools.get(toolName)
            if( tool == null )
                throw new IllegalArgumentException("Unknown agent tool `${toolName}`")

            final parsed = (argsJson != null && !argsJson.trim().isEmpty())
                ? new JsonSlurper().parseText(argsJson) as Map
                : Collections.emptyMap()

            // -- bind each input value (in declared order) onto its queue
            for( int i=0; i<tool.inputParamNames.size(); i++ ) {
                final paramName = tool.inputParamNames[i]
                tool.toolIns[i].bind(parsed.get(paramName))
            }

            // -- blocking pull: the task runs on the real executor
            final result = tool.toolOut.val

            final key = (tool.outputParamName == '$out') ? 'result' : tool.outputParamName
            return JsonOutput.toJson([(key): result])
        }
        catch( Exception e ) {
            log.warn("Agent tool `${toolName}` failed - ${e.message}", e)
            return JsonOutput.toJson([error: (e.message ?: e.toString())])
        }
    }

    /**
     * Poison every input queue so the pre-wired tool operators terminate and the session
     * can complete. Idempotent.
     */
    synchronized void close() {
        if( closed )
            return
        closed = true
        for( final tool : tools.values() ) {
            for( final in : tool.toolIns )
                in.bind(PoisonPill.instance)
        }
    }

    private static String descriptionOf(ProcessDef proc, String name) {
        return name
    }

    /**
     * The output property key: a named output uses its declared name; a bare typed value
     * lowers to the synthetic name {@code $out}.
     */
    private static String outputKey(ProcessOutput param) {
        final name = param.getName()
        return ( name == null ) ? '$out' : name
    }
}
