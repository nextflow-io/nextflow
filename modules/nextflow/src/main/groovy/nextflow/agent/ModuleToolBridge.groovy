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

import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Nextflow
import nextflow.extension.CH
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpec.ModuleParam
import nextflow.script.ChannelOut
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessOutput

/**
 * Bridges Nextflow processes/modules to the agent as LLM tools that execute as
 * real dataflow nodes.
 *
 * The bridge is built in the workflow body (inside {@link nextflow.script.AgentDef#run},
 * <b>before</b> the dataflow network is ignited). For every tool it derives a portable
 * {@link ToolDescriptor}, creates one persistent {@link CH#queue} per input channel and
 * pre-wires the cloned {@link ProcessDef} over those queues, capturing its output
 * {@link ChannelOut}.
 *
 * Two marshalling modes are supported:
 * <ul>
 *   <li><b>scalar</b> (no {@link ModuleSpec}): each typed scalar input param maps to one
 *       input queue; the single scalar output is serialized under its name. Schema from
 *       {@link ProcessToolSchema} (Phase 2 / 3.1).</li>
 *   <li><b>spec-driven</b> (a {@link ModuleSpec}, e.g. from a sibling {@code meta.yml}): one
 *       input queue per spec input channel. The LLM passes FLATTENED args; at dispatch the
 *       components of a tuple input channel are reassembled IN ORDER into a {@code List}
 *       (file/path → {@link Nextflow#file}, map → the arg Map, scalar → the value) and bound
 *       to that channel. Outputs (tuple or scalar) are serialized back to a JSON object keyed
 *       by emit name; file/path values become absolute path strings (the opaque-path
 *       contract). Schema from {@link ModuleSpecToolSchema} (Phase 3.2).</li>
 * </ul>
 *
 * At tool-call time (post-ignition, on the agent operator thread) {@link #call} binds the
 * marshalled args onto the input queues and blocks on the output channels' {@code .val} —
 * the task runs on the real executor in a real work dir. Calls are serialized
 * (synchronized) so inputs/outputs stay correlated. When the agent's input source is
 * exhausted {@link #close} poisons every input queue so the operators terminate.
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
        // scalar mode: the ordered input param names (one per queue)
        List<String> inputParamNames
        // scalar mode: the single output read channel + its key
        DataflowReadChannel toolOut
        String outputParamName
        // spec mode: the module spec + the captured ChannelOut (named outputs)
        ModuleSpec spec
        ChannelOut channelOut
    }

    private final Map<String,Tool> tools = new LinkedHashMap<String,Tool>()

    private final List<ToolDescriptor> descriptors = new ArrayList<ToolDescriptor>()

    private boolean closed

    /**
     * Build the bridge by pre-wiring each in-scope process tool into the dataflow
     * network (scalar mode). Must be invoked in the workflow body before ignition.
     *
     * @param resolved a name -> ProcessDef map of the in-scope tools the agent declared
     */
    ModuleToolBridge(Map<String,ProcessDef> resolved) {
        this(resolved, Collections.<String,ModuleSpec>emptyMap())
    }

    /**
     * Build the bridge, using spec-driven marshalling for any tool that has a
     * {@link ModuleSpec} (e.g. a sibling {@code meta.yml}) and the scalar typed-I/O
     * marshalling otherwise.
     *
     * @param resolved a name -> ProcessDef map of the resolved tools the agent declared
     * @param specs    a name -> ModuleSpec map; tools present here use spec-driven marshalling
     */
    ModuleToolBridge(Map<String,ProcessDef> resolved, Map<String,ModuleSpec> specs) {
        for( final entry : resolved.entrySet() ) {
            final spec = specs?.get(entry.key)
            if( spec != null )
                wireSpec(entry.key, entry.value, spec)
            else
                wireScalar(entry.key, entry.value)
        }
    }

    // -------------------------------------------------------------------------
    // SCALAR mode (Phase 2 / 3.1) - typed process, one queue per scalar input
    // -------------------------------------------------------------------------

    private void wireScalar(String name, ProcessDef proc) {
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
            name,
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

    // -------------------------------------------------------------------------
    // SPEC-driven mode (Phase 3.2) - meta.yml describes the tuple/path/map I/O
    // -------------------------------------------------------------------------

    private void wireSpec(String name, ProcessDef proc, ModuleSpec spec) {
        final inputs = spec.inputs ?: Collections.<ModuleParam>emptyList()

        // -- portable schema descriptor (flattened inputs + prose output shape)
        final inputSchema = ModuleSpecToolSchema.inputSchema(spec)
        final description = buildDescription(spec)
        final descriptor = new ToolDescriptor(name, description, inputSchema, null)

        // -- pre-wire: one queue per spec input channel
        final toolIns = new ArrayList<DataflowWriteChannel>(inputs.size())
        for( int i=0; i<inputs.size(); i++ )
            toolIns.add(CH.queue())
        final cloned = proc.cloneWithName(name)
        final ChannelOut out = (ChannelOut) cloned.run(toolIns as Object[])

        final tool = new Tool(
            name: name,
            toolIns: toolIns,
            spec: spec,
            channelOut: out )
        tools.put(name, tool)
        descriptors.add(descriptor)
    }

    private static String buildDescription(ModuleSpec spec) {
        final base = spec.description ?: spec.name ?: 'module tool'
        return "${base}\n\n${ModuleSpecToolSchema.outputDescription(spec)}".toString()
    }

    /**
     * The descriptors of the wired tools, for the agent runner request.
     */
    List<ToolDescriptor> descriptors() {
        return descriptors
    }

    /**
     * Execute a tool call. Parse the JSON args, marshal each input value onto the matching
     * pre-wired queue, then block on the output channel(s). Calls are serialized so that
     * inputs and outputs stay correlated.
     *
     * <p><b>Dispatch-level errors are returned as a tool result, not thrown.</b> An unknown
     * tool name, unparseable {@code argsJson}, or a malformed/mis-shaped argument all yield a
     * well-formed {@code {"error": "<message>"}} JSON object that names the tool and what went
     * wrong, so the LLM can see the failure and retry rather than the agent loop being aborted
     * by an exception escaping the dispatcher.
     *
     * <p><b>Limitation (documented, future hardening):</b> a failure of the underlying tool
     * <i>process task</i> (the Nextflow task the tool runs as) is NOT recovered here. Such a
     * failure surfaces through the process's own dataflow operator and aborts the session via
     * the standard error model; intercepting it as a recoverable tool result fights the runtime
     * and is deferred to a future phase. Only DISPATCH-level errors (unknown tool / argument
     * parsing / argument marshalling) are turned into a tool-result error here.
     *
     * @param toolName the name of the tool to invoke
     * @param argsJson the LLM-supplied arguments as a JSON object string
     * @return the tool result serialized as a JSON object string; on any dispatch-level failure
     *         a {@code {"error": "..."}} JSON object so the LLM can recover
     */
    @Override
    synchronized String call(String toolName, String argsJson) {
        try {
            final tool = tools.get(toolName)
            if( tool == null )
                throw new IllegalArgumentException("Unknown agent tool `${toolName}` - available tools: ${tools.keySet()}")

            final parsed = parseArgs(toolName, argsJson)

            return tool.spec != null
                ? callSpec(tool, parsed)
                : callScalar(tool, parsed)
        }
        catch( Exception e ) {
            // dispatch-level failure (unknown tool / arg parsing / arg marshalling): return it
            // as a tool result so the LLM can recover, rather than letting it abort the loop
            final message = "Agent tool `${toolName}` failed - ${e.message ?: e.toString()}".toString()
            log.warn(message, e)
            return JsonOutput.toJson([error: message])
        }
    }

    /**
     * Parse the LLM-supplied arguments into a {@code Map}. An empty/blank payload is treated as
     * no arguments; an unparseable payload, or one that is not a JSON object, raises an error
     * that names the offending tool so it round-trips to the LLM as a clear tool-result error.
     */
    private static Map parseArgs(String toolName, String argsJson) {
        if( argsJson == null || argsJson.trim().isEmpty() )
            return Collections.emptyMap()
        final Object parsed
        try {
            parsed = new JsonSlurper().parseText(argsJson)
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("could not parse the tool arguments as JSON (${e.message}); arguments must be a JSON object")
        }
        if( !(parsed instanceof Map) )
            throw new IllegalArgumentException("the tool arguments must be a JSON object, got ${parsed?.getClass()?.simpleName ?: 'null'}")
        return (Map) parsed
    }

    private String callScalar(Tool tool, Map parsed) {
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

    private String callSpec(Tool tool, Map parsed) {
        final spec = tool.spec
        final inputs = spec.inputs ?: Collections.<ModuleParam>emptyList()

        // -- marshal the flattened args into one channel value per input channel
        for( int i=0; i<inputs.size(); i++ ) {
            final param = inputs[i]
            final value = marshalInput(param, parsed)
            tool.toolIns[i].bind(value)
        }

        // -- block on each output channel, in declared order, serializing the records
        final result = new LinkedHashMap<String,Object>()
        final outputs = spec.outputs ?: Collections.<ModuleParam>emptyList()
        final emitNames = emitNamesByPosition(tool.channelOut)
        for( int i=0; i<outputs.size(); i++ ) {
            final param = outputs[i]
            final key = outputKey(param, emitNames, i)
            final ch = CH.getReadChannel(tool.channelOut[i])
            final value = ch.val
            result.put(key, serializeOutput(param, value))
        }
        return JsonOutput.toJson(result)
    }

    /**
     * Reassemble the channel value the process expects for the given input channel from the
     * flattened LLM args: a tuple channel becomes a {@code List} of its components in order;
     * a scalar channel becomes the single marshalled value.
     */
    private static Object marshalInput(ModuleParam param, Map parsed) {
        if( param.isTuple() ) {
            final list = new ArrayList<Object>(param.components.size())
            for( final comp : param.components )
                list.add(marshalComponent(comp, parsed.get(comp.name)))
            return list
        }
        return marshalComponent(param, parsed.get(param.name))
    }

    /**
     * Marshal a single flattened arg into the runtime value: file/path/directory args (path
     * strings) become a {@link Path} via {@link Nextflow#file}; everything else passes through.
     */
    private static Object marshalComponent(ModuleParam comp, Object arg) {
        if( arg == null )
            return null
        if( ModuleSpecToolSchema.isFileType(comp.type?.toLowerCase()) )
            return Nextflow.file(arg.toString())
        return arg
    }

    /**
     * Serialize an output channel value per the spec param shape: a tuple value (a List)
     * becomes an object keyed by component name; a scalar value is serialized directly.
     * File/path values become absolute path strings.
     */
    private static Object serializeOutput(ModuleParam param, Object value) {
        if( param.isTuple() ) {
            final record = new LinkedHashMap<String,Object>()
            final list = (value instanceof List) ? (List) value : [value]
            final comps = param.components
            for( int i=0; i<comps.size(); i++ ) {
                final comp = comps[i]
                final v = (i < list.size()) ? list[i] : null
                record.put(comp.name ?: "value${i}".toString(), serializeValue(comp, v))
            }
            return record
        }
        return serializeValue(param, value)
    }

    private static Object serializeValue(ModuleParam comp, Object value) {
        if( value == null )
            return null
        if( ModuleSpecToolSchema.isFileType(comp.type?.toLowerCase()) )
            return pathString(value)
        return value
    }

    private static String pathString(Object value) {
        if( value instanceof Path )
            return ((Path) value).toAbsolutePath().toString()
        return value.toString()
    }

    /**
     * Map each positional output channel of the {@link ChannelOut} to its emit name (when
     * the channel was declared with {@code emit:}). Used as the JSON key for tuple outputs
     * whose spec param carries no channel-level name.
     */
    private static Map<Integer,String> emitNamesByPosition(ChannelOut out) {
        final result = new LinkedHashMap<Integer,String>()
        for( final emitName : out.getNames() ) {
            final ch = out.getProperty(emitName)
            for( int i=0; i<out.size(); i++ ) {
                if( out[i].is(ch) ) {
                    result.put(i, emitName)
                    break
                }
            }
        }
        return result
    }

    private static String outputKey(ModuleParam param, Map<Integer,String> emitNames, int index) {
        if( param.name )
            return param.name
        final emit = emitNames.get(index)
        if( emit )
            return emit
        return "out${index}".toString()
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

    /**
     * The output property key: a named output uses its declared name; a bare typed value
     * lowers to the synthetic name {@code $out}.
     */
    private static String outputKey(ProcessOutput param) {
        final name = param.getName()
        return ( name == null ) ? '$out' : name
    }
}
