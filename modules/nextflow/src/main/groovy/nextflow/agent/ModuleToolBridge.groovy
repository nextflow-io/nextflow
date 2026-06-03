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
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import io.seqera.npr.api.schema.v1.ModuleMetadata
import nextflow.Nextflow
import nextflow.extension.CH
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpec.ModuleParam
import nextflow.script.ChannelOut
import nextflow.script.ProcessConfigV1
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.ProcessEntryHandler
import nextflow.script.params.OutParam
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
        // spec mode: the cloned ProcessDef the tool runs as, used to bind the LLM args
        // onto the input channels via the SAME logic as `module run` (ProcessEntryHandler)
        ProcessDef processDef
    }

    /**
     * The public registry {@link ModuleMetadata} for a tool plus its nf-core scope flag (for the
     * {@code meta.id} convention). When present for a spec-driven tool, the tool DESCRIPTOR
     * (description + input schema) is sourced from the metadata instead of the sibling
     * {@code meta.yml} {@link ModuleSpec}; marshalling/output stays on the spec + ProcessDef.
     */
    @TupleConstructor
    static class RegistryMeta {
        final ModuleMetadata metadata
        final boolean nfCore
    }

    private final Map<String,Tool> tools = new LinkedHashMap<String,Tool>()

    private final List<ToolDescriptor> descriptors = new ArrayList<ToolDescriptor>()

    private boolean closed

    /** The maximum size (in bytes) of a structured file output whose content is inlined to the
     * LLM; larger or non text-like/binary outputs stay opaque path handles. Defaults to 32 KB. */
    private long maxInlineBytes = ToolOutputReader.DEFAULT_INLINE_BYTES

    /** Set the cap (in bytes) for inlining structured file outputs; non-positive values reset
     * to the 32 KB default. Wired from the {@code agent.maxToolOutputInlineSize} config scope. */
    void setMaxInlineBytes(long bytes) { this.maxInlineBytes = bytes > 0 ? bytes : ToolOutputReader.DEFAULT_INLINE_BYTES }

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
        this(resolved, specs, Collections.<String,RegistryMeta>emptyMap())
    }

    /**
     * Build the bridge, sourcing the tool DESCRIPTOR from the public registry
     * {@link ModuleMetadata} when present (Task 3) and otherwise from the sibling {@code meta.yml}
     * {@link ModuleSpec}. Marshalling/output is always spec-driven (the spec + ProcessDef) for a
     * tool that has a {@link ModuleSpec}; tools with neither use the scalar typed-I/O path.
     *
     * @param resolved   a name -> ProcessDef map of the resolved tools the agent declared
     * @param specs      a name -> ModuleSpec map; tools present here use spec-driven marshalling
     * @param metadatas  a name -> {@link RegistryMeta} map; for a spec-driven tool present here the
     *                   descriptor is sourced from the registry metadata (description + schema)
     */
    ModuleToolBridge(Map<String,ProcessDef> resolved, Map<String,ModuleSpec> specs, Map<String,RegistryMeta> metadatas) {
        for( final entry : resolved.entrySet() ) {
            final spec = specs?.get(entry.key)
            final registryMeta = metadatas?.get(entry.key)
            if( spec != null )
                wireSpec(entry.key, entry.value, spec, registryMeta?.metadata, registryMeta != null ? registryMeta.nfCore : false)
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

    private void wireSpec(String name, ProcessDef proc, ModuleSpec spec, ModuleMetadata metadata = null, boolean nfCore = false) {
        final inputs = spec.inputs ?: Collections.<ModuleParam>emptyList()

        // -- portable schema descriptor (flattened inputs + prose output shape). When the public
        //    registry ModuleMetadata is available (Task 3) it is the canonical, richer source
        //    (per-field descriptions, patterns, enums, the nf-core meta.id convention, tool
        //    homepages); otherwise fall back to the sibling meta.yml ModuleSpec (offline /
        //    local-file / in-scope tools). Marshalling/output is UNCHANGED (always spec-driven).
        final Map inputSchema
        final String description
        if( metadata != null ) {
            inputSchema = ModuleMetadataToolSchema.inputSchema(metadata, nfCore)
            description = ModuleMetadataToolSchema.description(metadata)
            // consistency / silent-drift guard: the descriptor (from the registry `.latest`
            // metadata) and the executor (the installed meta.yml ModuleSpec used for marshalling)
            // are two views of the same module - warn if their flattened input names differ.
            warnOnInputDrift(name, metadata, spec)
        }
        else {
            inputSchema = ModuleSpecToolSchema.inputSchema(spec)
            description = buildDescription(spec)
        }
        final descriptor = new ToolDescriptor(name, description, inputSchema, null)

        // -- the number of input queues = the PROCESS's declared input-channel count;
        //    `ProcessEntryHandler.getProcessArguments` returns one value per declared input
        //    channel, so the queues must match that, NOT the spec input count. They should be
        //    equal (descriptor and executor are two views of the same meta.yml) -- warn on drift.
        final int nInputs = declaredInputChannelCount(proc)
        if( nInputs != inputs.size() )
            log.warn("Agent tool `${name}`: module spec declares ${inputs.size()} input(s) but the process declares ${nInputs} input channel(s) - using the process input count for marshalling")

        // -- pre-wire: one queue per process input channel
        final toolIns = new ArrayList<DataflowWriteChannel>(nInputs)
        for( int i=0; i<nInputs; i++ )
            toolIns.add(CH.queue())
        final cloned = proc.cloneWithName(name)
        final ChannelOut out = (ChannelOut) cloned.run(toolIns as Object[])

        final tool = new Tool(
            name: name,
            toolIns: toolIns,
            spec: spec,
            channelOut: out,
            processDef: cloned )
        tools.put(name, tool)
        descriptors.add(descriptor)
    }

    /**
     * The number of declared input channels of a process: each {@code val}/{@code path}/... is one
     * channel and a {@code tuple} counts as a single channel. Matches the per-input-channel arity
     * that {@link nextflow.script.ProcessEntryHandler#getProcessArguments} returns.
     */
    private static int declaredInputChannelCount(ProcessDef proc) {
        final config = proc.getProcessConfig()
        if( config instanceof ProcessConfigV2 )
            return ((ProcessConfigV2) config).getInputs().getParams().size()
        return ((ProcessConfigV1) config).getInputs().size()
    }

    /**
     * The ProcessDef's declared output params, in declaration order. These are positionally
     * aligned with the captured {@link ChannelOut} (the ProcessDef builds the {@code ChannelOut}
     * from these same outputs, in order), so {@code outParams[i]} corresponds to
     * {@code channelOut[i]}. Used to detect topic-routed outputs authoritatively.
     */
    private static List<OutParam> declaredOutputParams(ProcessDef proc) {
        final config = proc.getProcessConfig()
        if( config instanceof ProcessConfigV2 )
            return new ArrayList<OutParam>(((ProcessConfigV2) config).getOutputs().getParams())
        return new ArrayList<OutParam>(((ProcessConfigV1) config).getOutputs())
    }

    /**
     * True when the ProcessDef's declared output at position {@code i} is routed to a
     * {@code topic:} (nf-core {@code versions}-style bookkeeping). Such outputs use a
     * topic-source channel that never binds a readable per-invocation value, so the dispatcher
     * must not block reading them. This is authoritative -- the meta.yml {@code type} is
     * unreliable: some modules type the eval/version component as {@code string}, not
     * {@code eval} (e.g. nf-core/assemblyscan vs nf-core/skesa), so {@link #isEvalOutput} alone
     * would miss it and the dispatch would hang.
     */
    private static boolean isTopicOutput(List<OutParam> outParams, int i) {
        if( outParams == null || i >= outParams.size() )
            return false
        final p = outParams[i]
        return p != null && p.getChannelTopicName() != null
    }

    private static String buildDescription(ModuleSpec spec) {
        final base = spec.description ?: spec.name ?: 'module tool'
        return "${base}\n\n${ModuleSpecToolSchema.outputDescription(spec)}".toString()
    }

    /**
     * Consistency check: warn when the flattened input property names derived from the registry
     * {@link ModuleMetadata} (the descriptor source) differ from those of the executable
     * {@link ModuleSpec} (the marshalling source) - a silent-drift guard between the registry
     * {@code .latest} descriptor and the installed spec the tool actually runs as.
     */
    private static void warnOnInputDrift(String name, ModuleMetadata metadata, ModuleSpec spec) {
        final fromMeta = ModuleMetadataToolSchema.inputPropertyNames(metadata)
        final fromSpec = specInputNames(spec)
        if( fromMeta != fromSpec )
            log.warn("Agent tool `${name}`: registry metadata declares inputs ${fromMeta} but the installed module spec declares ${fromSpec} - the tool descriptor (from the registry) and its marshalling (from the local spec) may be out of sync")
    }

    /** The flattened input property names of a {@link ModuleSpec}, in declaration order. */
    private static List<String> specInputNames(ModuleSpec spec) {
        final names = new ArrayList<String>()
        final inputs = spec?.inputs ?: Collections.<ModuleParam>emptyList()
        for( final ModuleParam param : inputs ) {
            final comps = param.isTuple() ? param.components : Collections.singletonList(param)
            for( final ModuleParam comp : comps ) {
                if( comp?.name )
                    names.add(comp.name)
            }
        }
        return names
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

        // -- marshal the flattened LLM args into one channel value per input channel using the
        //    SAME logic as `nextflow module run` (ProcessEntryHandler.getProcessArguments): the
        //    args map is treated as the params map (dot-notation folded into nested maps), each
        //    declared input is looked up by name, coerced per the module spec types (file ->
        //    Nextflow.file, map -> Map, integer -> Integer, ...) and a tuple input is assembled
        //    into a List of its components (e.g. [[id:'s1'], file(reads)]).
        //    A missing/invalid arg throws IllegalArgumentException which the dispatcher in `call`
        //    turns into a {"error":...} tool result so the LLM can recover.
        final List args = ProcessEntryHandler.getProcessArguments(tool.processDef, parsed, spec)
        for( int i=0; i<tool.toolIns.size(); i++ )
            tool.toolIns[i].bind(args[i])

        // -- block on each output channel, in declared order, serializing the records
        final result = new LinkedHashMap<String,Object>()
        final outputs = spec.outputs ?: Collections.<ModuleParam>emptyList()
        final emitNames = emitNamesByPosition(tool.channelOut)
        final int nOut = tool.channelOut.size()
        // the ProcessDef's declared output params, positionally aligned with the output
        // channels (ProcessDef builds the ChannelOut from these in the same order); this is
        // the AUTHORITATIVE source for detecting topic-routed (`versions`-style) outputs
        final List<OutParam> outParams = declaredOutputParams(tool.processDef)
        for( int i=0; i<outputs.size(); i++ ) {
            final param = outputs[i]
            // Skip nf-core `versions`-style outputs: a tuple carrying an `eval`
            // component is pipeline bookkeeping routed to a `topic`, not a
            // per-invocation tool result -- and its channel never binds a readable
            // value, so reading `.val` here would block the dispatch forever.
            // Topic-routed outputs are detected authoritatively from the ProcessDef
            // (the meta.yml `type` is unreliable -- see isTopicOutput).
            if( isEvalOutput(param) || isTopicOutput(outParams, i) )
                continue
            // Defensively skip an output with no corresponding emitted data channel
            // (e.g. a topic-only output) rather than block on a non-existent channel.
            if( i >= nOut )
                continue
            final key = outputKey(param, emitNames, i)
            final ch = CH.getReadChannel(tool.channelOut[i])
            final value = ch.val
            result.put(key, serializeOutput(param, value))
        }
        return JsonOutput.toJson(result)
    }

    /**
     * True for an nf-core `versions`-style output: a tuple whose components include an
     * {@code eval} (or the param itself is an {@code eval}). Such outputs are computed
     * command captures routed to a `topic` for pipeline bookkeeping -- they are not a
     * per-invocation tool result and their channel does not bind a readable value, so
     * the dispatcher must not block reading them.
     */
    private static boolean isEvalOutput(ModuleParam param) {
        if( param == null )
            return false
        if( 'eval'.equalsIgnoreCase(param.type) )
            return true
        final comps = param.components
        if( comps != null ) {
            for( final ModuleParam c : comps )
                if( c != null && 'eval'.equalsIgnoreCase(c.type) )
                    return true
        }
        return false
    }

    /**
     * Serialize an output channel value per the spec param shape: a tuple value (a List)
     * becomes an object keyed by component name; a scalar value is serialized directly.
     * File/path values become absolute path strings, except small structured (text-like)
     * outputs whose contents are inlined for the LLM (see {@link ToolOutputReader}).
     */
    private Object serializeOutput(ModuleParam param, Object value) {
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

    private Object serializeValue(ModuleParam comp, Object value) {
        if( value == null )
            return null
        if( ModuleSpecToolSchema.isFileType(comp.type?.toLowerCase()) ) {
            // a file output: inline the contents of a small structured (text-like) file so the
            // LLM can reason over it, else return the opaque absolute path handle (the
            // opaque-path contract for bulk/binary/oversized data). Defensive fall-back to the
            // path string when the value is unexpectedly not a Path.
            return value instanceof Path
                ? ToolOutputReader.readOrHandle((Path) value, maxInlineBytes)
                : pathString(value)
        }
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
            for( final queue : tool.toolIns )
                queue.bind(PoisonPill.instance)
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
