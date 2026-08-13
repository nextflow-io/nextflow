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

import java.nio.channels.ClosedByInterruptException
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicBoolean

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import io.seqera.npr.api.schema.v1.ModuleMetadata
import nextflow.Global
import nextflow.Nextflow
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.extension.FilesEx
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpec.ModuleParam
import nextflow.script.ChannelOut
import nextflow.script.ProcessConfigV1
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.ProcessEntryHandler
import nextflow.script.params.BaseOutParam
import nextflow.script.params.OutParam
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessOutput

/**
 * Bridges Nextflow processes/modules to the agent as LLM tools that execute as
 * real, request-scoped dataflow nodes.
 *
 * The bridge is built in the workflow body (inside {@link nextflow.script.AgentDef#run},
 * <b>before</b> the dataflow network is ignited). For every tool it derives a portable
 * {@link ToolDescriptor} and starts a dataflow gateway over a request queue. Each tool call
 * carries its own reply variable; the gateway creates fresh input/output channels and invokes
 * a cloned {@link ProcessDef}. Correlation is therefore represented by dataflow variables
 * rather than by ordering on a shared process lane.
 *
 * Two marshalling modes are supported:
 * <ul>
 *   <li><b>scalar</b> (no {@link ModuleSpec}): each typed scalar input param maps to one
 *       request argument; the single scalar output is serialized under its name. Schema from
 *       {@link ProcessToolSchema} (Phase 2 / 3.1).</li>
 *   <li><b>spec-driven</b> (a {@link ModuleSpec}, e.g. from a sibling {@code meta.yml}): one
 *       input queue per spec input channel. The LLM passes FLATTENED args; at dispatch the
 *       components of a tuple input are reassembled IN ORDER into a {@code List}
 *       (file/path → {@link Nextflow#file}, map → the arg Map, scalar → the value) and bound
 *       to that channel. Outputs (tuple or scalar) are serialized back to a JSON object keyed
 *       by emit name; file/path values become absolute path strings (the opaque-path
 *       contract). Schema from {@link ModuleSpecToolSchema} (Phase 3.2).</li>
 * </ul>
 *
 * At tool-call time (post-ignition, on the agent task thread) {@link #call} submits a
 * {@link ToolCall} and waits on its reply. The gateway can service independent requests
 * concurrently; each resulting process still uses the normal executor, cache, work directory,
 * retry and tracing machinery. When the agent's input source is exhausted {@link #close}
 * poisons the request queue so the gateway terminates.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleToolBridge implements ToolDispatcher {

    /**
     * Holds the immutable definition and marshalling metadata for a tool.
     */
    private static class Tool {
        String name
        // scalar mode: the ordered input param names
        List<String> inputParamNames
        String outputParamName
        // spec mode: the module spec
        ModuleSpec spec
        // process template cloned for every request
        ProcessDef processDef
    }

    @TupleConstructor
    private static class ToolCall {
        Tool tool
        Map arguments
        DataflowVariable<String> reply
    }

    /**
     * Per-agent-invocation dispatch context: a sandbox work dir and the readable dirs that the
     * {@code fs:} tools may access. Stored in a {@link ThreadLocal} so the shared,
     * pre-ignition bridge instance is stateless across records: per-record state lives ONLY here,
     * never as an instance field.
     *
     * <p><b>Threading invariant</b>: AiServices dispatches tool calls sequentially on the calling
     * (agent operator) thread — {@code executeToolsConcurrently} is never enabled — so the
     * ThreadLocal value set by {@link #setContext} before a {@code runner.run(req)} call is
     * always the correct context when {@code call()} executes inside that run, regardless of the
     * number of tool calls within a single invocation.
     */
    private static final ThreadLocal<DispatchContext> CONTEXT = new ThreadLocal<DispatchContext>()

    /** Set the per-invocation dispatch context on the current thread. Call before {@code runner.run}. */
    static void setContext(DispatchContext ctx) { CONTEXT.set(ctx) }

    /** Clear the per-invocation dispatch context from the current thread. Call in a finally block after {@code runner.run}. */
    static void clearContext() { CONTEXT.remove() }

    /** Retrieve the per-invocation dispatch context for the current thread. May return {@code null} when called outside a dispatched invocation. */
    private static DispatchContext context() { return CONTEXT.get() }

    private final Map<String,Tool> tools = new LinkedHashMap<String,Tool>()

    private final List<ToolDescriptor> descriptors = new ArrayList<ToolDescriptor>()

    private final DataflowQueue requests = new DataflowQueue()

    private DataflowProcessor gateway

    private final AtomicBoolean closed = new AtomicBoolean()

    /** The wire names of the {@code fs:} tools this bridge serves in the driver JVM, a subset of
     * {@link FilesystemTools#NAMES} — exactly the leaves the agent's {@code tools} directive
     * selected. A name in here is routed to {@link FilesystemTools#call} and has a descriptor;
     * a name NOT in here is not a filesystem tool as far as this bridge is concerned, so a
     * process called {@code read} in an agent that never declared {@code fs:read} still reaches
     * its own module dispatch rather than being hijacked. */
    private final Set<String> filesystemTools

    /** True when the agent selected any {@code fs:} tool, i.e. when the task needs a sandbox
     * {@link DispatchContext}. Read by {@code AgentDef} to decide whether to bind one. */
    boolean isFilesystemEnabled() { !filesystemTools.isEmpty() }

    /** The maximum size (in bytes) of a structured file output whose content is inlined to the
     * LLM; larger or non text-like/binary outputs stay opaque path handles. Defaults to 32 KB. */
    private long maxInlineBytes = ToolOutputReader.DEFAULT_INLINE_BYTES

    /** Set the cap (in bytes) for inlining structured file outputs; non-positive values reset
     * to the 32 KB default. Wired from the {@code agent.maxToolOutputInlineSize} config scope. */
    void setMaxInlineBytes(long bytes) { this.maxInlineBytes = bytes > 0 ? bytes : ToolOutputReader.DEFAULT_INLINE_BYTES }

    /**
     * Build the bridge and its request gateway in the workflow body, before ignition. A tool with a
     * {@link ModuleSpec} (e.g. a sibling {@code meta.yml}) uses spec-driven marshalling; one without
     * uses the scalar typed-I/O path. When a {@link WiredModuleTool} carries a public registry
     * {@link ModuleMetadata}, that metadata is the descriptor source (richer schema); marshalling
     * stays spec-driven. Every wired module is advertised as its OWN tool whose {@code parameters}
     * schema IS that module's flattened input schema, so OpenAI function-calling enforces field
     * names/required-ness per module (an aggregate tool could only use a generic
     * {@code additionalProperties:true} object, which it cannot enforce).
     *
     * <p>The default-valued arg generates the convenience overload (tools-only) via Groovy; the
     * {@code fs:} tools named in {@code filesystemTools} have their calls routed to
     * {@link FilesystemTools#call}. They deliberately do NOT enter {@link #descriptors()} — see
     * {@link #filesystemDescriptors()}.
     *
     * @param wired           the brokered tools the agent declared, in the order they are advertised
     * @param filesystemTools the wire names of the {@code fs:} leaves this bridge must SERVE in the
     *                        driver JVM, a subset of {@link FilesystemTools#NAMES}. Empty for an
     *                        agent that declared no {@code fs:} ref AND for one running on a
     *                        containerized runner, which serves them with its own builtins
     */
    ModuleToolBridge(List<WiredModuleTool> wired,
                     Collection<String> filesystemTools = Collections.<String>emptyList()) {
        // preserve the canonical order of FilesystemTools.NAMES rather than the caller's, so the
        // descriptor list (and hence the tools fingerprint) does not depend on declaration order
        final Collection<String> selectedFs = filesystemTools != null ? filesystemTools : Collections.<String>emptyList()
        this.filesystemTools = new LinkedHashSet<String>(FilesystemTools.NAMES.findAll { selectedFs.contains(it) })
        for( final tool : wired ) {
            if( tool.spec != null )
                wireSpec(tool.name, tool.proc, tool.spec, tool.metadata, tool.nfCore)
            else
                wireScalar(tool.name, tool.proc)
        }
        if( !tools.isEmpty() )
            startGateway()
    }

    private void startGateway() {
        final session = Global.session as Session
        gateway = DataflowHelper.newOperator(
            [inputs: [requests]],
            { ToolCall request -> startInvocation(session, request) } )
    }

    /**
     * Construct each request-scoped graph on the single gateway operator. Script and DAG
     * construction are mutable, whereas the resulting channels are independent and can be
     * collected concurrently.
     */
    private void startInvocation(Session session, ToolCall request) {
        try {
            if( request.tool.spec != null )
                startSpecInvocation(session, request)
            else
                startScalarInvocation(session, request)
        }
        catch( Throwable e ) {
            request.reply.bindError(e)
        }
    }

    /**
     * Pull outputs away from the GPars operator so it remains free to construct the next
     * request and to run the TaskProcessor operators that produce those outputs.
     */
    private void collectAsync(Session session, ToolCall request, Closure<String> collector) {
        // ORCHESTRATION: this thread blocks reading the tool task's output, so it must not
        // come from the execution pool that same tool task needs to run on
        session.getAgentExecService().submit {
            try {
                request.reply.bind(collector.call())
            }
            catch( Throwable e ) {
                request.reply.bindError(e)
            }
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

        // -- portable schema descriptor. The tool name IS the process name: Nextflow process
        //    identifiers are valid OpenAI function-name identifiers (alphanumeric + underscore),
        //    so no sanitization is needed; uniqueness is guaranteed by the in-scope process
        //    namespace (and the `tools` map / `descriptors` keying dedups by name).
        final descriptor = new ToolDescriptor(
            name,
            name,
            ProcessToolSchema.inputSchema(proc),
            ProcessToolSchema.outputSchema(proc) )

        final tool = new Tool(
            name: name,
            inputParamNames: inputParamNames,
            outputParamName: outputParamName,
            processDef: proc )
        tools.put(name, tool)
        // each wired module is advertised as its OWN tool whose parameters schema is the
        // module's flattened input schema, so OpenAI function-calling enforces the field names
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
        // The tool name IS the process/module name (e.g. SKESA): Nextflow process identifiers are
        // valid OpenAI function-name identifiers (alphanumeric + underscore), so no sanitization is
        // needed; uniqueness is guaranteed by the in-scope process namespace (and the `tools` map /
        // `descriptors` keying dedups by name).
        final descriptor = new ToolDescriptor(name, description, inputSchema, null)

        // -- the number of input queues = the PROCESS's declared input-channel count;
        //    `ProcessEntryHandler.getProcessArguments` returns one value per declared input
        //    channel, so the queues must match that, NOT the spec input count. They should be
        //    equal (descriptor and executor are two views of the same meta.yml) -- warn on drift.
        final int nInputs = declaredInputChannelCount(proc)
        if( nInputs != inputs.size() )
            log.warn("Agent tool `${name}`: module spec declares ${inputs.size()} input(s) but the process declares ${nInputs} input channel(s) - using the process input count for marshalling")

        final tool = new Tool(
            name: name,
            spec: spec,
            processDef: proc )
        tools.put(name, tool)
        // each wired module is advertised as its OWN tool whose parameters schema is the
        // module's flattened input schema (required fields + additionalProperties:false), so
        // OpenAI function-calling validates the call against it and the model cannot omit/rename
        // fields — the per-module enforcement that an aggregate module_run tool cannot provide
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
        // Only classic-DSL2 (V1) outputs carry a `topic:` channel name on the param itself.
        // The typed (V2) ProcessOutput.getChannelTopicName() throws UnsupportedOperationException,
        // and V2 topics are a SEPARATE collection (not in getOutputs().getParams()), so a V2 data
        // output is never a topic-source here -- guard on BaseOutParam to avoid the throw.
        return (p instanceof BaseOutParam) && ((BaseOutParam) p).getChannelTopicName() != null
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
        final List<ModuleParam> inputs = spec?.inputs ?: Collections.<ModuleParam>emptyList()
        return inputs
            .collectMany { ModuleParam param -> (param.isTuple() ? param.components : [param]) as List<ModuleParam> }
            .findResults { ModuleParam comp -> comp?.name }
            .toList() as List<String>
    }

    /**
     * The descriptors of the BROKERED tools — the {@code nf:module_run} processes the driver
     * executes as real Nextflow tasks — and nothing else. This list becomes
     * {@code AgentRunnerRequest.toolSpecs}, whose partition invariant is that a runner-native
     * tool never appears in it: a name in {@code toolSpecs} is a name the runner is authorized to
     * call BACK into the driver with, and a {@code fs:}/{@code shell:} tool must never be that.
     * The {@code fs:} leaves this bridge serves in-JVM are exposed separately by
     * {@link #filesystemDescriptors()}.
     */
    List<ToolDescriptor> descriptors() {
        return descriptors
    }

    /**
     * The descriptors of the {@code fs:} leaves this bridge serves in the driver JVM, in the
     * canonical {@link FilesystemTools#NAMES} order. Empty unless the bridge was built for an
     * in-JVM runner with an {@code fs:} selection.
     *
     * <p>Kept OUT of {@link #descriptors()} on purpose: an in-JVM runner advertises these to the
     * model itself (from {@code AgentRunnerRequest.nativeToolNames}) and dispatches them straight
     * back here, so they never travel as brokered descriptors and never enter a broker allowlist.
     */
    List<ToolDescriptor> filesystemDescriptors() {
        return FilesystemTools.descriptors(filesystemTools)
    }

    /**
     * A tool name -> backing process script source map, i.e. the {@code BodyDef.source} that
     * {@code TaskHasher} folds into that process' own task hash. Feeds
     * {@code AgentDef.toolsFingerprint} so that editing a tool's script invalidates the resume
     * cache entry of every agent that can call it. The {@code fs:} tools have no backing
     * process and so do not appear here.
     */
    Map<String,String> toolSources() {
        final result = new LinkedHashMap<String,String>()
        for( final entry : tools.entrySet() )
            result.put(entry.key, entry.value.processDef?.getTaskBody()?.source)
        return result
    }

    /**
     * Execute a tool call. Parse the JSON args, submit a correlated request to the gateway,
     * then block on that request's reply variable.
     *
     * <p><b>Dispatch-level errors are returned as a tool result, not thrown.</b> An unknown
     * tool name, unparseable {@code argsJson}, or a malformed/mis-shaped argument all yield a
     * well-formed {@code {"error": "<message>"}} JSON object that names the tool and what went
     * wrong, so the LLM can see the failure and retry rather than the agent loop being aborted
     * by an exception escaping the dispatcher.
     *
     * <p><b>Task failure is fatal, not recoverable.</b> When the underlying tool <i>process task</i>
     * (the Nextflow task the tool runs as) hard-fails (exit ≠ 0), the session aborts the dataflow
     * network and interrupts the agent operator thread that is blocked on the tool's output channel
     * {@code .val} ({@link groovyx.gpars.dataflow.expression.DataflowExpression#getVal} throws
     * {@link InterruptedException}). This is NOT recoverable by the LLM: the run is already being
     * torn down. Such an {@link InterruptedException} is caught here, the thread's interrupt flag is
     * restored, and it is re-thrown as an {@link AgentToolFatalError} (an {@link Error}, NOT an
     * {@link Exception}) so it escapes langchain4j's tool-execution {@code try/catch(Exception)} and
     * propagates out of {@code agent.chat(...)} on the agent task-body thread and through
     * {@code TaskProcessor}'s exec-body failure handling, which aborts the run cleanly — instead of
     * being fed back to the model as an error tool result that loops to the iteration cap.
     *
     * @param toolName the name of the tool to invoke
     * @param argsJson the LLM-supplied arguments as a JSON object string
     * @return the tool result serialized as a JSON object string; on a dispatch-level failure
     *         a {@code {"error": "..."}} JSON object so the LLM can recover
     * @throws AgentToolFatalError when the underlying process task fails and the session aborts
     *         (the blocking output pull is interrupted) — the run is torn down, not retried
     */
    @Override
    String call(String toolName, String argsJson) {
        try {
            // Route the selected fs: tools BEFORE the module lookup so that a bridge with no
            // wired modules still dispatches them. Only the SELECTED names are routed here:
            // an unselected `read`/`find`/... is not a filesystem tool at all and falls through
            // to the module lookup, so a process of that name stays reachable.
            if( filesystemTools.contains(toolName) )
                return FilesystemTools.call(toolName, parseArgs(toolName, argsJson), context(), maxInlineBytes)

            final tool = tools.get(toolName)
            if( tool == null )
                throw new IllegalArgumentException("Unknown agent tool `${toolName}` - available tools: ${tools.keySet()}")

            final DataflowVariable<String> reply = new DataflowVariable<String>()
            if( closed.get() )
                throw new IllegalStateException("Agent tool bridge is closed")
            final parsed = parseArgs(toolName, argsJson)
            requests.bind(new ToolCall(tool, parsed, reply))
            final result = reply.val
            // after the module task completes, scan the result for file path strings and whitelist
            // their parent dirs in the dispatch context so the `fs:` tools can read module outputs.
            // ONLY whitelist on a non-error result: a failed module returns {"error":"..."}, and if
            // that error message contains an absolute path the parent dir must NOT be whitelisted —
            // that would silently widen the sandbox with data from an error string, not an output.
            final resultParsed = parseResultJson(result)
            if( !isErrorResult(resultParsed) )
                whitelistOutputDirs(resultParsed)
            return result
        }
        catch( Exception e ) {
            // A genuine underlying-task failure / session abort surfaces as an InterruptedException —
            // which GPars WRAPS in a RuntimeException (e.g. RuntimeException(cause=InterruptedException)
            // from DataflowStreamReadAdapter.getVal), so a bare `catch(InterruptedException)` misses it.
            // That case is FATAL: the run is being torn down — do NOT swallow it into a {"error":...}
            // tool result (that feeds the failure back to the model and loops to maxIterations).
            // Detect it via the cause chain, restore the interrupt flag, and re-throw as an Error so it
            // escapes langchain4j's tool-execution catch(Exception) and aborts the run cleanly.
            if( isInterrupted(e) ) {
                Thread.currentThread().interrupt()
                log.debug("Agent tool `${toolName}` aborted by session interrupt (underlying task failure)", e)
                throw new AgentToolFatalError("Agent tool `${toolName}` aborted - the underlying process task failed", e)
            }
            // dispatch-level failure (unknown tool / arg parsing / arg marshalling): return it
            // as a tool result so the LLM can recover, rather than letting it abort the loop
            final message = "Agent tool `${toolName}` failed - ${e.message ?: e.toString()}".toString()
            log.warn(message, e)
            return JsonOutput.toJson([error: message])
        }
    }

    /**
     * True when {@code t} is, or wraps anywhere in its cause chain, an {@link InterruptedException}.
     * A failed underlying task makes the session abort and interrupt this operator's blocking output
     * pull; GPars surfaces that as a {@link RuntimeException} wrapping the {@link InterruptedException}
     * (e.g. {@code DataflowStreamReadAdapter.getVal}), so a bare {@code instanceof} check misses it.
     */
    private static boolean isInterrupted(Throwable t) {
        for( Throwable c = t; c != null; c = c.getCause() ) {
            // ClosedByInterruptException is the IOException form raised when an abort interrupts
            // in-flight (interruptible) file I/O — also fatal, also outside InterruptedException
            if( c instanceof InterruptedException || c instanceof ClosedByInterruptException )
                return true
        }
        return false
    }

    /**
     * Predicate: true when the parsed result is an error-shaped object (a Map containing
     * an {@code "error"} key). Used to guard the {@code whitelistOutputDirs} call in
     * {@link #call} so that absolute paths appearing only in error messages
     * are not added to the filesystem sandbox whitelist.
     *
     * @param parsed the already-parsed result value (Map, List, String, null, …)
     * @return true when {@code parsed} is a Map with an {@code "error"} key, false otherwise
     */
    @PackageScope
    static boolean isErrorResult(Object parsed) {
        return parsed instanceof Map && ((Map) parsed).containsKey('error')
    }

    /**
     * Parse a JSON result string silently. Returns the parsed object (Map, List, String, …)
     * or {@code null} if the string is null/blank/unparseable.
     */
    private static Object parseResultJson(String resultJson) {
        if( resultJson == null || resultJson.isBlank() )
            return null
        try {
            return new JsonSlurper().parseText(resultJson)
        }
        catch( Exception e ) {
            log.trace("Agent tool result is not parsable JSON, treating as no result: ${e.message}")
            return null
        }
    }

    /**
     * Scan an already-parsed JSON result value for absolute file path strings and add each
     * file's parent directory to the dispatch context's readable-dirs whitelist. Called after a
     * module tool task completes successfully (callers must skip error results).
     * No-op when no dispatch context is active or the parsed value is null.
     *
     * <p>Package-visible for unit testing; not part of the public API.</p>
     */
    @PackageScope
    static void whitelistOutputDirs(Object parsed) {
        final DispatchContext ctx = context()
        if( ctx == null || parsed == null )
            return
        collectPathsFromValue(parsed, ctx)
    }

    private static void collectPathsFromValue(Object value, DispatchContext ctx) {
        if( value instanceof Map ) {
            ((Map) value).values().each { collectPathsFromValue(it, ctx) }
        }
        else if( value instanceof List ) {
            ((List) value).each { collectPathsFromValue(it, ctx) }
        }
        else if( value instanceof String ) {
            final s = (String) value
            if( s.startsWith('/') ) {
                try {
                    final p = Path.of(s)
                    // Whitelist the produced path ITSELF, never its parent. The parent is every
                    // sibling of the output too, and for an output that lands outside the work tree
                    // — a `publishDir` target, a path on a shared filesystem — that is a directory
                    // of content no cache key covers, which is what forces `cache false` on a
                    // filesystem agent. A file entry matches only that file; a directory output
                    // still admits its contents, since containment is by path prefix.
                    // Guarded on a REAL produced path, so an arbitrary path-shaped string in the
                    // tool result grants nothing, and on a non-empty name count, so `/` cannot be
                    // whitelisted.
                    if( p.getNameCount() > 0 && Files.exists(p) )
                        ctx.addReadablePath(p)
                }
                catch( Exception e ) {
                    log.trace("Agent filesystem whitelist: skipping non-path output value `${s}`: ${e.message}")
                }
            }
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

    private void startScalarInvocation(Session session, ToolCall request) {
        final tool = request.tool
        final parsed = request.arguments
        final args = tool.inputParamNames.collect { parsed.get(it) }
        final invocation = tool.processDef.clone()
        final ChannelOut out = (ChannelOut) invocation.run(args as Object[])
        final DataflowReadChannel resultChannel = CH.getReadChannel(out[0])

        collectAsync(session, request) {
            final result = resultChannel.val
            final key = (tool.outputParamName == '$out') ? 'result' : tool.outputParamName
            return JsonOutput.toJson([(key): result])
        }
    }

    private void startSpecInvocation(Session session, ToolCall request) {
        final tool = request.tool
        final parsed = request.arguments
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
        final ProcessDef invocation = tool.processDef.clone()
        final ChannelOut channelOut = (ChannelOut) invocation.run(args as Object[])
        // Capture all readers before pulling any value. This matters if a future component
        // invocation returns broadcast outputs rather than singleton variables.
        final List<DataflowReadChannel> outReadChannels = (0..<channelOut.size()).collect { CH.getReadChannel(channelOut[it]) }

        collectAsync(session, request) {
            // -- block on each output channel, in declared order, serializing the records
            final result = new LinkedHashMap<String,Object>()
            final outputs = spec.outputs ?: Collections.<ModuleParam>emptyList()
            final emitNames = emitNamesByPosition(channelOut)
            final int nOut = channelOut.size()
            // The ProcessDef's declared output params are positionally aligned with ChannelOut.
            final List<OutParam> outParams = declaredOutputParams(invocation)
            for( int i=0; i<outputs.size(); i++ ) {
                final param = outputs[i]
                // Topic/eval outputs are bookkeeping and never bind a per-invocation value.
                if( isEvalOutput(param) || isTopicOutput(outParams, i) )
                    continue
                if( i >= nOut )
                    continue
                final key = outputKey(param, emitNames, i)
                final value = outReadChannels[i].val
                result.put(key, serializeOutput(param, value))
            }
            return JsonOutput.toJson(result)
        }
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
        return 'eval'.equalsIgnoreCase(param.type) ||
            (param.components?.any { it != null && 'eval'.equalsIgnoreCase(it.type) } ?: false)
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
        if( ToolSchema.isFileType(comp.type?.toLowerCase()) ) {
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
            return FilesEx.toUriString(((Path) value).toAbsolutePath())
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
     * Poison the gateway request queue so its operator terminates. Idempotent.
     */
    void close() {
        // close() is triggered only after the agent output completes, so no valid call can
        // race with the poison pill. The atomic flag provides idempotence and visibility.
        if( !closed.compareAndSet(false, true) )
            return
        if( gateway != null )
            requests.bind(PoisonPill.instance)
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
