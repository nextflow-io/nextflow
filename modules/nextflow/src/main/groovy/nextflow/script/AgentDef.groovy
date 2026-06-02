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
import nextflow.Session
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

    /** Maps a compiled external-module {@link ProcessDef} back to its source file path,
     * so the sibling {@code meta.yml} can be located (Phase 3.2). */
    private transient Map<ProcessDef,Path> compiledModulePaths = new IdentityHashMap<ProcessDef,Path>()

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
        // Phase 2 is "tools XOR structured output": when tools are declared the plugin's
        // tool loop drives the conversation and ignores any responseFormat, so the final
        // answer is free text. Binding that text to a record (structured) output would fail
        // at JSON-parse time - reject the combination up front with a clear message.
        if( this.tools && structured )
            throw new ScriptRuntimeException("Agent `${name}`: combining tools with a record (structured) output is not yet supported - use a plain output type (e.g. String) when declaring tools")
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
     * Resolve the agent's declared {@code tools} entries to their {@link ProcessDef}s and
     * pre-wire them as agent tools. Each entry is a String that is resolved as follows:
     * <ul>
     *   <li>an <b>in-scope process name</b> (resolvable via {@link ScriptMeta#getProcess})
     *       is used directly (Phase 2 behavior);</li>
     *   <li>a <b>local module file path</b> (ends with {@code .nf}, or starts with {@code ./},
     *       {@code ../} or {@code /}) is compiled to a {@link ProcessDef} at agent-run time
     *       on the owner's live {@link nextflow.Session} (Phase 3.1);</li>
     *   <li>a <b>registry reference</b> (e.g. {@code nf-core/fastqc}) is not yet supported
     *       (Phase 3.3) and fails loudly.</li>
     * </ul>
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
        final specs = new LinkedHashMap<String, nextflow.module.ModuleSpec>()
        for( final entry : declared ) {
            final toolRef = entry?.toString()
            // 1) an in-scope process name (Phase 2, unchanged)
            final inScope = meta?.getProcess(toolRef)
            if( inScope != null ) {
                resolved.put(toolRef, inScope)
                continue
            }
            // 2) a local module file path -> compile to a ProcessDef (Phase 3.1);
            //    when a sibling meta.yml is present, drive the schema/marshalling from
            //    its ModuleSpec (Phase 3.2)
            if( looksLikeModulePath(toolRef) ) {
                final entry2 = resolveModuleTool(meta, toolRef)
                final toolName = entry2.key
                resolved.put(toolName, entry2.value)
                final spec = loadSiblingSpec(entry2.value)
                if( spec != null )
                    specs.put(toolName, spec)
                continue
            }
            // 3) a registry reference (e.g. `nf-core/fastqc`) -> not yet supported (Phase 3.3)
            if( toolRef ==~ /[^\/\s]+\/[^\/\s]+/ )
                throw new ScriptRuntimeException("Agent `${name}` tool `${toolRef}`: registry module tools are not yet supported (Phase 3.3); use an in-scope process or a local .nf module file")
            // 4) otherwise: not a process in scope, not a path
            throw new ScriptRuntimeException("Agent `${name}` tool `${toolRef}` is not a process in scope")
        }
        return new ModuleToolBridge(resolved, specs)
    }

    /**
     * Look for a sibling {@code meta.yml} / {@code meta.yaml} next to the compiled module file
     * and, when present, load it into a {@link nextflow.module.ModuleSpec} to drive spec-driven
     * tool schema and tuple/path/map marshalling (Phase 3.2). Returns {@code null} when no
     * sibling spec is found.
     */
    @CompileDynamic
    private nextflow.module.ModuleSpec loadSiblingSpec(ProcessDef proc) {
        final modPath = compiledModulePaths.get(proc)
        if( modPath == null )
            return null
        final dir = modPath.getParent()
        for( final candidate : ['meta.yml', 'meta.yaml'] ) {
            final specPath = dir.resolve(candidate)
            if( specPath.toFile().exists() )
                return nextflow.module.ModuleSpecFactory.fromYaml(specPath)
        }
        return null
    }

    /**
     * Whether a {@code tools} entry looks like a local module file path: it ends with
     * {@code .nf}, or starts with {@code ./}, {@code ../} or {@code /}.
     */
    private static boolean looksLikeModulePath(String ref) {
        if( !ref )
            return false
        return ref.endsWith('.nf') || ref.startsWith('./') || ref.startsWith('../') || ref.startsWith('/')
    }

    /**
     * Compile a local module file into a runnable {@link ProcessDef} on the owner's live
     * {@link nextflow.Session} (pre-ignition). The path is resolved relative to the owner
     * script's directory when relative; the compiled module's single process becomes the
     * tool (the tool name is the process name).
     *
     * @return a {@code toolName -> ProcessDef} entry
     */
    @CompileDynamic
    private Map.Entry<String,ProcessDef> resolveModuleTool(ScriptMeta meta, String ref) {
        final session = Global.session as Session
        if( session == null )
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: cannot compile a module file - no active session")

        // -- resolve the module path relative to the owner script's dir when relative
        final base = ownerBaseDir(meta, session)
        Path modPath = Path.of(ref)
        if( !modPath.isAbsolute() )
            modPath = base.resolve(ref).normalize()
        if( !modPath.toFile().exists() )
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: module file not found - resolved to `${modPath}`")

        // -- compile the module on the live session (spike-validated recipe);
        //    setModule(true) is mandatory so the standalone process does not auto-build
        //    an entry workflow
        final BaseScript modScript
        try {
            final loader = new nextflow.script.parser.v2.ScriptLoaderV2(session)
            loader.setModule(true)
            loader.parse(modPath)
            loader.runScript()
            modScript = loader.getScript()
        }
        catch( Exception e ) {
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: failed to compile module file `${modPath}` - ${e.message}", e)
        }

        // -- pick the module's single process; the tool name is the process name
        final modMeta = ScriptMeta.get(modScript)
        final procNames = new ArrayList<String>(modMeta.getProcessNames())
        if( procNames.isEmpty() )
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: module file `${modPath}` declares no process")
        if( procNames.size() > 1 )
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: module file `${modPath}` declares more than one process (${procNames.sort()}); naming a specific process is not yet supported")
        final procName = procNames[0]
        final proc = modMeta.getProcess(procName)
        // remember the source path so a sibling meta.yml can be located (Phase 3.2)
        compiledModulePaths.put(proc, modPath)
        return new AbstractMap.SimpleImmutableEntry<String,ProcessDef>(procName, proc)
    }

    /**
     * The directory used to resolve relative module-file tool paths: the owner script's
     * directory when known, otherwise the session base dir, otherwise the launch (current)
     * directory.
     */
    private static Path ownerBaseDir(ScriptMeta meta, Session session) {
        final moduleDir = meta?.getModuleDir()
        if( moduleDir != null )
            return moduleDir
        if( session?.baseDir != null )
            return session.baseDir
        return Path.of('.').toAbsolutePath().normalize()
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
