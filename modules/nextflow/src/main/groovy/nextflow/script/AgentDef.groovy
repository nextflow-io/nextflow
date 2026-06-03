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
import nextflow.agent.AgentConfig
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

    /**
     * Read the {@code agent} config scope from the live {@link nextflow.Session}.
     * Returns an empty (all-null) {@link AgentConfig} when there is no active
     * session or the scope is not defined, so callers can always read defaults.
     */
    protected AgentConfig agentConfig() {
        final session = Global.session as Session
        final opts = (session?.config?.get('agent') as Map) ?: Collections.emptyMap()
        return new AgentConfig(opts)
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
        // read the `agent` config scope from the live session; the agent's own
        // directives always take precedence - the scope only fills in defaults
        // for a value the agent did not declare
        final agentConfig = agentConfig()
        // effective model: directive first, else the configured default; when both
        // are null the existing required-model error fires downstream (in the runner)
        final agentModel = this.model ?: agentConfig.defaultModel
        final agentInstruction = this.instruction
        final agentTools = this.tools
        // effective max iterations: directive first, else the configured default, else 20
        final agentMaxIter = (this.maxIterations != null
            ? this.maxIterations
            : (agentConfig.maxIterationsDefault != null ? agentConfig.maxIterationsDefault : 20)) as int
        // effective request timeout (seconds): the configured value, else 120
        final agentTimeoutSecs = (agentConfig.requestTimeout != null
            ? agentConfig.requestTimeout.seconds
            : 120L) as int
        final promptDef = this.prompt

        // resolve declared `tools` to in-scope processes and pre-wire them into the
        // dataflow network (before ignition); the bridge is poisoned on completion
        final ModuleToolBridge bridge = createToolBridge()
        if( bridge != null )
            bridge.setMaxInlineBytes(agentConfig.maxToolOutputInlineBytes())
        final List<ToolDescriptor> toolSpecs = bridge != null ? bridge.descriptors() : null

        final mapper = { Object item ->
            final cl = (Closure) promptDef.closure.clone()
            cl.setDelegate([(inputName): item])
            cl.setResolveStrategy(Closure.DELEGATE_FIRST)
            final promptText = cl.call()?.toString()
            final inputJson = toJson(item)
            final req = new AgentRunnerRequest(agentModel, agentInstruction, promptText, agentMaxIter, agentTools, outputSchema, inputJson, toolSpecs, (bridge as ToolDispatcher), agentTimeoutSecs)
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
     *   <li>a <b>registry reference</b> (e.g. {@code nf-core/fastqc}, optionally {@code @version})
     *       is resolved through the module registry via {@link nextflow.module.ModuleResolver}
     *       (local install first, else auto-install/download) and compiled to a {@link ProcessDef}
     *       at agent-run time (Phase 3.3).</li>
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
        // registry-sourced tool metadata (Phase 3.3 / Task 3): when a registry tool resolves and
        // its public `ModuleMetadata` is fetched, the descriptor (description + input schema) is
        // sourced from it; the value is wrapped together with the nf-core scope flag (for the
        // `meta.id` convention). Tools missing here fall back to the sibling meta.yml ModuleSpec.
        final metadatas = new LinkedHashMap<String, ModuleToolBridge.RegistryMeta>()
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
            // 3) a registry reference (e.g. `nf-core/fastqc`, optionally `@version`) -> resolve
            //    via the module registry (ModuleResolver), compile the resolved main.nf to a
            //    ProcessDef (the Phase 3.1 recipe) and, when an installed sibling meta.yml is
            //    present, drive the schema/marshalling from its ModuleSpec (Phase 3.2). The
            //    descriptor is sourced from the public registry `ModuleMetadata` when available.
            if( looksLikeRegistryRef(toolRef) ) {
                final entry3 = resolveRegistryTool(toolRef, metadatas)
                final toolName = entry3.key
                resolved.put(toolName, entry3.value)
                final spec = loadSiblingSpec(entry3.value)
                if( spec != null )
                    specs.put(toolName, spec)
                continue
            }
            // 4) otherwise: not a process in scope, not a path, not a registry ref
            throw new ScriptRuntimeException("Agent `${name}` tool `${toolRef}` is not a process in scope")
        }
        return new ModuleToolBridge(resolved, specs, metadatas)
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

        // -- compile the resolved module file to a ProcessDef; the tool name is the process name
        final proc = compileModuleProcess(ref, modPath, session)
        return new AbstractMap.SimpleImmutableEntry<String,ProcessDef>(proc.getName(), proc)
    }

    /**
     * Compile a resolved module {@code main.nf} into a runnable {@link ProcessDef} on the
     * owner's live {@link nextflow.Session} (pre-ignition). Shared by the local-file (Phase 3.1)
     * and registry (Phase 3.3) tool-resolution paths.
     *
     * <p>{@code setModule(true)} is mandatory so the standalone process does not auto-build an
     * entry workflow. The resolved {@code modPath} is remembered (keyed by the returned
     * {@link ProcessDef}) so a sibling {@code meta.yml} can be located later (Phase 3.2).
     *
     * @param ref     the original {@code tools} entry (used for error messages)
     * @param modPath the resolved module {@code main.nf} path
     * @param session the live session
     * @return the module's single {@link ProcessDef}
     */
    @CompileDynamic
    private ProcessDef compileModuleProcess(String ref, Path modPath, Session session) {
        final BaseScript modScript
        try {
            final loader = new nextflow.script.parser.v2.ScriptLoaderV2(session)
            // -- isolate the module's binding (seeded with the owner's params) so loading it
            //    does not overwrite the shared session binding's `moduleDir` (and `container`/
            //    `conda` resolved lazily against it). Mirrors IncludeDef.loadModuleV2, which runs
            //    each included module on its own ScriptBinding. Without this, declaring multiple
            //    file/registry tools makes every tool resolve `moduleDir` to the LAST-compiled
            //    module dir (container bleed).
            final modBinding = new ScriptBinding()
            final ownerParams = owner?.getBinding()?.getParams()
            if( ownerParams != null )
                modBinding.setParams(ownerParams)
            loader.setMainBinding(modBinding)
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
        final proc = modMeta.getProcess(procNames[0])
        // remember the source path so a sibling meta.yml can be located (Phase 3.2)
        compiledModulePaths.put(proc, modPath)
        return proc
    }

    /**
     * Whether a {@code tools} entry looks like a module <b>registry</b> reference, i.e. a
     * {@code scope/name} (optionally suffixed with {@code @version}) that parses as a
     * {@link nextflow.module.ModuleReference}. A leading path (handled by
     * {@link #looksLikeModulePath}) is excluded by checking the path form first at the call site.
     */
    private static boolean looksLikeRegistryRef(String ref) {
        if( !ref )
            return false
        final refNoVer = stripVersion(ref).first
        try {
            nextflow.module.ModuleReference.parse(refNoVer)
            return true
        }
        catch( Exception ignored ) {
            return false
        }
    }

    /** Split a {@code scope/name@version} ref into {@code [scope/name, version|null]}. */
    private static Tuple2<String,String> stripVersion(String ref) {
        final at = ref.indexOf('@')
        if( at < 0 )
            return new Tuple2<String,String>(ref, null)
        return new Tuple2<String,String>(ref.substring(0, at), ref.substring(at + 1) ?: null)
    }

    /**
     * Resolve a registry module reference (e.g. {@code nf-core/fastqc} or
     * {@code nf-core/fastqc@1.0.0}) to a runnable {@link ProcessDef} (Phase 3.3).
     *
     * <p>The reference is resolved through the module registry via {@link nextflow.module.ModuleResolver}
     * built from the session base dir and the {@code registry} config scope (mirroring
     * {@code CmdModuleRun.resolveRemoteModule}). {@code resolve(ref, version, autoInstall=true)}
     * checks the local install FIRST (no network when already installed) and otherwise downloads
     * the module. The resolved {@code main.nf} is then compiled via {@link #compileModuleProcess};
     * the tool name is the sanitized reference (so it is a valid LLM-facing identifier).
     *
     * @return a {@code sanitizedToolName -> ProcessDef} entry
     */
    @CompileDynamic
    private Map.Entry<String,ProcessDef> resolveRegistryTool(String ref, Map<String,ModuleToolBridge.RegistryMeta> metadatas) {
        final session = Global.session as Session
        if( session == null )
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: cannot resolve a registry module - no active session")

        final parts = stripVersion(ref)
        final refNoVer = parts.first
        final version = parts.second

        final nextflow.module.ModuleReference moduleRef
        try {
            moduleRef = nextflow.module.ModuleReference.parse(refNoVer)
        }
        catch( Exception e ) {
            throw new ScriptRuntimeException("Agent `${name}` tool `${ref}`: invalid registry module reference - ${e.message}", e)
        }

        // -- the modules/ dir lives at the project (base) root; mirror CmdModuleRun
        final baseDir = session.baseDir ?: Path.of('.').toAbsolutePath().normalize()
        final registryConfig = new nextflow.config.RegistryConfig((session.config?.registry as Map) ?: Collections.emptyMap())

        // -- resolve (local install first, else auto-install via the registry client)
        final client = nextflow.module.RegistryClientFactory.forConfig(registryConfig)
        final Path modPath
        try {
            final resolver = new nextflow.module.ModuleResolver(baseDir, client)
            modPath = resolver.resolve(moduleRef, version, /*autoInstall*/ true)
        }
        catch( Exception e ) {
            throw new ScriptRuntimeException("Unable to resolve agent tool module `${ref}` from the registry: ${e.message}. Check registry access/credentials.", e)
        }
        if( modPath == null || !modPath.toFile().exists() )
            throw new ScriptRuntimeException("Unable to resolve agent tool module `${ref}` from the registry: resolved path `${modPath}` does not exist. Check registry access/credentials.")

        // -- compile the resolved main.nf to a ProcessDef (the Phase 3.1 recipe)
        final proc = compileModuleProcess(ref, modPath, session)
        // the LLM-facing tool name must be a valid identifier: sanitize `/`, `-`, `@`, `.` to `_`
        final toolName = sanitizeToolName(refNoVer)

        // -- fetch the PUBLIC registry metadata (description + input schema source). Reuse the
        //    SAME RegistryClient (CmdModuleView pattern): `getModule(fullName).latest?.metadata`,
        //    or `getModuleRelease(fullName, version)?.metadata` when a version is pinned. A registry
        //    hiccup / private module degrades to null -> the bridge falls back to the meta.yml spec.
        final metadata = fetchModuleMetadata(client, moduleRef, version)
        if( metadata != null )
            metadatas.put(toolName, new ModuleToolBridge.RegistryMeta(metadata, moduleRef.scope == 'nf-core'))

        return new AbstractMap.SimpleImmutableEntry<String,ProcessDef>(toolName, proc)
    }

    /**
     * Fetch the public {@link io.seqera.npr.api.schema.v1.ModuleMetadata} for a resolved registry
     * module via the SAME {@link io.seqera.npr.client.RegistryClient} already built to resolve the
     * module ({@code GET /api/v1/modules/{name}} is anonymous/public). Any failure (network hiccup,
     * private module without metadata, missing release) is logged and yields {@code null} so the
     * caller degrades gracefully to the sibling {@code meta.yml} {@link nextflow.module.ModuleSpec}.
     */
    @CompileDynamic
    private io.seqera.npr.api.schema.v1.ModuleMetadata fetchModuleMetadata(io.seqera.npr.client.RegistryClient client, nextflow.module.ModuleReference moduleRef, String version) {
        try {
            if( version )
                return client.getModuleRelease(moduleRef.fullName, version)?.metadata
            return client.getModule(moduleRef.fullName)?.latest?.metadata
        }
        catch( Exception e ) {
            log.warn("Agent `${name}` tool `${moduleRef.fullName}`: could not fetch registry metadata (${e.message}) - falling back to the local module spec for the tool descriptor")
            return null
        }
    }

    /**
     * Turn a module reference into a valid identifier usable as an LLM-facing tool name by
     * replacing any non-word character (e.g. {@code /}, {@code -}, {@code .}) with {@code _}.
     */
    private static String sanitizeToolName(String ref) {
        return ref.replaceAll(/[^A-Za-z0-9_]/, '_')
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
