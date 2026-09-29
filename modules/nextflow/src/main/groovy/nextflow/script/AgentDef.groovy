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
import java.util.regex.Pattern

import groovy.json.JsonOutput
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel

import nextflow.Global
import nextflow.Session
import nextflow.agent.AgentCallInfo
import nextflow.agent.AgentConfig
import nextflow.agent.AgentLaunchConditions
import nextflow.agent.AgentLaunchSpec
import nextflow.agent.AgentOutputMode
import nextflow.agent.AgentOutputPlan
import nextflow.agent.rpc.AgentRpcHost
import nextflow.agent.rpc.AgentRpcRegistration
import nextflow.agent.AgentRunner
import nextflow.agent.AgentRunnerProvider
import nextflow.agent.AgentRunnerRequest
import nextflow.agent.AgentTaskInfo
import nextflow.agent.DispatchContext
import nextflow.agent.ModuleToolBridge
import nextflow.agent.ModuleToolResolver
import nextflow.agent.RecordSchema
import nextflow.agent.SkillDescriptor
import nextflow.agent.SkillResolver
import nextflow.agent.SkillResource
import nextflow.agent.ToolDescriptor
import nextflow.agent.ToolDispatcher
import nextflow.agent.ToolRefResolver
import nextflow.agent.ToolSchema
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.FilesEx
import nextflow.plugin.Plugins
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskConfig
import nextflow.util.CacheHelper
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.processor.TaskPath
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.dsl.ProcessConfigBuilder
import nextflow.script.params.v2.ProcessFileInput
import nextflow.script.params.v2.ProcessFileOutput
import nextflow.script.types.Record
import org.pf4j.PluginWrapper

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

    /** The runner named in the error raised when {@code shell:} is declared against a runner
     * that cannot serve it. The family needs the runner's own container boundary, which only a
     * canonical launch-spec runner provides. */
    private static final String SHELL_FAMILY_RUNNER = 'pi'

    /**
     * What a wire name — {@code ToolDescriptor.name}, i.e. the name the LLM sees — must match:
     * the OpenAI function-name charset (§4). Deliberately checked rather than sanitized, because
     * rewriting `my$proc` to `my_proc` would silently merge it with a process already called
     * that; see {@link #checkWireNames}.
     */
    private static final Pattern WIRE_NAME = ~/[a-zA-Z0-9_-]{1,64}/

    /** Max length of a wire name, restated so the error can quote the number it enforces. */
    private static final int WIRE_NAME_MAX = 64

    /** Tool the runner injects when the agent declares `skills` (both runners). */
    private static final String SKILL_ACTIVATE_TOOL = 'activate_skill'

    /** Its companion, injected under the same condition. */
    private static final String SKILL_RESOURCE_TOOL = 'read_skill_resource'

    /** Tool the canonical runner injects when the agent declares a structured output. */
    private static final String FINAL_ANSWER_TOOL = 'final_answer'

    /** Iteration cap applied when neither the agent nor the `agent` scope declares one. Restated
     * in {@link nextflow.agent.AgentConfig}'s `maxIterations` description as the user-facing default. */
    private static final int DEFAULT_MAX_ITERATIONS = 20

    /** Per-request timeout, in seconds, applied when the `agent` scope declares none. Restated in
     * {@link nextflow.agent.AgentConfig}'s `requestTimeout` description as the user-facing default. */
    private static final long DEFAULT_REQUEST_TIMEOUT_SECONDS = 120L

    private BaseScript owner
    private String name
    private String simpleName
    /**
     * The name the agent was declared with in its own script or module. Unlike {@link #name}
     * and {@link #simpleName} it is stable and cannot be changed by aliasing the include, so
     * a `withName:'<declared name>'` config selector matches an aliased agent — mirroring
     * {@link ProcessDef#baseName}.
     */
    private String baseName
    private Map<String,Object> directives
    private List<AgentInput> inputs
    private List<AgentOutput> outputs
    private PromptDef prompt
    /**
     * The implicit file stagers the compiler inferred from the declared input types, replayed
     * into every invocation's {@link ProcessConfigV2}. A {@link ProcessFileInput} is stateless
     * under {@code resolve} (it clones the closure per task), so holding one here has exactly
     * the lifetime a process's has.
     */
    private List<ProcessFileInput> fileInputs
    /** The unstagers inferred from `file(...)`/`files(...)` in the output expressions. */
    private Map<String,ProcessFileOutput> fileOutputs

    /** Immutable result of lowering an agent before the processor is started. */
    private static class BuiltAgentTask {
        final TaskProcessor processor
        final ModuleToolBridge bridge

        BuiltAgentTask(TaskProcessor processor, ModuleToolBridge bridge) {
            this.processor = processor
            this.bridge = bridge
        }
    }

    /** The runner picked for an agent, with the two facts read off it once. */
    private static class SelectedRunner {
        final AgentRunner runner
        /** Its stable user-facing identifier, which is part of the cache identity. */
        final String name
        /** Its canonical launch description, or {@code null} for a legacy in-JVM runner. */
        final AgentLaunchSpec launchSpec

        SelectedRunner(AgentRunner runner, String name, AgentLaunchSpec launchSpec) {
            this.runner = runner
            this.name = name
            this.launchSpec = launchSpec
        }
    }

    /** The agent's process config, with the executor it was resolved to. */
    private static class ResolvedProcessConfig {
        final ProcessConfigV2 config
        /**
         * The RESOLVED executor, read off the config ONCE. It travels beside the config rather than
         * being re-read downstream so the admission check ({@link AgentDef#buildProcessConfig}, which
         * decides whether this runner may be offloaded at all) and the launch check
         * ({@link AgentDef#resolveLaunch}) cannot come to disagree about which executor was admitted
         * -- which is what a step inserted between them that touched `executor` would otherwise do.
         */
        final String executor

        ResolvedProcessConfig(ProcessConfigV2 config, String executor) {
            this.config = config
            this.executor = executor
        }
    }

    /** The resolved tool selection, and everything lowered from it before ignition. */
    private static class ResolvedTools {
        /** The declared `tools` refs expanded into the two halves of the §5 partition. */
        final ToolRefResolver.Selection selection
        /** The brokered tools' dataflow request gateway; {@code null} when nothing is brokered. */
        final ModuleToolBridge bridge
        /** The BROKERED half of the partition, as descriptors. */
        final List<ToolDescriptor> toolSpecs
        /** The RUNNER-NATIVE half of the partition, as bare wire names. */
        final List<String> nativeToolNames
        /** Whether the body must build a per-task sandbox context for the in-JVM `fs:` tools. */
        final boolean needsSandbox

        ResolvedTools(ToolRefResolver.Selection selection, ModuleToolBridge bridge,
                List<ToolDescriptor> toolSpecs, List<String> nativeToolNames, boolean needsSandbox) {
            this.selection = selection
            this.bridge = bridge
            this.toolSpecs = toolSpecs
            this.nativeToolNames = nativeToolNames
            this.needsSandbox = needsSandbox
        }
    }

    /** Effective, immutable values shared by the canonical and in-JVM request paths. */
    private static class ResolvedAgentSettings {
        final String model
        final String instruction
        final String goal
        final int maxIterations
        final int requestTimeoutSeconds
        final boolean trace
        final String agentName
        final List tools
        final List<ToolDescriptor> toolSpecs
        /** The RUNNER-NATIVE half of the resolved selection (§5): the bare wire names of the
         * {@code fs:}/{@code shell:} leaves, travelling BESIDE {@code toolSpecs} and never inside
         * it. A containerized runner enables its own builtins from these; an in-JVM runner turns
         * them into descriptors itself and dispatches them back through the bridge. */
        final List<String> nativeToolNames
        final List<SkillDescriptor> skills
        final Map outputSchema
        /** LLM provider credential resolved by the core ladder; NEVER folded into the cache
         * key and never written to the canonical source or the task info. */
        final String apiKey
        /** OpenAI-compatible endpoint resolved by the core ladder; part of the agent's
         * identity, so it DOES enter the cache key. */
        final String baseUrl
        /** Provider namespace the pair above was resolved from; carried so a runner can name the
         * variables the ladder actually consulted when it has to report a missing credential. */
        final String apiProvider
        /** Whether a provider credential resolved and was withheld by the endpoint gate -- NOT the
         * same as none resolving; see {@link nextflow.agent.AgentConfig#credentialWithheldFor}. */
        final boolean credentialWithheld
        /** The address THIS agent's task dials the driver's broker on, resolved by the pre-ignition
         * guard. Null on the in-JVM path, which dials nothing. */
        final AgentRpcHost brokerHost

        ResolvedAgentSettings(String model, String instruction, String goal, int maxIterations,
                int requestTimeoutSeconds, boolean trace, String agentName, List tools,
                List<ToolDescriptor> toolSpecs, List<String> nativeToolNames,
                List<SkillDescriptor> skills, Map outputSchema,
                String apiKey, String baseUrl, String apiProvider, boolean credentialWithheld,
                AgentRpcHost brokerHost = null) {
            this.model = model
            this.instruction = instruction
            this.goal = goal
            this.maxIterations = maxIterations
            this.requestTimeoutSeconds = requestTimeoutSeconds
            this.trace = trace
            this.agentName = agentName
            this.tools = tools
            this.toolSpecs = toolSpecs
            this.nativeToolNames = nativeToolNames
            this.skills = skills
            this.outputSchema = outputSchema
            this.apiKey = apiKey
            this.baseUrl = baseUrl
            this.apiProvider = apiProvider
            this.credentialWithheld = credentialWithheld
            this.brokerHost = brokerHost
        }

        AgentRunnerRequest createRequest(String prompt, String inputJson, ToolDispatcher dispatch, String workDir) {
            return new AgentRunnerRequest(
                model: model,
                instruction: instruction,
                prompt: prompt,
                maxIterations: maxIterations,
                outputSchema: outputSchema,
                inputJson: inputJson,
                requestTimeoutSeconds: requestTimeoutSeconds,
                goal: goal,
                agentName: agentName,
                trace: trace,
                tools: tools,
                toolSpecs: toolSpecs,
                // §5: the two halves of the selection are PARTITIONED across two fields, and the
                // partition is what keeps a runner-native tool out of the broker's allowlist
                nativeToolNames: nativeToolNames,
                dispatch: dispatch,
                skills: skills,
                // Reasoning models can reject an explicit temperature. Resume remains input-keyed.
                temperature: null,
                workDir: workDir,
                apiKey: apiKey,
                baseUrl: baseUrl,
                apiProvider: apiProvider,
                credentialWithheld: credentialWithheld,
                // each agent definition dials the address the guard resolved for IT, so a run that
                // mixes a local-docker agent with, say, a k8s one advertises each task its own
                brokerHost: brokerHost)
        }
    }

    AgentDef(BaseScript owner, String name, Map<String,Object> directives, List<AgentInput> inputs,
            List<AgentOutput> outputs, PromptDef prompt,
            List<ProcessFileInput> fileInputs = Collections.<ProcessFileInput>emptyList(),
            Map<String,ProcessFileOutput> fileOutputs = Collections.<String,ProcessFileOutput>emptyMap()) {
        this.owner = owner
        this.name = name
        this.simpleName = name
        this.baseName = name
        this.directives = directives
        this.inputs = inputs
        this.outputs = outputs
        this.prompt = prompt
        this.fileInputs = fileInputs
        this.fileOutputs = fileOutputs
    }

    @Override String getType() { TYPE }
    @Override String getName() { name }
    String getSimpleName() { simpleName }
    String getBaseName() { baseName }
    BaseScript getOwner() { owner }

    String getModel() { directives.get('model') as String }
    String getInstruction() { directives.get('instruction') as String }
    String getGoal() { directives.get('goal') as String }
    private List directiveList(String key) {
        final value = directives.get(key)
        return value == null ? [] : (value instanceof List ? (List) value : [value])
    }
    List getTools() { directiveList('tools') }
    List getSkills() { directiveList('skills') }
    /** The agent's declared `label` values, matched by `agent { withLabel: ... }` selectors. */
    List getLabels() { directiveList('label') }
    Integer getMaxIterations() { directives.get('maxIterations') as Integer }
    List<AgentInput> getInputs() { inputs }
    List<AgentOutput> getOutputs() { outputs }
    PromptDef getPrompt() { prompt }

    @Override
    ComponentDef cloneWithName(String name) {
        // register the alias and the workflow-scoped name so an `agent` scope `withName:` selector
        // targeting them is not reported as unmatched by Session#checkConfig (mirrors ProcessDef);
        // the agent name set is kept separate because a `process` selector never matches an agent
        ScriptMeta.addResolvedAgentName(name)
        def copy = (AgentDef) this.clone()
        copy.@name = name
        copy.@simpleName = ProcessDef.stripScope(name)
        // NOTE: the shallow clone deliberately preserves `owner` (so module-local skills and
        // relative tool paths keep resolving from the DEFINING module dir under an alias) and
        // `baseName` (so declared-name config selectors keep matching)
        return copy
    }

    /**
     * The raw {@code agent} config scope of the live {@link nextflow.Session}, i.e. both the
     * agent-only options and the task directives, plus any `withName:`/`withLabel:` selector
     * blocks. Empty when there is no active session or the scope is not defined.
     */
    private Map<String,Object> agentScope() {
        final session = Global.session as Session
        return (session?.config?.get('agent') as Map<String,Object>) ?: Collections.<String,Object>emptyMap()
    }

    /**
     * Resolve the AGENT-ONLY options of the {@code agent} config scope for THIS agent, applying
     * the same selector ladder as the task directives (see {@link AgentConfig#resolveOptions}).
     * Returns an empty (all-null) {@link AgentConfig} when there is no active session or the
     * scope is not defined, so callers can always read defaults.
     */
    protected AgentConfig agentConfig() {
        return new AgentConfig(AgentConfig.resolveOptions(agentScope(), getLabels() as List<String>, baseName, simpleName, name))
    }

    private DataflowReadChannel createSourceChannel(Object value) {
        if( value instanceof DataflowReadChannel || value instanceof DataflowBroadcast )
            return CH.getReadChannel(value)
        final result = CH.value()
        result.bind(value)
        return result
    }

    /**
     * Single lowering path (design §4.9, M-Tools): ALL agents — free-text, structured,
     * skills-only AND tool agents — lower to a real {@link TaskProcessor}/{@link TaskRun}
     * on the task path ({@link #runAsTask}) so every agent inherits work dir,
     * parallelism, progress table and lineage natively. A tool agent's shared
     * {@link ModuleToolBridge} gateway is created in {@link #buildAgentTask} (before ignition)
     * and invoked from the task-body thread; the legacy serial GPars operator path has
     * been removed.
     */
    @Override
    Object run(Object[] args0) {
        final args = ChannelOut.spread(args0)
        return runAsTask(args)
    }

    /**
     * Task path (design §4.1/§4.5, M-Tools): lower ANY agent — free-text, structured,
     * skills-only or tool — to a {@link ProcessConfigV2} + GROOVY ({@code exec}) {@link BodyDef}
     * and drive it through the unmodified {@link TaskProcessor} pipeline, mirroring
     * {@link ProcessDef#runV2}. The LLM call is the task body and stays in-JVM behind the
     * {@link AgentRunner} SPI; a tool agent's shared {@link ModuleToolBridge} is invoked from
     * that same task-body thread. Supports multiple structured inputs and multiple named
     * structured outputs; a queue input runs the agent per-item (map) while a value/singleton
     * input runs it once (fan-in).
     */
    private Object runAsTask(List args) {
        final built = buildAgentTaskWithBridge(args)
        final processor = built.processor
        final bridge = built.bridge
        // start the processor (progress table, lineage, events, work dir all for free)
        processor.run()
        final channels = outputChannels(processor)
        // A tool agent's request gateway consumes a persistent queue that stays alive until
        // poisoned; if never closed, session.await() hangs. Subscribe to the agent's
        // (first) output channel — captured pre-ignition so it sees the terminal poison —
        // and close the shared bridge once the agent's map completes (poisoning the request
        // queue so its operator terminates). On abort the session terminates all operators
        // directly (and afterStop still fires here), so close() is reached on both paths.
        if( bridge != null ) {
            final firstOut = channels.values().iterator().next()
            DataflowHelper.subscribeImpl(CH.getReadChannel(firstOut), [onComplete: { bridge.close() }])
        }
        return new ChannelOut(channels)
    }

    /**
     * The processor's declared output channels, keyed by param name, in declaration order.
     *
     * <p>The one dynamic hop in {@link #runAsTask}, isolated here so the rest of that method stays
     * statically checked: {@code TaskProcessor.getConfig()} is declared as the V1
     * {@link ProcessConfig}, which has no typed output-param accessor, and an agent's config is
     * always a {@link ProcessConfigV2}.
     */
    @CompileDynamic
    private static Map<String,DataflowWriteChannel> outputChannels(TaskProcessor processor) {
        final channels = new LinkedHashMap<String,DataflowWriteChannel>()
        for( final param : processor.getConfig().getOutputs().getParams() )
            channels.put(param.getName(), param.getChannel())
        return channels
    }

    /** Render the prompt closure against the task context (delegate-first). */
    private static String renderPrompt(PromptDef promptDef, Object ctx) {
        final Closure pc = (Closure) promptDef.closure.clone()
        pc.setDelegate(ctx)
        pc.setResolveStrategy(Closure.DELEGATE_FIRST)
        return pc.call()?.toString()
    }

    /** Bare single value for 1 input (byte-identical to the legacy path); {name:value} for N>1. */
    private static String buildInputJson(List<AgentInput> ins, Map ctx) {
        return ins.size() == 1
            ? toJson(ctx.get(ins[0].name))
            : toJson(ins.collectEntries { [(it.name): ctx.get(it.name)] })
    }

    /**
     * @param outputs the MODEL-ANSWERED outputs only (see the {@code modelOuts} partition in
     *        {@link #buildAgentTaskWithBridge}); an output with an explicit right-hand side is
     *        not part of the contract the model is given.
     */
    private static AgentOutputPlan resolveOutputPlan(String agentName, List<AgentOutput> outputs, List tools) {
        // an agent whose every output is a work-dir collection asks the model for nothing: its
        // observable result is the files it wrote, and its final text is discarded
        if( !outputs )
            return new AgentOutputPlan(AgentOutputMode.TEXT, null)
        if( outputs.size() > 1 )
            return new AgentOutputPlan(AgentOutputMode.WRAPPED, buildWrapperSchema(agentName, outputs))
        final output = outputs[0]
        if( Record.isAssignableFrom(output.type as Class) )
            return new AgentOutputPlan(AgentOutputMode.RECORD, RecordSchema.of(output.type as Class))
        if( tools )
            return new AgentOutputPlan(AgentOutputMode.SCALAR_CONTRACT, scalarOutputSchema(output))
        return new AgentOutputPlan(AgentOutputMode.TEXT, null)
    }

    @CompileDynamic
    private static Closure createCanonicalBody(PromptDef promptDef, List<AgentInput> inputs,
            ResolvedAgentSettings settings, AgentLaunchSpec launchSpec,
            AgentRunner runner, ModuleToolBridge bridge, boolean needsSandbox) {
        final String runnerName = runner.getName()
        return { ->
            final ctx = getDelegate()
            final String promptText = AgentDef.renderPrompt(promptDef, ctx)
            final String inputJson = AgentDef.buildInputJson(inputs, ctx)
            final TaskConfig taskConfig = ctx.get('task') as TaskConfig
            final Path taskWorkDir = taskConfig?.workDir as Path
            final DispatchContext dispatchContext = needsSandbox ? new DispatchContext(taskWorkDir) : null
            final ToolDispatcher contextualDispatch = bridge == null ? null : ({ String toolName, String argsJson ->
                if( dispatchContext != null )
                    ModuleToolBridge.setContext(dispatchContext)
                try {
                    return bridge.call(toolName, argsJson)
                }
                finally {
                    if( dispatchContext != null )
                        ModuleToolBridge.clearContext()
                }
            } as ToolDispatcher)
            final request = settings.createRequest(promptText, inputJson, contextualDispatch, '.')
            // A canonical launch command is made of paths that exist ONLY inside the runner image,
            // so the task always runs in a container (validated pre-ignition by
            // AgentLaunchConditions.requireContainerized) and the broker endpoint is therefore ALWAYS dialled from
            // outside the driver's network namespace -- hence remote=true unconditionally.
            // The image is re-checked HERE because `container` may be a lazy value (a closure or a
            // GString) that is truthy pre-ignition yet resolves to null per task.
            AgentLaunchConditions.requireTaskContainer(settings.agentName, runnerName, taskConfig?.getContainer())
            final AgentRpcRegistration registration = runner.register(request, true)
            return launchSpec.shellCommand([
                '--endpoint', registration.endpoint,
                '--invocation', registration.invocationId,
                *registration.transportArgs(),
                '--token', registration.token ])
        }
    }

    @CompileDynamic
    private static Closure createInJvmBody(PromptDef promptDef, List<AgentInput> inputs,
            List<AgentOutput> outputs, ResolvedAgentSettings settings, AgentRunner runner,
            ModuleToolBridge bridge, boolean needsSandbox, ProcessConfigV2 config,
            AgentOutputPlan outputPlan) {
        return { ->
            final ctx = getDelegate()
            final String promptText = AgentDef.renderPrompt(promptDef, ctx)
            final String inputJson = AgentDef.buildInputJson(inputs, ctx)
            final String workDir = (ctx.get('task')?.workDir as Path)?.toString()
            final request = settings.createRequest(promptText, inputJson, bridge as ToolDispatcher, workDir)
            // Clear stale snapshots before invoking a runner on this pooled task thread.
            AgentCallInfo.clear()
            if( needsSandbox )
                ModuleToolBridge.setContext(AgentDef.createSandboxContext(inputs, ctx))
            final Object result
            try {
                result = runner.run(request)
            }
            catch( Throwable error ) {
                // Consume a fatal tool abort's interrupt so it cannot leak to the next pooled task.
                Thread.interrupted()
                throw error
            }
            finally {
                if( needsSandbox )
                    ModuleToolBridge.clearContext()
            }
            final resolvedModel = AgentCallInfo.consumeResolvedModel()
            if( config.isCacheable() && resolvedModel != null )
                ctx.put('$agentResolvedModel', resolvedModel)
            outputPlan.bind(ctx, result, outputs)
            return null
        }
    }

    /**
     * The sandbox context for an in-JVM agent: the task work dir, plus the SOURCE of every input
     * that was staged into it.
     *
     * <p>Staging materializes an input as a symlink in the work dir ({@link
     * nextflow.executor.local.AgentTaskHandler}), while {@link nextflow.agent.SandboxGuard}
     * resolves symlinks before testing containment — so without this the agent would be handed a
     * name the {@code fs:} tools then refuse as "path outside sandbox", and `ls` would report it
     * as an opaque link. Whitelisting the staged source is the in-JVM analogue of the read-only
     * bind mount a canonical (containerized) agent gets for the same input: a readable entry is a
     * path-prefix test, so a file entry admits exactly that file and nothing beside it.
     */
    @CompileDynamic
    private static DispatchContext createSandboxContext(List<AgentInput> inputs, Object ctx) {
        final sandbox = new DispatchContext(ctx.get('task')?.workDir as Path)
        for( final inp : inputs )
            addStagedSources(sandbox, ctx.get(inp.name))
        return sandbox
    }

    /**
     * Walk an input value for staged paths. The shapes are exactly the ones
     * {@code TaskInputResolver.normalizeValue} produces: a bare {@link TaskPath}, a collection of
     * them, or a record (a {@code RecordMap} at runtime) holding them.
     */
    private static void addStagedSources(DispatchContext sandbox, Object value) {
        // toRealPath() on a TaskPath is a pure accessor for the store path -- no file I/O
        if( value instanceof TaskPath )
            sandbox.addReadablePath(value.toRealPath())
        else if( value instanceof Map )
            ((Map) value).values().each { addStagedSources(sandbox, it) }
        else if( value instanceof Collection )
            ((Collection) value).each { addStagedSources(sandbox, it) }
    }

    /**
     * Build (but do NOT run) the agent {@link TaskProcessor}: synthesize the
     * {@link ProcessConfigV2}, the canonical GROOVY {@link BodyDef}, create the tool bridge
     * (for a tool agent), wire the input/output channels and create the processor via
     * {@link ProcessDef#createTaskProcessorResolved}. Extracted from {@link #runAsTask} so resume tests
     * can inspect the built artifacts (canonical {@code BodyDef.source}, folded prompt
     * {@code valRefs}) and drive {@code checkCachedOutput} without igniting the dataflow network.
     */
    protected TaskProcessor buildAgentTask(List args) {
        return buildAgentTaskWithBridge(args).processor
    }

    private BuiltAgentTask buildAgentTaskWithBridge(List args) {
        requireInvocable(args)

        // -- resolve effective directives (identical to the legacy path)
        final agentConfig = agentConfig()
        final SelectedRunner selected = resolveRunner(agentConfig)
        final AgentLaunchSpec launchSpec = selected.launchSpec

        // -- build the V2 process config and resolve it from the `agent` scope BEFORE any
        //    skill/module resolution, so a misconfiguration fails fast (no git clone, no
        //    registry download, no gateway operator).
        final resolvedConfig = buildProcessConfig(selected)
        final config = resolvedConfig.config
        // the address THIS definition's task dials the driver's broker on; it travels on the
        // request so the broker advertises the address resolved for this agent rather than one
        // resolved for whichever definition happened to be built first
        final AgentRpcHost brokerHost = resolveLaunch(resolvedConfig, agentConfig, selected)

        final promptDef = this.prompt
        final List agentTools = this.tools
        // -- the output partition (two independent facts, not one rule):
        //    an output with an explicit right-hand side takes its value FROM that expression, so
        //    the model is neither asked for it (no schema entry) nor allowed to bind it; a
        //    `file(...)`/`files(...)` call inside that expression ADDITIONALLY registered an
        //    unstager, which is what makes it a work-dir collection
        final List<AgentOutput> modelOuts = outputs.findAll { it.value == null }
        final AgentOutputPlan outputPlan = resolveOutputPlan(name, modelOuts, agentTools)

        // -- capture read-only locals for the body closure (resolve lexically under
        //    DELEGATE_ONLY; do NOT reference `this.name`/`this.inputs`/etc. in the body)
        final String agentName = this.name
        final List<AgentInput> ins = this.inputs
        final List<AgentOutput> outs = this.outputs
        final Session session = Global.session as Session
        // resolve declared skills ONCE, pre-ignition (portable descriptors; no dataflow
        // coupling). Null when no skills are declared, so a tool-free/skill-free agent
        // carries a null `skills` on the request exactly as before.
        final List<SkillDescriptor> skillDescriptors = resolveSkills()

        final ResolvedTools resolvedTools = resolveTools(agentConfig, launchSpec, agentName, outputPlan, skillDescriptors)
        final settings = resolveSettings(agentConfig, agentName, agentTools, resolvedTools, skillDescriptors,
                outputPlan, brokerHost)

        declareParams(config, ins, outs, launchSpec, outputPlan)

        final body = createBody(promptDef, ins, modelOuts, settings, config, outputPlan, skillDescriptors,
                selected, resolvedTools, agentConfig)

        attachTaskInfo(config, selected, settings, promptDef, resolvedTools, skillDescriptors)

        wireChannels(config, args)

        applyExecutionPolicy(config, launchSpec)

        // build the processor (progress table, lineage, events, work dir all for free);
        // The config is fully resolved from the `agent` scope above. Use the canonical
        // processor/executor pipeline without applying the unrelated `process` scope.
        final processor = ProcessDef.createTaskProcessorResolved(session, owner, name, config, body)
        return new BuiltAgentTask(processor, resolvedTools.bridge)
    }

    /** Arity + zero-output guards (generalized; mirror ProcessDef.runV2). */
    private void requireInvocable(List args) {
        if( args.size() != inputs.size() )
            throw new ScriptRuntimeException("Agent `${name}` expects ${inputs.size()} input channel(s) but received ${args.size()}")
        if( !outputs )
            throw new ScriptRuntimeException("Agent `${name}` must declare exactly one output - zero outputs are not yet supported")
    }

    /**
     * Resolve the concrete runner before constructing the task. This makes the
     * selected implementation part of the cache identity even when the user
     * relies on the backwards-compatible single-runner default.
     */
    private static SelectedRunner resolveRunner(AgentConfig agentConfig) {
        final selectedRunner = AgentRunnerProvider.get(agentConfig.runner)
        return new SelectedRunner(selectedRunner, selectedRunner.getName(), selectedRunner.getLaunchSpec())
    }

    /**
     * Synthesize the agent's {@link ProcessConfigV2} and resolve it from the `agent` scope.
     * Called before any skill or module resolution, so a misconfiguration fails fast.
     */
    private ResolvedProcessConfig buildProcessConfig(SelectedRunner selected) {
        final config = new ProcessConfigV2(owner, name)
        final builder = new ProcessConfigBuilder(config, AgentDef.TYPE, AgentConfig.AGENT_ONLY_OPTIONS)
        // declared labels must reach the config BEFORE applyConfig so `withLabel:` can match;
        // route through the builder method so the value is a validated, de-duplicated ConfigList
        for( final lbl : getLabels() )
            builder.label(lbl as String)
        // Agent task placement and resources live in the `agent` scope and use the same
        // selector semantics as `process`. The `process` scope is never applied: process and
        // agent tasks are configured independently even though both reuse ProcessConfigV2.
        builder.applyConfig(agentScope(), baseName, simpleName, name)
        if( !config.containsKey('executor') )
            config.put('executor', AgentConfig.DEFAULT_EXECUTOR)
        // an agent is admitted without a cpu/capacity throttle, so absent a cap it fans out as wide
        // as its input channel -- bound concurrent LLM calls unless the user asked for more
        if( !config.containsKey('maxForks') )
            config.put('maxForks', AgentConfig.DEFAULT_MAX_FORKS)
        // legacy in-JVM runners cannot be offloaded; read the RESOLVED executor so a
        // selector-provided value is rejected loudly instead of silently downgraded below
        final resolvedExecutor = config.get('executor')?.toString()
        if( selected.launchSpec == null && resolvedExecutor != AgentConfig.DEFAULT_EXECUTOR )
            throw new ScriptRuntimeException("Agent runner `${selected.name}` does not support executor `${resolvedExecutor}`; only canonical launch-spec runners can be offloaded")
        return new ResolvedProcessConfig(config, resolvedExecutor)
    }

    /**
     * Admit a canonical launch and return the broker address this definition's task dials,
     * or {@code null} for an in-JVM runner, which dials nothing.
     */
    private AgentRpcHost resolveLaunch(ResolvedProcessConfig resolved, AgentConfig agentConfig, SelectedRunner selected) {
        if( selected.launchSpec == null )
            return null
        final config = resolved.config
        // a runner whose runtime lives IN an image knows which image that is -- it generates the
        // coordinate from its own VERSION -- so `agent.container` is optional for it. Put the value
        // in the config exactly as the user would have, so everything downstream (the guard below,
        // the per-task re-check in createCanonicalBody, the container fingerprint TaskHasher adds)
        // sees one kind of value and needs to know nothing about the default. This must stay AHEAD
        // of every `config.get('container')` read below.
        // `containsKey` and NOT truthiness: `agent.container = false` is the documented opt-out
        // (see AgentLaunchConditions.hasContainer) and must keep meaning "no container", i.e. must
        // keep failing below.
        if( !config.containsKey('container') ) {
            final runnerImage = selected.runner.getDefaultContainer()
            if( runnerImage )
                config.put('container', runnerImage)
        }
        // the RESOLVED executor, i.e. the same value buildProcessConfig admitted this runner on --
        // carried over from there rather than re-read, see ResolvedProcessConfig#executor
        final resolvedExecutor = resolved.executor
        // a canonical agent MUST be containerized on every executor: its launch command is built
        // from the absolute paths the runner's proxy and harness have INSIDE the runner image, which
        // exist nowhere else. Fail here rather than let the task fail with `No such file`.
        // `agent.runner` unset means the runner above was picked for the user, so a failure
        // naming it has to say where it came from
        // the task's `containerOptions` go IN, because the ladder reads them to decide whether
        // the container shares the driver's network namespace -- and they must be read BEFORE
        // AgentLaunchConditions.withDockerHostGateway appends to them, or that append would be inspected as if the
        // user had written it
        final launch = AgentLaunchConditions.requireCanonicalLaunch(name, selected.name, resolvedExecutor,
                config.get('container'), config.get('containerOptions'), agentConfig.rpc,
                Global.session as Session, !agentConfig.runner)
        // the task dials the broker back, so it may need a run option to be able to resolve the
        // driver's host at all -- see AgentLaunchConditions.withDockerHostGateway. Keyed off the RESOLVED address, not
        // the configured one, which is null for every inferred row
        final containerOptions = AgentLaunchConditions.withDockerHostGateway(config.get('containerOptions'),
                launch.containerEngine, launch.brokerHost?.host)
        if( containerOptions != null )
            config.put('containerOptions', containerOptions)
        return launch.brokerHost
    }

    /**
     * Expand the declared `tools` refs and lower them into the artifacts the request and the
     * cache key are built from: the bridge, the brokered descriptors, the runner-native names.
     */
    private ResolvedTools resolveTools(AgentConfig agentConfig, AgentLaunchSpec launchSpec, String agentName,
            AgentOutputPlan outputPlan, List<SkillDescriptor> skillDescriptors) {
        // -- expand the declared `tools` refs into the resolved selection ONCE: the brokered
        //    `nf:module_run:X` processes wired into the bridge below, and the runner-native
        //    `fs:`/`shell:` names the runner serves itself. A canonical launch-spec runner owns a
        //    container, which is what the `shell:` family requires.
        final ToolRefResolver.Selection toolSelection = resolveToolSelection(launchSpec != null)
        // -- §4: the model sees ONE flat namespace, so every selected tool must carry a legal wire
        //    name and no two sources may claim the same one. Checked HERE, before anything is
        //    wired: this is the only place that can see the whole namespace at once — the declared
        //    tools AND the names the runner injects for `skills` and structured output.
        checkWireNames(agentName, toolSelection, skillDescriptors != null && !skillDescriptors.isEmpty(),
                outputPlan.schema != null, launchSpec != null)
        // -- create the brokered tools' dataflow request gateway HERE (buildAgentTask runs before
        //    ignition), so the shared bridge is a captured final local the body closure can invoke
        //    from the task-body thread. Null for a tool-free/skills-only agent.
        final ModuleToolBridge bridge = new ModuleToolResolver(agentName, owner)
                .createToolBridge(toolSelection, launchSpec != null)
        if( bridge != null )
            bridge.setMaxInlineBytes(agentConfig.maxToolOutputInlineBytes())
        // §5: `toolSpecs` carries the BROKERED half only. The runner-native names travel beside it
        // on their own field, so a `fs:`/`shell:` name never enters the broker's allowlist and can
        // never be called back into the driver JVM.
        final List<ToolDescriptor> toolSpecs = bridge?.descriptors()
        final List<String> nativeToolNames = toolSelection?.nativeNames ?: null
        // the in-JVM fs: tools need a per-task sandbox context (the real task work dir). A
        // containerized runner serves them itself, so its bridge has none and needs no context.
        final boolean needsSandbox = bridge != null && bridge.filesystemEnabled
        return new ResolvedTools(toolSelection, bridge, toolSpecs, nativeToolNames, needsSandbox)
    }

    /** The effective, immutable values both request paths are built from. */
    private ResolvedAgentSettings resolveSettings(AgentConfig agentConfig, String agentName, List agentTools,
            ResolvedTools resolvedTools, List<SkillDescriptor> skillDescriptors, AgentOutputPlan outputPlan,
            AgentRpcHost brokerHost) {
        // the effective model, resolved BEFORE the settings because the endpoint and the
        // credential are scoped to its provider (see below)
        final String effectiveModel = this.model ?: agentConfig.model
        return new ResolvedAgentSettings(
            effectiveModel,
            this.instruction,
            this.goal,
            (this.maxIterations != null
                ? this.maxIterations
                : (agentConfig.maxIterations != null ? agentConfig.maxIterations : DEFAULT_MAX_ITERATIONS)) as int,
            (agentConfig.requestTimeout != null ? agentConfig.requestTimeout.seconds : DEFAULT_REQUEST_TIMEOUT_SECONDS) as int,
            agentConfig.traceEnabled(),
            agentName,
            agentTools,
            resolvedTools.toolSpecs,
            resolvedTools.nativeToolNames,
            skillDescriptors,
            outputPlan.schema,
            // endpoint and credential are resolved ONCE, in core, by the AgentConfig ladder
            // (config, then NXF_AGENT_*, then the provider's own variable) so no runner reads the
            // environment. Both are SCOPED to the model's API provider, and the provider tier of
            // the credential additionally requires an endpoint that provider owns: a runner
            // installs what it is given as the credential of the model's provider, ahead of
            // anything it could resolve itself -- see AgentConfig.apiKeyFor.
            agentConfig.apiKeyFor(effectiveModel),
            agentConfig.baseUrlFor(effectiveModel),
            agentConfig.apiProviderFor(effectiveModel),
            // "a key resolved but the gate refused to send it" is carried SEPARATELY from the null
            // credential above: only that case may be reported as an error by a runner with no
            // other source, and neither case may become the no-credential placeholder.
            agentConfig.credentialWithheldFor(effectiveModel),
            brokerHost)
    }

    /**
     * Complete the V2 process config (already resolved from the `agent` scope): one ProcessInput
     * per AgentInput, one ProcessOutput per AgentOutput with a synthetic lazy value closure that
     * reads the context slot the body writes.
     */
    @CompileDynamic
    private void declareParams(ProcessConfigV2 config, List<AgentInput> ins, List<AgentOutput> outs,
            AgentLaunchSpec launchSpec, AgentOutputPlan outputPlan) {
        for( final inp : ins )
            config.getInputs().addParam(inp.name, inp.type as Class, inp.optional)
        // replay the compiler-inferred stagers: this is what puts the declared Path inputs into
        // `task.inputFiles`, hence into the stage-in script AND the container bind mounts
        for( final f : fileInputs )
            config.getInputs().addFile(f)
        for( final out : outs ) {
            final String outName = out.name   // capture per-iteration (avoid loop-var capture)
            final Class outType = out.type as Class
            if( out.value != null ) {
                // the compiler's RHS closure -- for a file output `{ _file([:], '$path0') }`,
                // which resolves against TaskOutputResolver like any process output would
                config.getOutputs().addParam(outName, outType, out.value)
            }
            else if( launchSpec != null )
                config.getOutputs().addParam(outName, outType, {
                    outputPlan.decode(stdout(), outName, outType)
                })
            else
                config.getOutputs().addParam(outName, outType, { getProperty(outName) })
        }
        // replay the compiler-inferred unstagers, so `_file`/`_files` can resolve their key
        for( final entry : fileOutputs )
            config.getOutputs().addFile(entry.key, entry.value)
        // Collected outputs make `storeDir` newly valid for an agent (TaskProcessor.isInvalidStoreDir
        // tests the file-output set), and a task with a store dir has its outputs collected FROM it
        // -- which only works when something copies them there. That something is the wrapper's
        // unstage step, which an in-JVM agent does not have, so the task would fail with a missing
        // output on every run. Refuse the combination instead of failing per task.
        if( launchSpec == null && fileOutputs && config.get('storeDir') )
            throw new ScriptRuntimeException("Agent `${name}` cannot combine `storeDir` with a collected `path` output on an in-JVM runner - use `publishDir`, or select a containerized runner")
    }

    /** The task body closure, wrapped in the {@link BodyDef} that carries the cache identity. */
    private BodyDef createBody(PromptDef promptDef, List<AgentInput> ins, List<AgentOutput> modelOuts,
            ResolvedAgentSettings settings, ProcessConfigV2 config, AgentOutputPlan outputPlan,
            List<SkillDescriptor> skillDescriptors, SelectedRunner selected, ResolvedTools resolvedTools,
            AgentConfig agentConfig) {
        final AgentLaunchSpec launchSpec = selected.launchSpec
        final Closure bodyClosure = launchSpec != null
            ? createCanonicalBody(promptDef, ins, settings, launchSpec, selected.runner, resolvedTools.bridge, resolvedTools.needsSandbox)
            : createInJvmBody(promptDef, ins, modelOuts, settings, selected.runner, resolvedTools.bridge, resolvedTools.needsSandbox, config, outputPlan)

        // -- §7: a runner-native tool has no descriptor to hash, so the ONLY thing that can stand
        //    for its behaviour in the cache key is the runner that implements it. Resolved lazily
        //    (the plugin lookup is skipped entirely for an agent with no native tools), which also
        //    keeps the fingerprint of a brokered-only agent byte-identical to before.
        final List<String> nativeToolRefs = resolvedTools.selection?.nativeRefs
        final String runnerIdentity = nativeToolRefs ? runnerIdentity(selected.runner) : null

        // canonical BodyDef.source built from the EFFECTIVE resolved values (design §7.2/D2)
        // so a config-default model enters the cache key; fold the prompt's free-variable
        // refs into the synthetic BodyDef so params.* prompt-globals also enter the key (D3).
        return new BodyDef(
            bodyClosure,
            canonicalAgentSource(settings.model, settings.maxIterations, outputPlan.schema, skillDescriptors, selected.name,
                toolsFingerprint(resolvedTools.toolSpecs, resolvedTools.bridge?.toolSources(), nativeToolRefs, runnerIdentity),
                settings.baseUrl, agentConfig.apiProvider),
            launchSpec != null ? 'script' : 'exec',
            promptDef.valRefs)
    }

    /**
     * Attach the resolved agent identity to the config so an observer (lineage) can tell an
     * agent task from a process task and record what the agent actually was.
     */
    private void attachTaskInfo(ProcessConfigV2 config, SelectedRunner selected, ResolvedAgentSettings settings,
            PromptDef promptDef, ResolvedTools resolvedTools, List<SkillDescriptor> skillDescriptors) {
        // Immutable POJO, never a Map/Closure: TaskConfig is a LazyMap and would deep-copy/invoke
        // those on every read. Not one of the hashed directive names, so it stays out of the task hash.
        final agentInfo = new AgentTaskInfo(
            selected.name,
            settings.model,
            settings.instruction,
            settings.goal,
            promptDef?.source,
            settings.maxIterations,
            settings.outputSchema != null ? canonicalJson(settings.outputSchema) : null,
            // §7: lineage records the RESOLVED WIRE NAMES, so it must span both halves of the
            // partition -- a runner-native tool has no descriptor, and recording only `toolSpecs`
            // would silently drop every `fs:`/`shell:` tool from the record
            resolvedWireNames(resolvedTools.toolSpecs, resolvedTools.nativeToolNames),
            skillDescriptors ? skillDescriptors.collect { it.name } : null )
        config.put(AgentTaskInfo.CONFIG_KEY, agentInfo)
    }

    /** Invocation wiring (mirror ProcessDef.runV2). */
    private void wireChannels(ProcessConfigV2 config, List args) {
        final declaredInputs = config.getInputs().getParams()
        for( int i = 0; i < declaredInputs.size(); i++ )
            declaredInputs[i].setChannel(createSourceChannel(args[i]))
        final singleton = config.getInputs().isSingleton()
        for( final param : config.getOutputs().getParams() )
            param.setChannel(CH.create(singleton))
    }

    /**
     * Run agents on the dedicated `agent` executor: an agent body is an in-JVM orchestrator that
     * BLOCKS on its tool sub-tasks (and the LLM call), so it must not consume the compute
     * executor's cpu/capacity slot or a bounded worker thread — otherwise concurrent tool agents
     * deadlock the run. The agent executor uses an unthrottled monitor + unbounded thread pool;
     * the tool sub-tasks it dispatches still run throttled on the standard local executor.
     */
    private static void applyExecutionPolicy(ProcessConfigV2 config, AgentLaunchSpec launchSpec) {
        if( launchSpec == null )
            config.put('executor', 'agent')
    }

    /**
     * The complete set of wire names the model is offered, across BOTH halves of the §5
     * partition: the brokered descriptors plus the runner-native names. Order is descriptors
     * first, then natives in inventory order — the same order the model is given them in.
     *
     * @return the names, or {@code null} for a tool-free agent (so the lineage record keeps
     *         carrying {@code null} rather than an empty list, as it always has)
     */
    private static List<String> resolvedWireNames(List<ToolDescriptor> toolSpecs, List<String> nativeToolNames) {
        final List<String> names = new ArrayList<String>()
        if( toolSpecs )
            for( final ToolDescriptor d : toolSpecs )
                names.add(d.name)
        if( nativeToolNames )
            names.addAll(nativeToolNames)
        return names ?: null
    }

    /**
     * Canonical, deterministic identity string of the agent, set as the synthetic
     * {@link BodyDef#source} so {@link nextflow.processor.TaskHasher} folds the static agent
     * identity into the resume cache key with no hashing-code change (design §7.2/D2). Built
     * from the EFFECTIVE resolved values (so a config-default model invalidates the cache when
     * changed). The temperature line is a stable literal (`temperature=default`) because the
     * task path deliberately leaves temperature UNSET (M2 reconciliation).
     *
     * <p>Everything after {@code schema} defaults to {@code null}, and each of those five is
     * APPEND-ONLY: a null one contributes no line at all, so the string an agent that does not use
     * it produces is BYTE-FOR-BYTE what the narrower call always produced. That is what lets a
     * capability be added here without invalidating any existing agent's stored runs, and it is
     * pinned by {@code AgentDefTest}'s byte-for-byte equalities between adjacent arities.
     *
     * @param skills      folds the declared skills' identity in, so a changed {@code SKILL.md} (or
     *                    bundled resource) invalidates the cache instead of replaying a stale
     *                    generation (M-Skills correctness item). Non-empty appends a trailing
     *                    {@code skills=<fingerprint>} line. The fingerprint is order-independent:
     *                    the descriptors are sorted by name, then hashed over a stable list of
     *                    their identity fields (name/description/content/resources) via
     *                    {@link CacheHelper}
     * @param runner      the explicitly selected runner; prepends a leading {@code agentRunner=}
     *                    line, the one value written BEFORE the model
     * @param tools       the declared tools' fingerprint ({@link #toolsFingerprint}); appends a
     *                    trailing {@code tools=<fingerprint>} line
     * @param baseUrl     the RESOLVED endpoint: a different endpoint serves a different model under
     *                    the same id, so a replay must not be shared across endpoints (design D5).
     *                    Folding in the resolved value matches the rule that this string is built
     *                    from effective values -- so the same pipeline run with a different
     *                    {@code NXF_AGENT_BASE_URL} has a different key. The credential is
     *                    deliberately NOT folded in: it is not part of an agent's identity, and
     *                    hashing it would invalidate every entry on key rotation
     * @param apiProvider an EXPLICIT {@code agent.apiProvider}, which selects which environment
     *                    variables the endpoint and the credential come from (design D1/D6), so it
     *                    is part of how this agent was configured. Only the explicit value is folded
     *                    in -- an INFERRED provider is a pure function of {@code baseUrl}, which is
     *                    already in the key, so it adds nothing, and leaving it out means a later
     *                    addition to the inference table cannot silently invalidate anyone's stored
     *                    runs. Setting it -- even redundantly to the value that was already inferred
     *                    or taken from the model prefix -- appends a trailing
     *                    {@code apiProvider=<value>} line and invalidates that agent's entries once
     */
    protected String canonicalAgentSource(String model, int maxIter, Map schema, List<SkillDescriptor> skills = null, String runner = null, String tools = null, String baseUrl = null, String apiProvider = null) {
        final sb = new StringBuilder()
        if( runner )
            sb.append('agentRunner=').append(runner).append('\n')
        sb.append('agentModel=').append(model ?: '').append('\n')
        sb.append('temperature=default').append('\n')
        sb.append('instruction=').append(this.instruction ?: '').append('\n')
        sb.append('goal=').append(this.goal ?: '').append('\n')
        sb.append('maxIterations=').append(maxIter).append('\n')
        sb.append('prompt=').append(this.prompt?.source ?: '').append('\n')
        sb.append('outputSchema=').append(canonicalJson(schema))
        if( skills )
            sb.append('\n').append('skills=').append(skillsFingerprint(skills))
        if( tools )
            sb.append('\n').append('tools=').append(tools)
        if( baseUrl )
            sb.append('\n').append('baseUrl=').append(baseUrl)
        if( apiProvider )
            sb.append('\n').append('apiProvider=').append(apiProvider)
        return sb.toString()
    }

    /**
     * Deterministic, order-independent fingerprint of the declared tools' identity, so an agent
     * that can call a tool resumes only while that tool is unchanged. Sort the descriptors by
     * name, then hash each tool's identity: the descriptor the LLM actually sees (name,
     * description, input/output schema) plus the backing process' {@code BodyDef.source} — the
     * very string {@link nextflow.processor.TaskHasher} folds into that process' own task hash, so
     * the agent's key is exactly as sensitive to a tool edit as the tool task itself is.
     *
     * <p>Returns {@code null} for no tools, which keeps the canonical source byte-identical to the
     * tool-free form.
     *
     * <p>A {@code fs:}/{@code shell:} tool is served by the runner itself and therefore has NO
     * descriptor to hash — it would contribute <b>nothing</b> here, so an agent declaring
     * {@code fs:*} would share a cache key with the same agent declaring no tools at all, and a
     * runner upgrade that changes what {@code edit} or {@code grep} does would replay a stale
     * generation. What is folded in for those tools is the resolved <b>ref</b> plus the <b>runner
     * identity and version</b> ({@link #runnerIdentity}), not a schema: a native tool's behaviour
     * changes with the runner image, which is the granularity the plugin version already pins. It
     * also makes the key correctly runner-dependent — the same {@code fs:read} is a different
     * implementation on each runner (§7).
     *
     * <p>INVARIANT: with {@code nativeRefs} null or empty this returns exactly what it returned
     * before the runner-native pair existed (including {@code null} for an agent with no tools at
     * all), so no existing agent's cache key moves.
     *
     * @param nativeRefs the canonical refs of the runner-native tools, e.g. {@code [fs:read, shell:bash]};
     *                   order-independent, they are sorted before hashing
     * @param runnerId   the runner identity the native tools are served by, e.g. {@code pi@0.5.0}
     */
    static String toolsFingerprint(List<ToolDescriptor> tools, Map<String,String> sources, List<String> nativeRefs = null, String runnerId = null) {
        if( !tools && !nativeRefs )
            return null
        final List<Object> canonical = new ArrayList<>()
        for( final ToolDescriptor d : (tools ?: Collections.<ToolDescriptor>emptyList()).toSorted { it.name } )
            canonical.add([d.name, d.description, canonicalJson(d.inputSchema), canonicalJson(d.outputSchema), sources?.get(d.name)])
        if( nativeRefs )
            canonical.add(['runner-native', runnerId, new ArrayList<String>(new TreeSet<String>(nativeRefs))])
        return CacheHelper.hasher(canonical).hash().toString()
    }

    /**
     * The identity a runner-native tool's behaviour hangs off: the runner's stable name plus the
     * version of the plugin that supplies it — for {@code pi} the {@code nf-agent-pi} version,
     * which pins the runner image, which pins the SDK that implements the tool.
     *
     * <p>The version is recovered from the plugin that owns the runner class rather than from the
     * SPI, so no runner has to remember to report it. It is absent for a runner injected without a
     * plugin (the test seam of {@link AgentRunnerProvider}) and for an embedded distribution with
     * no plugin manager, in which case the name alone is used: a missing version must never make
     * the key non-deterministic within one installation.
     */
    static String runnerIdentity(AgentRunner runner) {
        if( runner == null )
            return null
        final String version = runnerVersion(runner)
        return version ? "${runner.getName()}@${version}".toString() : runner.getName()
    }

    private static String runnerVersion(AgentRunner runner) {
        try {
            final PluginWrapper wrapper = Plugins.getManager()?.whichPlugin(runner.getClass())
            return wrapper?.getDescriptor()?.getVersion()
        }
        catch( Exception e ) {
            log.debug("Unable to resolve the plugin version of agent runner `${runner.getName()}` - ${e.message}")
            return null
        }
    }

    /**
     * Deterministic, order-independent fingerprint of the declared skills' identity: sort the
     * descriptors by name, then hash a stable list of each skill's identity fields
     * (name, description, content, and each bundled resource's relativePath + content) via
     * {@link CacheHelper}. A changed {@code SKILL.md} body or bundled resource changes the hash.
     */
    protected static String skillsFingerprint(List<SkillDescriptor> skills) {
        final List<Object> canonical = new ArrayList<>()
        for( final SkillDescriptor d : skills.toSorted { it.name } ) {
            final List<Object> res = new ArrayList<>()
            for( final SkillResource r : (d.resources ?: Collections.<SkillResource>emptyList()) )
                res.add([r.relativePath, r.content])
            canonical.add([d.name, d.description, d.content, res])
        }
        return CacheHelper.hasher(canonical).hash().toString()
    }

    /**
     * Deterministic key-sorted JSON serialization of a (possibly nested) schema Map, used for
     * the output-schema fingerprint in {@link #canonicalAgentSource}. Recursively sorts Map
     * keys (TreeMap) so insertion order does not affect the fingerprint; List order is preserved
     * (it is semantically meaningful and built deterministically). Uses {@code JsonOutput.toJson}
     * rather than {@code Map.toString()} to avoid hash-order variance (design §7.2/D2).
     */
    protected static String canonicalJson(Object obj) {
        return JsonOutput.toJson(canonicalize(obj))
    }

    private static Object canonicalize(Object obj) {
        if( obj instanceof Map ) {
            final sorted = new TreeMap<String,Object>()
            for( final e : (obj as Map).entrySet() )
                // coerce a null key to '' so the TreeMap's natural ordering never NPEs
                // (defensive: RecordSchema.of/buildWrapperSchema only produce String keys)
                sorted.put(e.key != null ? e.key.toString() : '', canonicalize(e.value))
            return sorted
        }
        if( obj instanceof Collection )
            return (obj as Collection).collect { canonicalize(it) }
        return obj
    }

    /**
     * Synthesize the wrapper object schema for a multi-output agent (design §4.5/§5.3b):
     * one object whose {@code properties[out.name]} is the record schema (for record
     * outputs) or the scalar fragment (for supported scalar outputs), all names
     * {@code required}, {@code additionalProperties:false}. A top-level output whose
     * type is neither a record nor a supported scalar (e.g. {@code Path}, a top-level
     * collection) is rejected with a clear message.
     */
    static Map buildWrapperSchema(String agentName, List<AgentOutput> outs) {
        final props = new LinkedHashMap<String,Object>()
        final required = new ArrayList<String>()
        for( final o : outs ) {
            final Class t = o.type as Class
            final Map frag = (t != null && Record.isAssignableFrom(t))
                ? RecordSchema.of(t)
                : RecordSchema.scalarFragment(t)
            if( frag == null )
                throw new ScriptRuntimeException("Agent `${agentName}` output `${o.name}` has unsupported type ${t?.name} - supported: String, integer, number, boolean, or a record type")
            props.put(o.name, frag)
            required.add(o.name)
        }
        return ToolSchema.object(props, required)
    }

    /** Machine-readable wrapper for a single scalar output from a tool agent. */
    static Map scalarOutputSchema(AgentOutput out) {
        final Class type = out.type as Class
        final Map fragment = type != null && Path.isAssignableFrom(type)
            ? [type: 'string', description: 'Absolute path returned by the tool']
            : RecordSchema.scalarFragment(type)
        if( fragment == null )
            throw new ScriptRuntimeException("Agent output `${out.name}` has unsupported tool-result type ${type?.name}")
        final Map<String,Object> props = new LinkedHashMap<String,Object>()
        props.put(out.name, fragment)
        return ToolSchema.object(props, [out.name])
    }

    /**
     * Resolve the declared {@code skills} entries to portable {@link SkillDescriptor}s once,
     * pre-ignition. Each entry is either a remote GitHub reference (cloned + cached) or a local
     * skill name resolved under the {@code skills/} directory beside the script. Returns {@code null}
     * when no skills are declared; rejects duplicate skill names across all entries.
     */
    private List<SkillDescriptor> resolveSkills() {
        final declared = this.skills
        if( !declared )
            return null
        final session = Global.session as Session
        final meta = ScriptMeta.get(owner)
        final Path skillsRoot = ownerBaseDir(meta, session).resolve(SkillResolver.SKILLS_DIR)
        final List<SkillDescriptor> result = new ArrayList<>()
        final Set<String> seen = new HashSet<>()
        for( final entry : declared ) {
            final ref = entry?.toString()
            if( !ref )
                continue
            final List<SkillDescriptor> resolved = SkillResolver.isRemoteRef(ref)
                ? SkillResolver.loadRemote(skillsRoot, ref)
                : SkillResolver.loadLocal(skillsRoot, ref)
            for( final SkillDescriptor d : resolved ) {
                if( !seen.add(d.name) )
                    throw new ScriptRuntimeException("Agent `${name}`: duplicate skill name `${d.name}` - skills must have unique names")
                result.add(d)
            }
        }
        return result
    }

    /**
     * Expand the declared {@code tools} entries into the resolved selection, applying the
     * declaration grammar. Every entry is a namespaced ref — {@code nf:module_run[:PROCESS]},
     * {@code fs:<op>}, {@code shell:bash} — with no fallthrough: a bare process name, a module
     * path and a registry reference are all errors, and the capability they used to carry is
     * expressed by {@code include}ing the module and naming its process under {@code nf:module_run}.
     *
     * <p>The members of {@code nf:module_run} are the processes in scope for the owner script,
     * which is the same enumeration {@link nextflow.agent.ModuleToolResolver} then wires; the
     * {@code fs:} and {@code shell:} members are fixed by the release. Resolution is deliberately kept in
     * {@link nextflow.agent.ToolRefResolver}, a pure function of (refs, available members), so
     * the grammar is testable without a session.
     *
     * @param containerized whether the selected runner executes the agent inside its own
     *                      container. The {@code shell:} family needs that boundary: with an
     *                      in-JVM runner a shell tool would run LLM-authored commands on the
     *                      driver host, so the family is refused rather than served unsafely.
     * @return the resolved selection, or {@code null} when no tools are declared
     */
    private ToolRefResolver.Selection resolveToolSelection(boolean containerized) {
        final declared = this.tools
        if( !declared )
            return null
        final meta = ScriptMeta.get(owner)
        final procNames = meta != null ? meta.getProcessNames() : Collections.<String>emptySet()
        final shellUnavailable = containerized
                ? null
                : "the `${SHELL_FAMILY_RUNNER}` runner is required - an in-JVM runner would execute LLM-authored commands on the driver host with no container boundary (set `agent.runner = '${SHELL_FAMILY_RUNNER}'`)".toString()
        return ToolRefResolver.standard("Agent `${name}`".toString(), procNames, shellUnavailable).resolve(declared)
    }

    /**
     * Check the <b>wire</b> namespace — the flat list of tool names the LLM actually sees (§4) —
     * against its two requirements. Both are checked in this one pass because this is the only
     * point that can see the whole namespace at once: the resolved selection plus the names the
     * runner injects on the agent's behalf, which the {@link ModuleToolBridge} constructor cannot
     * know about.
     *
     * <ol>
     *   <li><b>Validate, never sanitize.</b> Every wire name must match the OpenAI function-name
     *       charset {@code [a-zA-Z0-9_-]{1,64}}. A Nextflow process name is a
     *       {@code JavaLetter JavaLetterOrDigit*} (ScriptLexer.g4), so {@code process my$proc},
     *       {@code process Σ_SORT} and names over 64 characters are all legal Nextflow and illegal
     *       on the wire — reachable here through a glob or a bare {@code nf:module_run}, since an
     *       explicitly-named ref could not carry those characters through the declaration grammar.
     *       Such a process is a hard error rather than a silent rename, because rewriting
     *       {@code my$proc} to {@code my_proc} would merge it with a process already called that,
     *       and the model would then call one and get the other.</li>
     *   <li><b>One wire name, one source.</b> Two <i>different</i> sources claiming the same name
     *       is a hard error naming both. Same-source duplicates are impossible by construction
     *       (the resolver returns a set, G9). The sources are sorted, never listed in declaration
     *       order, so the message is stable however the directive was written.</li>
     * </ol>
     *
     * <p>The injected names are not optional extras: the model cannot tell them apart from a tool,
     * so a process named {@code activate_skill} in a skills agent, or {@code final_answer} in a
     * structured-output agent on the canonical runner, is exactly the collision this rejects.
     *
     * @param agentName        label the errors are raised against
     * @param selection        the resolved tools, or {@code null} when none are declared
     * @param skillsDeclared   whether the agent declares {@code skills}, which makes the runner
     *                         inject {@code activate_skill}/{@code read_skill_resource}
     * @param structuredOutput whether the agent declares a structured output
     * @param containerized    whether the selected runner is a canonical launch-spec one; only that
     *                         one injects {@code final_answer} (the in-JVM runner decodes the
     *                         structured answer without a tool)
     */
    static void checkWireNames(String agentName, ToolRefResolver.Selection selection,
                               boolean skillsDeclared, boolean structuredOutput, boolean containerized) {
        // wire name -> the sources claiming it. A sorted set is what makes the message
        // deterministic under declaration order, and collapses a same-source duplicate.
        final Map<String,Set<String>> claims = new LinkedHashMap<String,Set<String>>()
        if( selection != null ) {
            for( final tool : selection.getTools() ) {
                checkWireName(agentName, tool)
                claim(claims, tool.name, "`${tool.ref}`".toString())
            }
        }
        if( skillsDeclared ) {
            claim(claims, SKILL_ACTIVATE_TOOL, 'the `skills` directive')
            claim(claims, SKILL_RESOURCE_TOOL, 'the `skills` directive')
        }
        if( structuredOutput && containerized )
            claim(claims, FINAL_ANSWER_TOOL, 'the agent output declaration')
        for( final entry : claims.entrySet() ) {
            if( entry.value.size() < 2 )
                continue
            throw new ScriptRuntimeException("Agent `${agentName}`: the tool name `${entry.key}` is claimed by ${entry.value.size()} different sources: ${entry.value.join(' and ')} - the model sees a single flat namespace, so rename the process or drop one of the refs")
        }
    }

    private static void claim(Map<String,Set<String>> claims, String name, String source) {
        // TreeSet: the sources of a collision are reported in their own sorted order, so the
        // message does not change when the directive entries are reordered
        Set<String> sources = claims.get(name)
        if( sources == null )
            claims.put(name, sources = new TreeSet<String>())
        sources.add(source)
    }

    /** Reject a resolved tool whose wire name the LLM API cannot carry; see {@link #checkWireNames}. */
    private static void checkWireName(String agentName, ToolRefResolver.ResolvedTool tool) {
        final String name = tool.name
        if( WIRE_NAME.matcher(name).matches() )
            return
        final List<String> reasons = new ArrayList<String>()
        final illegal = illegalWireChars(name)
        if( illegal )
            reasons.add("the illegal character(s) ${illegal.collect { "`${it}`" }.join(', ')}".toString())
        if( name.length() > WIRE_NAME_MAX )
            reasons.add("${name.length()} characters (the limit is ${WIRE_NAME_MAX})".toString())
        if( reasons.isEmpty() )
            reasons.add('a shape the wire namespace cannot carry')
        throw new ScriptRuntimeException("Agent `${agentName}`: tool `${tool.ref}` cannot be exposed to the model as `${name}` - the name has ${reasons.join(' and ')}, and a tool name must match `[a-zA-Z0-9_-]` with at most ${WIRE_NAME_MAX} characters. Rename it: the name is never rewritten automatically, because a rewrite could silently merge it with another tool")
    }

    /** The distinct characters of a name that the wire charset forbids, in order of first use. */
    private static List<String> illegalWireChars(String name) {
        final out = new LinkedHashSet<String>()
        for( int i = 0; i < name.length(); i++ ) {
            final String ch = name.substring(i, i + 1)
            if( !WIRE_NAME.matcher(ch).matches() )
                out.add(ch)
        }
        return new ArrayList<String>(out)
    }

    /**
     * The directory relative paths declared by the agent (e.g. its {@code skills/} root) are
     * resolved against: the owner script's directory when known, otherwise the session base dir,
     * otherwise the launch (current) directory.
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
     * Serialize the input record (a Map at runtime) to JSON, rendering any
     * {@link Path} value as an absolute, scheme-preserving string so the model receives
     * a portable representation (for example, {@code s3://bucket/key}).
     */
    protected static String toJson(Object item) {
        return JsonOutput.toJson(normalizeForJson(item))
    }

    private static Object normalizeForJson(Object value) {
        // A STAGED input is a TaskPath, whose alias is its identity inside the task dir. Render it
        // by toString() so the JSON agrees with the prompt's `${x}` interpolation of the same
        // input -- and so the model is given a name it can actually open in its runner. A TaskPath
        // is precisely and only the shape staging produces, so this check IS "was this staged";
        // it also cannot take toAbsolutePath(), which it throws on by design.
        if( value instanceof TaskPath )
            return value.toString()
        if( value instanceof Path )
            return FilesEx.toUriString(value.toAbsolutePath())
        if( value instanceof Map )
            return value.collectEntries { k, v -> [(k): normalizeForJson(v)] }
        if( value instanceof Collection )
            return value.collect { normalizeForJson(it) }
        return value
    }

}
