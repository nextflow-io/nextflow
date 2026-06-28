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

import dev.langchain4j.agent.tool.ToolExecutionRequest
import dev.langchain4j.agent.tool.ToolSpecification
import dev.langchain4j.data.message.ChatMessage
import dev.langchain4j.data.message.SystemMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.invocation.InvocationContext
import dev.langchain4j.memory.ChatMemory
import dev.langchain4j.memory.chat.MessageWindowChatMemory
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.listener.ChatModelListener
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import dev.langchain4j.service.AiServices
import dev.langchain4j.service.tool.ToolExecutionResult
import dev.langchain4j.service.tool.ToolExecutor
import dev.langchain4j.service.tool.ToolProvider
import dev.langchain4j.service.tool.ToolProviderRequest
import dev.langchain4j.service.tool.ToolProviderResult
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.pf4j.Extension

/**
 * langchain4j-backed {@link AgentRunner}.
 *
 * When the request declares no tools, the runner does a single-shot chat
 * (optionally constraining the model to a structured-output JSON schema). When
 * the request declares one or more tools, the runner uses a langchain4j {@code AiServices}
 * proxy which registers the tool specifications with their executors. For each
 * tool-execution request the model emits, the executor delegates to the dispatch
 * callback (which runs the real module), and the proxy feeds the result back to
 * the model, looping until the model returns a final text answer or the iteration
 * cap is reached.
 *
 * For Phase 2 tools and structured-output are mutually exclusive: when tools are
 * in play no {@code responseFormat} is forced on the model.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Extension
@CompileStatic
class LangChainAgentRunner implements AgentRunner {

    private static final int DEFAULT_TIMEOUT_SECONDS = 120

    private static final int DEFAULT_MAX_ITERATIONS = 10

    ChatModelFactory modelFactory = new ChatModelFactory()

    @Override
    String run(AgentRunnerRequest request) {
        if( !request.model )
            throw new IllegalArgumentException("Agent `model` directive is required")

        return (request.toolSpecs || request.skills)
            ? runWithTools(request)
            : runSingleShot(request)
    }

    /**
     * The per-request LLM chat timeout (seconds): the configured value carried by
     * the request, or the built-in default when none was configured.
     */
    private static int timeoutSeconds(AgentRunnerRequest request) {
        return request.requestTimeoutSeconds > 0 ? request.requestTimeoutSeconds : DEFAULT_TIMEOUT_SECONDS
    }

    /**
     * Single-shot chat (no tool calls). When an output record type is declared
     * the model is constrained to a structured-output JSON schema.
     */
    private String runSingleShot(AgentRunnerRequest request) {
        // derive the structured-output schema (when an output record type is declared)
        final JsonSchema schema = request.outputSchema
            ? JsonSchemaMapper.toJsonSchema('Output', request.outputSchema)
            : null

        final model = modelFactory.createModel(request.model, timeoutSeconds(request), schema)

        final List<ChatMessage> messages = composeMessages(request)

        log.debug "Running agent model=${request.model}; messages=${messages.size()}; structured=${schema != null}"
        final ChatResponse response = model.chat(messages)
        return response.aiMessage().text()
    }

    /**
     * Tool-call loop driven by a langchain4j {@code AiServices} proxy. The tool
     * specifications are registered with their executors; for each tool-execution
     * request the model emits, the executor delegates to the dispatch callback
     * which runs the real module, and the proxy feeds the result back to the
     * model — looping until the model returns a final text answer or the
     * iteration cap is reached.
     *
     * Tool calls are dispatched sequentially on the calling thread (the
     * AiServices default); {@code executeToolsConcurrently} is never enabled.
     *
     * Prompt composition: only the optional {@link SystemMessage} is seeded into
     * the chat memory; the full composed user text (prompt plus input JSON) is
     * passed to {@code chat(...)} so the conversation holds exactly one
     * {@link UserMessage}.
     */
    private String runWithTools(AgentRunnerRequest request) {
        // optional execution tracer: a ChatModelListener that logs the turns and
        // model reasoning, paired with per-tool input/output lines below
        final AgentTrace trace = request.trace ? new AgentTrace(request.agentName) : null
        final List<ChatModelListener> listeners = trace != null ? Collections.<ChatModelListener>singletonList(trace) : null

        // Phase 2: tools XOR structured-output — do NOT force a responseFormat
        // when tools are in play (pass a null schema to the model factory). The model is
        // built with listeners only when tracing, keeping the default path on the 3-arg call.
        final ChatModel model = listeners
            ? modelFactory.createModel(request.model, timeoutSeconds(request), null, listeners)
            : modelFactory.createModel(request.model, timeoutSeconds(request), null)

        final int maxIterations = request.maxIterations > 0 ? request.maxIterations : DEFAULT_MAX_ITERATIONS

        // Build one ToolSpecification->ToolExecutor entry per declared tool.
        // Each executor delegates to the core dispatch callback (which runs the
        // real module) using the langchain4j 2-arg executor signature. May be empty
        // for a skills-only agent.
        final Map<ToolSpecification,ToolExecutor> tools = new LinkedHashMap<>()
        if( request.toolSpecs ) {
            for( final descriptor : request.toolSpecs ) {
                final ToolSpecification spec = ModuleToolAdapter.toToolSpecification(descriptor)
                final ToolExecutor executor = { ToolExecutionRequest ter, Object memoryId ->
                    log.debug "Agent tool call name=${ter.name()}; args=${ter.arguments()}"
                    final String result = request.dispatch.call(ter.name(), ter.arguments())
                    if( trace != null )
                        trace.tool(ter.name(), ter.arguments(), result)
                    return result
                } as ToolExecutor
                tools.put(spec, executor)
            }
        }

        // Map declared skills onto a langchain4j Skills container (Tool Mode); its tool provider
        // adds the activate_skill/read_skill_resource tools alongside the static module tools.
        final dev.langchain4j.skills.Skills skills = request.skills
            ? SkillAdapter.toSkills(request.skills as List<SkillDescriptor>)
            : null

        // Seed memory with ONLY the system message (instruction + optional goal); for a skills run
        // append the available-skills catalog so the model knows which skills it can activate (the
        // activate_skill `skill_name` parameter is free text, not an enum). The user text is passed
        // to chat().
        final ChatMemory memory = MessageWindowChatMemory.withMaxMessages(Integer.MAX_VALUE)
        String systemMessage = composeSystemMessage(request)
        if( skills != null ) {
            final String catalog = skills.formatAvailableSkills()
            systemMessage = systemMessage ? (systemMessage + '\n\n' + catalog) : catalog
        }
        if( systemMessage )
            memory.add(SystemMessage.from(systemMessage))

        // The full composed user text (prompt + input JSON) becomes the single user message.
        final String userText = composeUserText(request)

        // The skills tool provider (activate_skill / read_skill_resource). When tracing, wrap it so
        // each skill-tool invocation is reported like a module tool; also enumerate its tool names
        // (best-effort) so they appear in the run header alongside the static module tools.
        final List<String> skillToolNames = new ArrayList<>()
        ToolProvider skillProvider = null
        if( skills != null ) {
            skillProvider = skills.toolProvider()
            if( trace != null ) {
                skillToolNames.addAll(enumerateSkillTools(skillProvider, userText))
                skillProvider = tracingSkillProvider(skillProvider, trace)
            }
        }

        log.debug "Running agent model=${request.model}; tools=${tools.keySet()*.name()}; skills=${skillToolNames}; maxIterations=${maxIterations}"
        if( trace != null )
            trace.begin(request.model, (tools.keySet()*.name() as List<String>) + skillToolNames)

        // A cap of N permits only N-1 model round-trips, so pass maxIterations+1.
        // The system message is seeded directly into the chat memory above;
        // no systemMessageProvider is needed (that would double-add it).
        // In langchain4j 1.17 static .tools(map) and .toolProvider(...) coexist and merge, so the
        // skills provider is additive; .tools(...) is skipped when there are no static module tools.
        final AiServices<AgentService> builder = AiServices.builder(AgentService)
            .chatModel(model)
            .chatMemory(memory)
            .maxSequentialToolsInvocations(maxIterations + 1)
        if( !tools.isEmpty() )
            builder.tools(tools)
        if( skillProvider != null )
            builder.toolProvider(skillProvider)
        final AgentService agent = builder.build()

        try {
            final String answer = agent.chat(userText)
            if( trace != null )
                trace.end(answer)
            return answer
        }
        catch( IllegalStateException e ) {
            // Let any genuine IllegalStateException (e.g. from ChatModelFactory or
            // our own code) propagate unchanged.
            throw e
        }
        catch( RuntimeException e ) {
            // langchain4j 1.16.3 throws a plain RuntimeException on cap exceed
            // (dev.langchain4j.internal.Exceptions.runtime, message "Something is
            // wrong, exceeded %s tool calling round trips ..."). Re-throw with the
            // historical IllegalStateException message/shape.
            throw new IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})", e)
        }
    }

    /**
     * Wrap a skills {@link ToolProvider} so each skill-tool invocation is reported to the
     * {@link AgentTrace} (the {@code → activate_skill(...)} line), the same way module-tool executors
     * report theirs. langchain4j-skills' executors carry their real logic in {@code executeWithContext}
     * (it reads activation state from the chat-message history via the {@link InvocationContext}), and
     * AiServices invokes {@code executeWithContext} — so the wrapper MUST override and delegate that
     * method (not just {@code execute}), or skill activation would break. {@code isDynamic()} is also
     * delegated so the provider's invocation cadence is preserved.
     */
    @PackageScope
    static ToolProvider tracingSkillProvider(ToolProvider delegate, AgentTrace trace) {
        return new ToolProvider() {
            @Override
            boolean isDynamic() { return delegate.isDynamic() }

            @Override
            ToolProviderResult provideTools(ToolProviderRequest request) {
                final ToolProviderResult res = delegate.provideTools(request)
                final Map<ToolSpecification,ToolExecutor> wrapped = new LinkedHashMap<>()
                for( final Map.Entry<ToolSpecification,ToolExecutor> entry : res.tools().entrySet() ) {
                    final ToolExecutor inner = entry.value
                    wrapped.put(entry.key, new ToolExecutor() {
                        @Override
                        String execute(ToolExecutionRequest req, Object memoryId) {
                            return inner.execute(req, memoryId)
                        }
                        @Override
                        ToolExecutionResult executeWithContext(ToolExecutionRequest req, InvocationContext ctx) {
                            final ToolExecutionResult result = inner.executeWithContext(req, ctx)
                            trace.tool(req.name(), req.arguments(), result?.resultText())
                            return result
                        }
                    })
                }
                return new ToolProviderResult(wrapped)
            }
        }
    }

    /**
     * Best-effort enumeration of a skills provider's tool names (e.g. {@code activate_skill},
     * {@code read_skill_resource}) for the trace header. Failures degrade to an empty list — the
     * per-call itemization above does not depend on this.
     */
    private static List<String> enumerateSkillTools(ToolProvider provider, String userText) {
        try {
            final ToolProviderRequest req = new ToolProviderRequest(null, UserMessage.from(userText))
            return new ArrayList<String>(provider.provideTools(req).tools().keySet()*.name())
        }
        catch( Exception e ) {
            log.debug("Agent trace: could not enumerate skill tool names: ${e.message}")
            return Collections.<String>emptyList()
        }
    }

    /**
     * Compose the system message: the {@code instruction} (role/persona),
     * followed by an optional {@code goal} section that steers the multi-turn
     * loop toward an objective. Returns {@code null} when neither is set, so the
     * caller seeds no {@link SystemMessage}.
     */
    private static String composeSystemMessage(AgentRunnerRequest request) {
        final StringBuilder sb = new StringBuilder()
        if( request.instruction )
            sb.append(request.instruction)
        if( request.goal ) {
            if( sb.length() > 0 )
                sb.append('\n\n')
            sb.append('Goal:\n').append(request.goal)
              .append('\nYou are done when the goal is met; produce your final answer as plain text.')
        }
        // when tracing a TOOL run, ask the model to narrate its reasoning so the trace can surface
        // it: OpenAI emits no reasoning on the wire for tool-call turns, so a one-line rationale in
        // the message content is the only way to make the "thinking" visible. Gated on toolSpecs so
        // it never touches the single-shot / structured-output path; trace-only, so normal runs are
        // unaffected. It nudges the model to think out loud, which generally aids tool selection.
        if( request.trace && (request.toolSpecs || request.skills) ) {
            if( sb.length() > 0 )
                sb.append('\n\n')
            sb.append('Before each tool call, briefly state your reasoning for the call in one short sentence.')
        }
        return sb.length() > 0 ? sb.toString() : null
    }

    /**
     * Compose the user text: the rendered prompt, plus the input record
     * serialized as JSON when present. This is the single user message passed
     * to the AiServices proxy.
     */
    private static String composeUserText(AgentRunnerRequest request) {
        String userText = request.prompt
        if( request.inputJson )
            userText += "\n\nInput (JSON):\n" + request.inputJson
        return userText
    }

    /**
     * Compose the chat messages: an optional system instruction followed by the
     * user message (the rendered prompt, plus the input record serialized as JSON
     * when present).
     */
    private static List<ChatMessage> composeMessages(AgentRunnerRequest request) {
        final List<ChatMessage> messages = new ArrayList<ChatMessage>()
        final String systemMessage = composeSystemMessage(request)
        if( systemMessage )
            messages.add(SystemMessage.from(systemMessage))
        messages.add(UserMessage.from(composeUserText(request)))
        return messages
    }
}
