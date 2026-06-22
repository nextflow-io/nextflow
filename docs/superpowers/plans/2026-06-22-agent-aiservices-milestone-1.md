# Nextflow Agent Refactor — Milestone 1 (AiServices Engine Swap) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (- [ ]) syntax for tracking.

**Goal:** Replace the hand-rolled chat→tool→chat loop in `LangChainAgentRunner.runWithTools()` with a langchain4j `AiServices` proxy while preserving byte-identical observable behaviour (same tool names, same args JSON, same results fed back, same final text, same iteration-cap exception).

**Architecture:** The nf-agent plugin keeps two paths: the single-shot structured-output path (`runSingleShot`, unchanged) and the tool path (`runWithTools`, rewritten). The tool path now builds a `Map<ToolSpecification,ToolExecutor>` (each executor delegates to the core `ToolDispatcher` callback), seeds a `MessageWindowChatMemory` with only the `SystemMessage`, and drives the loop through an `AiServices`-generated proxy of a new minimal `AgentService` interface. Core (`modules/nextflow`) stays langchain4j-free; all langchain4j coupling remains plugin-local.

**Tech Stack:** Groovy 4.0.x (`@CompileStatic`), Spock, Gradle (`./gradlew`), langchain4j 1.16.3 (`langchain4j-open-ai` + `langchain4j` aggregator), pf4j plugin system.

## Global Constraints

- langchain4j must resolve to **exactly one** `langchain4j-core` at version **1.16.3** (bump `langchain4j-open-ai` to 1.16.3 AND add aggregator `dev.langchain4j:langchain4j:1.16.3`).
- Core module `modules/nextflow` must **never** import `dev.langchain4j.*` (boundary enforced by grep gate).
- The plugin must **never** touch `ProcessDef`, `Channel`, or `Path`.
- **Never** call `executeToolsConcurrently()` — tool calls stay sequential on the calling thread (AiServices default behaviour).
- All plugin production and test code is `@CompileStatic` Groovy where the existing files are (note: existing Spock specs are NOT `@CompileStatic`; keep them as-is — only `AgentService` and `LangChainAgentRunner` carry `@CompileStatic`).
- Tools-XOR-structured-output is enforced **upstream** in `run()` (not in this milestone); `runWithTools()` always passes a `null` schema to `ChatModelFactory.createModel`.
- Seed **only** the `SystemMessage` into memory; pass the full composed user text (prompt + inputJson) to `chat()`. The conversation must contain **exactly one** `UserMessage`.
- Iteration cap: pass `maxIterations + 1` to `maxSequentialToolsInvocations(int)` (cap of N permits N-1 round-trips), and wrap `agent.chat(...)` in try/catch to re-throw `IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})")`.
- All commits use `git commit -s` (DCO sign-off required).

---

## File Structure

| File | Status | Responsibility |
|------|--------|----------------|
| `plugins/nf-agent/build.gradle` | Modify (lines 53) | Bump `langchain4j-open-ai` to 1.16.3; add `langchain4j` aggregator 1.16.3 so `langchain4j-core` resolves to a single 1.16.3. |
| `plugins/nf-agent/src/main/nextflow/agent/AgentService.groovy` | **Create** | Minimal AiServices proxy interface: `String chat(String)`. No annotations on methods; named `AgentService` to avoid clashing with langchain4j-agentic `@Agent`. |
| `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy` | Modify (imports + `runWithTools` body, lines 104-138) | Replace manual loop with AiServices proxy; keep `run`, `runSingleShot`, `composeMessages`, `timeoutSeconds`, constants. |
| `plugins/nf-agent/src/test/nextflow/agent/LangChainAgentToolLoopTest.groovy` | Rewrite | Drive the real AiServices builder with a Spock-mocked `ChatModel`; assert on memory contents + dispatch invocation + null-schema factory call + cap exception. |
| `plugins/nf-agent/src/main/nextflow/agent/ChatModelFactory.groovy` | No change | `createModel(modelId, timeoutSeconds, schema)` consumed with `schema=null`. |
| `plugins/nf-agent/src/main/nextflow/agent/ModuleToolAdapter.groovy` | No change | `toToolSpecification(ToolDescriptor)` builds the `ToolSpecification`. |
| `modules/nextflow/src/main/groovy/nextflow/agent/{AgentRunnerRequest,ToolDispatcher,ToolDescriptor}.groovy` | No change | Core SPI; stays langchain4j-free. |
| `plugins/nf-agent/src/test/nextflow/agent/{LangChainAgentRunnerTest,AgentEndToEndTest,AgentToolEndToEndTest}.groovy` | No change | Single-shot and e2e tests; must still pass post-refactor. |

---

## Task 1: Bump langchain4j dependencies to a single 1.16.3 core

**Files:**
- Modify: `plugins/nf-agent/build.gradle` (lines 48-58, specifically line 53)
- Test: dependency-resolution check via `./gradlew :plugins:nf-agent:dependencies`

**Interfaces:**
- Consumes: nothing (build config).
- Produces: a classpath where `dev.langchain4j:langchain4j-core` resolves to exactly `1.16.3`, and AiServices classes (`dev.langchain4j.service.AiServices`) and memory classes (`dev.langchain4j.memory.chat.MessageWindowChatMemory`) are on the compile classpath via the aggregator.

### Steps

- [ ] **Read** the current dependencies block to confirm line 53 still reads `api 'dev.langchain4j:langchain4j-open-ai:1.15.1'`.

- [ ] **Apply the edit.** In `plugins/nf-agent/build.gradle`, replace:
  ```groovy
      api 'dev.langchain4j:langchain4j-open-ai:1.15.1'
  ```
  with:
  ```groovy
      api 'dev.langchain4j:langchain4j-open-ai:1.16.3'
      api 'dev.langchain4j:langchain4j:1.16.3'
  ```
  The full resulting dependencies block:
  ```groovy
  dependencies {
      compileOnly project(':nextflow')
      compileOnly 'org.slf4j:slf4j-api:2.0.17'
      compileOnly 'org.pf4j:pf4j:3.14.1'

      api 'dev.langchain4j:langchain4j-open-ai:1.16.3'
      api 'dev.langchain4j:langchain4j:1.16.3'

      testImplementation(testFixtures(project(":nextflow")))
      testImplementation "org.apache.groovy:groovy:4.0.31"
      testImplementation "org.apache.groovy:groovy-nio:4.0.31"
  }
  ```

- [ ] **Run the dependency resolution check.** Command:
  ```bash
  ./gradlew :plugins:nf-agent:dependencies --configuration runtimeClasspath | grep langchain4j-core
  ```
  **Expected output:** every `langchain4j-core` line resolves to `1.16.3` (lines like `+--- dev.langchain4j:langchain4j-core:1.15.1 -> 1.16.3` or `\--- dev.langchain4j:langchain4j-core:1.16.3`). There must be **no** resolved version other than `1.16.3`. Confirm with:
  ```bash
  ./gradlew :plugins:nf-agent:dependencies --configuration runtimeClasspath | grep 'langchain4j-core' | grep -oE '1\.[0-9]+\.[0-9]+$' | sort -u
  ```
  **Expected output:** a single line: `1.16.3`

- [ ] **Verify the aggregator is on the compile classpath** (AiServices and MessageWindowChatMemory classes are resolvable; both live in the aggregator jar, not in `langchain4j-core`). Command:
  ```bash
  ./gradlew :plugins:nf-agent:dependencies --configuration compileClasspath | grep -E 'dev.langchain4j:langchain4j:'
  ```
  **Expected output:** a line containing `dev.langchain4j:langchain4j:1.16.3` (no `-open-ai`, no `-core` suffix).

- [ ] **Commit.** Command:
  ```bash
  git add plugins/nf-agent/build.gradle && git commit -s -m "Bump nf-agent langchain4j to 1.16.3 single core

Bump langchain4j-open-ai to 1.16.3 and add the langchain4j aggregator
1.16.3 so AiServices and MessageWindowChatMemory classes resolve and
langchain4j-core pins to a single 1.16.3 version."
  ```

---

## Task 2: Create the minimal AiServices proxy interface `AgentService`

**Files:**
- Create: `plugins/nf-agent/src/main/nextflow/agent/AgentService.groovy`
- Test: `plugins/nf-agent/src/test/nextflow/agent/AgentServiceTest.groovy`

**Interfaces:**
- Consumes: nothing.
- Produces: `interface AgentService { String chat(String userMessage) }` — the type passed to `AiServices.builder(AgentService.class)`.

### Steps

- [ ] **Write a failing test** that asserts the interface shape (one method `String chat(String)`), so the file's existence and signature are pinned. Create `plugins/nf-agent/src/test/nextflow/agent/AgentServiceTest.groovy`:
  ```groovy
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

  import java.lang.reflect.Method

  import spock.lang.Specification

  class AgentServiceTest extends Specification {

      def 'should declare a single chat(String):String proxy method'() {
          when:
          Method[] methods = AgentService.getDeclaredMethods()

          then: 'exactly one declared method'
          methods.length == 1

          and: 'named chat, returning String, taking a single String argument'
          final m = methods[0]
          m.name == 'chat'
          m.returnType == String
          m.parameterTypes.toList() == [String]

          and: 'it is an interface with no annotations on the method (avoids @Agent clash)'
          AgentService.isInterface()
          m.annotations.length == 0
      }
  }
  ```

- [ ] **Run the failing test** (the interface does not exist yet). Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.AgentServiceTest"
  ```
  **Expected:** FAIL — compilation error `unable to resolve class AgentService` (or `BUILD FAILED` with `AgentServiceTest` not compiling).

- [ ] **Create the interface** `plugins/nf-agent/src/main/nextflow/agent/AgentService.groovy`:
  ```groovy
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

  import groovy.transform.CompileStatic

  /**
   * Minimal proxy interface implemented at runtime by langchain4j
   * {@code AiServices}. A single {@code chat} method drives one agent turn-set:
   * the proxy advertises the registered tools, dispatches tool-execution
   * requests, feeds results back to the model and returns the model's final
   * text answer.
   *
   * Named {@code AgentService} (not {@code Agent}) to avoid clashing with the
   * langchain4j-agentic {@code @Agent} annotation. The method carries no
   * annotations so the full user text is passed verbatim as the user message.
   *
   * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
   */
  @CompileStatic
  interface AgentService {

      String chat(String userMessage)

  }
  ```

- [ ] **Run the test** — now passes. Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.AgentServiceTest"
  ```
  **Expected:** `BUILD SUCCESSFUL`, `AgentServiceTest` 1 test passed.

- [ ] **Commit.** Command:
  ```bash
  git add plugins/nf-agent/src/main/nextflow/agent/AgentService.groovy plugins/nf-agent/src/test/nextflow/agent/AgentServiceTest.groovy && git commit -s -m "Add AgentService proxy interface for AiServices

Minimal interface 'String chat(String userMessage)' implemented at
runtime by langchain4j AiServices. Named AgentService to avoid clashing
with langchain4j-agentic @Agent; carries no method annotations."
  ```

---

## Task 3: Rewrite `runWithTools()` to use the AiServices proxy (test-first)

**Files:**
- Rewrite: `plugins/nf-agent/src/test/nextflow/agent/LangChainAgentToolLoopTest.groovy` (the behavioral test that drives this rewrite)
- Modify: `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy`
  - imports (lines 18-31): drop loop-only imports, add AiServices/memory/tool imports
  - `runWithTools(...)` body (lines 97-138): full replacement
  - KEEP unchanged: `run()` (60-68), `timeoutSeconds()` (74-76), `runSingleShot()` (82-95), constants (54-56)
  - ADD `composeUserText()` helper; refactor `composeMessages()` (145-155) to reuse it (byte-identical output)

**Interfaces:**
- Consumes:
  - `AgentRunnerRequest` (core): fields `model`, `instruction`, `prompt`, `maxIterations`, `inputJson`, `toolSpecs : List<ToolDescriptor>`, `dispatch : ToolDispatcher`, `requestTimeoutSeconds`.
  - `ToolDispatcher.call(String toolName, String argsJson) : String` (core SPI).
  - `ModuleToolAdapter.toToolSpecification(ToolDescriptor) : ToolSpecification`.
  - `ChatModelFactory.createModel(String modelId, int timeoutSeconds, JsonSchema schema) : ChatModel` (called with `schema = null`).
  - langchain4j 1.16.3: `ToolExecutor.execute(ToolExecutionRequest, Object) : String` (2-arg SAM); `AiServices.builder(Class).chatModel(ChatModel).tools(Map<ToolSpecification,ToolExecutor>).chatMemory(ChatMemory).systemMessageProvider(Function<Object,String>).maxSequentialToolsInvocations(int).build() : T`; `MessageWindowChatMemory.withMaxMessages(int) : ChatMemory`.
- Produces: `runWithTools(AgentRunnerRequest) : String` — the model's final text answer, or throws `IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})")`.

### Steps

- [ ] **Read** `LangChainAgentRunner.groovy` to confirm current imports (18-28) and the `runWithTools` method+javadoc (97-138), the `runSingleShot` use of `JsonSchema` (line 84), and `composeMessages` (145-155) match the grounding. Also **read** the current `LangChainAgentToolLoopTest.groovy` to confirm it matches the grounding (closure-coerced `ChatModel`, `capturedRequests`, the `capturedRequests[1].messages()` assertion, the 9-arg request ctor, the cap test).

- [ ] **Write the new behavioral test FIRST** (red). Overwrite `plugins/nf-agent/src/test/nextflow/agent/LangChainAgentToolLoopTest.groovy` with the full rewritten spec below. This spec asserts on the AiServices-seeded `MessageWindowChatMemory` contents (single `SystemMessage`, exactly one `UserMessage`, the assistant tool-request turn, the fed-back tool result), the null-schema factory call, and the cap exception:
  ```groovy
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
  import dev.langchain4j.data.message.AiMessage
  import dev.langchain4j.data.message.ChatMessage
  import dev.langchain4j.data.message.SystemMessage
  import dev.langchain4j.data.message.ToolExecutionResultMessage
  import dev.langchain4j.data.message.UserMessage
  import dev.langchain4j.model.chat.ChatModel
  import dev.langchain4j.model.chat.request.ChatRequest
  import dev.langchain4j.model.chat.request.json.JsonSchema
  import dev.langchain4j.model.chat.response.ChatResponse
  import spock.lang.Specification

  /**
   * Exercises {@link LangChainAgentRunner#runWithTools} which is now driven by a
   * langchain4j {@code AiServices} proxy. AiServices routes the LLM call through
   * {@code ChatModel.chat(ChatRequest)}, so a Groovy-closure-coerced ChatModel
   * mock overriding that overload remains the injection seam. AiServices builds
   * each ChatRequest from the shared (seeded) chat memory, so each request's
   * {@code messages()} snapshot reflects the accumulating memory contents.
   */
  class LangChainAgentToolLoopTest extends Specification {

      private static final Map GREET_INPUT_SCHEMA = [
          type: 'object',
          properties: [name: [type: 'string']],
          required: ['name'],
          additionalProperties: false,
      ]

      private static final Map GREET_OUTPUT_SCHEMA = [
          type: 'object',
          properties: [greeting: [type: 'string']],
          required: ['greeting'],
          additionalProperties: false,
      ]

      def 'should drive the AiServices tool loop: call dispatch then return the final text'() {
          given: 'a mock model that requests the greet tool first, then answers'
          List<ChatRequest> capturedRequests = []
          int calls = 0
          // AiServices invokes chat(ChatRequest); override only that overload.
          ChatModel model = [
              chat: { ChatRequest req ->
                  // snapshot the memory-derived messages as the proxy sees them
                  capturedRequests << req
                  calls++
                  if( calls == 1 ) {
                      // first turn: ask to run the greet tool
                      final ter = ToolExecutionRequest.builder()
                          .id('call-1')
                          .name('greet')
                          .arguments('{"name":"Ada"}')
                          .build()
                      return ChatResponse.builder().aiMessage(AiMessage.from([ter])).build()
                  }
                  // second turn: final plain-text answer (no tool requests)
                  return ChatResponse.builder().aiMessage(AiMessage.from('Final: greeted Ada')).build()
              }
          ] as ChatModel

          and: 'a stub dispatcher recording the call and returning a canned JSON result'
          List<List<String>> dispatched = []
          ToolDispatcher dispatch = { String name, String args ->
              dispatched << [name, args]
              return '{"greeting":"Hello Ada!"}'
          } as ToolDispatcher

          and: 'the runner wired to a factory that never forces a responseFormat'
          JsonSchema capturedSchema = null
          boolean factoryCalled = false
          def factory = Stub(ChatModelFactory) {
              createModel(_, _, _) >> { String id, int timeout, JsonSchema schema ->
                  factoryCalled = true
                  capturedSchema = schema
                  model
              }
          }
          def runner = new LangChainAgentRunner(modelFactory: factory)

          and: 'a request carrying the greet tool spec and the dispatch callback'
          def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, GREET_OUTPUT_SCHEMA)
          def req = new AgentRunnerRequest(
              'openai/gpt-5-mini',
              'Use the greet tool.',
              'greet Ada',
              5,
              [],
              null,
              null,
              [descriptor],
              dispatch)

          when:
          def answer = runner.run(req)

          then: 'the final text is returned'
          answer == 'Final: greeted Ada'

          and: 'the dispatcher was called exactly once with the right name and args'
          dispatched == [['greet', '{"name":"Ada"}']]

          and: 'the model was chatted twice (one round-trip per tool turn plus the final)'
          calls == 2

          and: 'tools were advertised on every request, including the greet spec'
          capturedRequests.size() == 2
          capturedRequests.every { it.toolSpecifications() && it.toolSpecifications()*.name().contains('greet') }

          and: 'no structured-output schema was forced (tools XOR responseFormat)'
          factoryCalled
          capturedSchema == null

          and: 'the FIRST request memory held exactly: system message + one user message (full prompt)'
          List<ChatMessage> first = capturedRequests[0].messages()
          first.count { it instanceof SystemMessage } == 1
          first.count { it instanceof UserMessage } == 1
          (first.find { it instanceof SystemMessage } as SystemMessage).text() == 'Use the greet tool.'
          (first.find { it instanceof UserMessage } as UserMessage).singleText() == 'greet Ada'

          and: 'the SECOND request memory grew with the assistant tool-request turn and the tool result'
          List<ChatMessage> second = capturedRequests[1].messages()
          // still exactly one user message — no prompt duplication
          second.count { it instanceof UserMessage } == 1
          // the assistant turn carrying the tool request is present
          second.find { it instanceof AiMessage && (it as AiMessage).hasToolExecutionRequests() } != null
          // the tool result message was fed back with the right name and JSON payload
          def resultMsg = second.find { it instanceof ToolExecutionResultMessage } as ToolExecutionResultMessage
          resultMsg != null
          resultMsg.toolName() == 'greet'
          resultMsg.text() == '{"greeting":"Hello Ada!"}'

          and: 'the original user prompt persisted unchanged across the loop'
          (second.find { it instanceof UserMessage } as UserMessage).singleText() == 'greet Ada'
      }

      def 'should compose user text with the input JSON as a single user message'() {
          given: 'a model that answers immediately on the first turn'
          List<ChatRequest> capturedRequests = []
          ChatModel model = [
              chat: { ChatRequest req ->
                  capturedRequests << req
                  return ChatResponse.builder().aiMessage(AiMessage.from('done')).build()
              }
          ] as ChatModel

          and:
          ToolDispatcher dispatch = { String name, String args -> '{}' } as ToolDispatcher
          def factory = Stub(ChatModelFactory) { createModel(_, _, _) >> model }
          def runner = new LangChainAgentRunner(modelFactory: factory)
          def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, null)
          def req = new AgentRunnerRequest(
              'openai/gpt-5-mini',
              'sys',
              'do it',
              5,
              [],
              null,
              '{"k":"v"}',
              [descriptor],
              dispatch)

          when:
          def answer = runner.run(req)

          then: 'the answer is returned'
          answer == 'done'

          and: 'exactly one user message carrying prompt + input JSON'
          List<ChatMessage> msgs = capturedRequests[0].messages()
          msgs.count { it instanceof UserMessage } == 1
          (msgs.find { it instanceof UserMessage } as UserMessage).singleText() == 'do it\n\nInput (JSON):\n{"k":"v"}'
      }

      def 'should throw IllegalStateException when the iteration cap is exceeded without a final answer'() {
          given: 'a model that always requests the tool, never answering'
          ChatModel model = [
              chat: { ChatRequest req ->
                  final ter = ToolExecutionRequest.builder()
                      .id('loop')
                      .name('greet')
                      .arguments('{"name":"Ada"}')
                      .build()
                  ChatResponse.builder().aiMessage(AiMessage.from([ter])).build()
              }
          ] as ChatModel

          and:
          int dispatchCalls = 0
          ToolDispatcher dispatch = { String name, String args -> dispatchCalls++; '{"greeting":"Hello"}' } as ToolDispatcher
          def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
              createModel(_, _, _) >> model
          })
          def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, null)
          def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'greet Ada', 2, [], null, null, [descriptor], dispatch)

          when:
          runner.run(req)

          then: 'the cap RuntimeException is surfaced as the historical IllegalStateException'
          def e = thrown(IllegalStateException)
          e.message == 'Agent exceeded the maximum number of tool-call iterations (2)'
      }
  }
  ```

  > **Cap-test note:** with `maxIterations = 2`, `runWithTools` passes `maxSequentialToolsInvocations(3)`. A cap of 3 permits 2 round-trips before AiServices throws its plain `RuntimeException`, which `runWithTools` re-throws as `IllegalStateException("Agent exceeded the maximum number of tool-call iterations (2)")`. The test asserts on the exception type and the exact historical message — it does NOT assert `dispatchCalls` count, because the AiServices round-trip-to-dispatch ratio is an internal langchain4j detail and is not part of the observable contract.

- [ ] **Run the rewritten spec against the OLD runner to confirm RED** (do this BEFORE applying the `runWithTools` rewrite below). Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.LangChainAgentToolLoopTest"
  ```
  **Expected:** `BUILD FAILED` — the memory-content assertions (single `SystemMessage`, exactly one `UserMessage` seeded into `MessageWindowChatMemory`, the fed-back `ToolExecutionResultMessage`) fail because the old hand-rolled loop builds the message list locally and never uses chat memory.

- [ ] **Replace the import block** in `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy`. Apply this exact keep/drop/add directive:
  - **Remove these imports:** `AiMessage`, `ToolExecutionResultMessage`, `ChatRequest`.
  - **Add these imports:** `java.util.function.Function`, `dev.langchain4j.memory.ChatMemory`, `dev.langchain4j.memory.chat.MessageWindowChatMemory`, `dev.langchain4j.service.AiServices`, `dev.langchain4j.service.tool.ToolExecutor`.
  - **Keep these imports:** `ToolExecutionRequest`, `ToolSpecification`, `ChatMessage`, `SystemMessage`, `UserMessage`, `ChatModel`, `ChatResponse`, **`JsonSchema`** (still used by `runSingleShot` at line 84).

  Replace the current import lines (18-28):
  ```groovy
  import dev.langchain4j.agent.tool.ToolExecutionRequest
  import dev.langchain4j.agent.tool.ToolSpecification
  import dev.langchain4j.data.message.AiMessage
  import dev.langchain4j.data.message.ChatMessage
  import dev.langchain4j.data.message.SystemMessage
  import dev.langchain4j.data.message.ToolExecutionResultMessage
  import dev.langchain4j.data.message.UserMessage
  import dev.langchain4j.model.chat.ChatModel
  import dev.langchain4j.model.chat.request.ChatRequest
  import dev.langchain4j.model.chat.request.json.JsonSchema
  import dev.langchain4j.model.chat.response.ChatResponse
  ```
  with the final import block:
  ```groovy
  import java.util.function.Function

  import dev.langchain4j.agent.tool.ToolExecutionRequest
  import dev.langchain4j.agent.tool.ToolSpecification
  import dev.langchain4j.data.message.ChatMessage
  import dev.langchain4j.data.message.SystemMessage
  import dev.langchain4j.data.message.UserMessage
  import dev.langchain4j.memory.ChatMemory
  import dev.langchain4j.memory.chat.MessageWindowChatMemory
  import dev.langchain4j.model.chat.ChatModel
  import dev.langchain4j.model.chat.request.json.JsonSchema
  import dev.langchain4j.model.chat.response.ChatResponse
  import dev.langchain4j.service.AiServices
  import dev.langchain4j.service.tool.ToolExecutor
  ```

- [ ] **Replace the `runWithTools` method body.** Replace the entire current method (lines 97-138, javadoc + body):
  ```groovy
      /**
       * Manual tool-call loop. The tool specifications are advertised on every
       * chat request; for each tool-execution request the model emits, the dispatch
       * callback runs the real module and its result is fed back to the model. The
       * loop ends when the model returns a final text answer (no tool requests) or
       * the iteration cap is reached.
       */
      private String runWithTools(AgentRunnerRequest request) {
          // Phase 2: tools XOR structured-output — do NOT force a responseFormat
          // when tools are in play (pass a null schema to the model factory).
          final model = modelFactory.createModel(request.model, timeoutSeconds(request), null)

          final List<ToolSpecification> specs = request.toolSpecs.collect { ModuleToolAdapter.toToolSpecification(it) }
          final List<ChatMessage> messages = composeMessages(request)
          final int maxIterations = request.maxIterations > 0 ? request.maxIterations : DEFAULT_MAX_ITERATIONS

          log.debug "Running agent model=${request.model}; tools=${specs*.name()}; maxIterations=${maxIterations}"

          AiMessage ai = null
          for( int i = 0; i < maxIterations; i++ ) {
              final ChatRequest req = ChatRequest.builder()
                  .messages(messages)
                  .toolSpecifications(specs)
                  .build()
              final ChatResponse response = model.chat(req)
              ai = response.aiMessage()
              if( !ai.hasToolExecutionRequests() ) {
                  // final answer
                  return ai.text()
              }
              // append the assistant turn (carrying the tool requests) ...
              messages.add(ai)
              // ... then run each requested tool and append its result
              for( ToolExecutionRequest ter : ai.toolExecutionRequests() ) {
                  log.debug "Agent tool call name=${ter.name()}; args=${ter.arguments()}"
                  final String result = request.dispatch.call(ter.name(), ter.arguments())
                  messages.add(ToolExecutionResultMessage.from(ter, result))
              }
          }

          throw new IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})")
      }
  ```
  with:
  ```groovy
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
          // Phase 2: tools XOR structured-output — do NOT force a responseFormat
          // when tools are in play (pass a null schema to the model factory).
          final ChatModel model = modelFactory.createModel(request.model, timeoutSeconds(request), null)

          final int maxIterations = request.maxIterations > 0 ? request.maxIterations : DEFAULT_MAX_ITERATIONS

          // Build one ToolSpecification->ToolExecutor entry per declared tool.
          // Each executor delegates to the core dispatch callback (which runs the
          // real module) using the langchain4j 2-arg executor signature.
          final Map<ToolSpecification,ToolExecutor> tools = new LinkedHashMap<>()
          for( final descriptor : request.toolSpecs ) {
              final ToolSpecification spec = ModuleToolAdapter.toToolSpecification(descriptor)
              final ToolExecutor executor = { ToolExecutionRequest ter, Object memoryId ->
                  log.debug "Agent tool call name=${ter.name()}; args=${ter.arguments()}"
                  return request.dispatch.call(ter.name(), ter.arguments())
              } as ToolExecutor
              tools.put(spec, executor)
          }

          // Seed memory with ONLY the system message; the user text is passed to chat().
          final ChatMemory memory = MessageWindowChatMemory.withMaxMessages(Integer.MAX_VALUE)
          if( request.instruction )
              memory.add(SystemMessage.from(request.instruction))

          // The full composed user text (prompt + input JSON) becomes the single user message.
          final String userText = composeUserText(request)

          log.debug "Running agent model=${request.model}; tools=${tools.keySet()*.name()}; maxIterations=${maxIterations}"

          // A cap of N permits only N-1 model round-trips, so pass maxIterations+1.
          final Function<Object,String> systemMessageProvider = { Object memoryId -> request.instruction } as Function<Object,String>
          final AgentService agent = AiServices.builder(AgentService)
              .chatModel(model)
              .tools(tools)
              .chatMemory(memory)
              .systemMessageProvider(systemMessageProvider)
              .maxSequentialToolsInvocations(maxIterations + 1)
              .build()

          try {
              return agent.chat(userText)
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
  ```

- [ ] **Add the `composeUserText` helper** and refactor `composeMessages` to reuse it. Insert `composeUserText` immediately before the existing `composeMessages`:
  ```groovy
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
  ```
  Then replace the body of `composeMessages` (lines 145-155):
  ```groovy
      private static List<ChatMessage> composeMessages(AgentRunnerRequest request) {
          String userText = request.prompt
          if( request.inputJson )
              userText += "\n\nInput (JSON):\n" + request.inputJson

          final List<ChatMessage> messages = new ArrayList<ChatMessage>()
          if( request.instruction )
              messages.add(SystemMessage.from(request.instruction))
          messages.add(UserMessage.from(userText))
          return messages
      }
  ```
  with:
  ```groovy
      private static List<ChatMessage> composeMessages(AgentRunnerRequest request) {
          final List<ChatMessage> messages = new ArrayList<ChatMessage>()
          if( request.instruction )
              messages.add(SystemMessage.from(request.instruction))
          messages.add(UserMessage.from(composeUserText(request)))
          return messages
      }
  ```

- [ ] **Compile the plugin.** Command:
  ```bash
  ./gradlew :plugins:nf-agent:compileGroovy
  ```
  **Expected:** `BUILD SUCCESSFUL`. (The 2-arg coercion `{ ToolExecutionRequest ter, Object memoryId -> ... } as ToolExecutor` and the `Function<Object,String>` coercion are valid under `@CompileStatic`: `ToolExecutor.execute(ToolExecutionRequest, Object):String` is the single abstract method and the interface is `@FunctionalInterface`. The `JsonSchema` import is retained, so `runSingleShot` still compiles.)

- [ ] **Run the rewritten tool-loop spec to confirm GREEN.** Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.LangChainAgentToolLoopTest"
  ```
  **Expected:** `BUILD SUCCESSFUL`, `3 tests completed, 0 failed`. (`UserMessage.singleText()` is present in 1.16.3, so the prompt assertions pass as written.)

- [ ] **Run the single-shot regression spec** (it does not touch `runWithTools`; `composeMessages` output is unchanged). Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.LangChainAgentRunnerTest"
  ```
  **Expected:** `BUILD SUCCESSFUL`, all `LangChainAgentRunnerTest` tests pass.

- [ ] **Commit.** Command:
  ```bash
  git add plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy plugins/nf-agent/src/test/nextflow/agent/LangChainAgentToolLoopTest.groovy && git commit -s -m "Drive agent tool loop via langchain4j AiServices

Replace the hand-rolled chat->tool->chat loop in runWithTools() with an
AiServices proxy: build a ToolSpecification->ToolExecutor map (each
executor delegates to the core ToolDispatcher), seed memory with only the
system message, pass the full user text to chat(). Cap is maxIterations+1
with the RuntimeException re-thrown as the historical IllegalStateException.
Tools stay sequential; structured-output single-shot path unchanged.
Rewrite the tool-loop test to assert on the seeded MessageWindowChatMemory
contents (single system message, one user message, fed-back tool result),
the null-schema factory call and the cap exception."
  ```

---

## Task 4: Milestone check-in — full module gate, boundary grep, gated E2E

**Files:** none modified — verification only.

**Interfaces:**
- Consumes: the complete nf-agent test suite, core source tree, optional `OPENAI_API_KEY`-gated e2e tests.
- Produces: green build, proven core langchain4j-free boundary, and (when credentials present) green e2e.

### Exit-criteria steps

- [ ] **Single-core dependency resolution check (re-confirm from Task 1).** Command:
  ```bash
  ./gradlew :plugins:nf-agent:dependencies --configuration runtimeClasspath | grep 'langchain4j-core' | grep -oE '1\.[0-9]+\.[0-9]+$' | sort -u
  ```
  **Expected output:** a single line `1.16.3`.

- [ ] **Boundary grep: core must never import langchain4j.** Command:
  ```bash
  grep -RIl 'dev\.langchain4j' modules/nextflow/src/main || echo "CLEAN: no langchain4j imports in core"
  ```
  **Expected output:** `CLEAN: no langchain4j imports in core`

- [ ] **Boundary grep: plugin must not touch ProcessDef/Channel/Path-as-workflow types.** Command:
  ```bash
  grep -RIn -E 'import nextflow\.(script\.ProcessDef|extension\.|Channel)|import java\.nio\.file\.Path' plugins/nf-agent/src/main || echo "CLEAN: plugin does not touch ProcessDef/Channel/Path"
  ```
  **Expected output:** `CLEAN: plugin does not touch ProcessDef/Channel/Path`

- [ ] **Guard: executeToolsConcurrently is never called.** Command:
  ```bash
  grep -RIn 'executeToolsConcurrently' plugins/nf-agent/src/main || echo "CLEAN: executeToolsConcurrently never called"
  ```
  **Expected output:** `CLEAN: executeToolsConcurrently never called`

- [ ] **Full nf-agent unit test suite (excludes credential-gated e2e via smoke).** Command:
  ```bash
  NXF_SMOKE=1 ./gradlew :plugins:nf-agent:test
  ```
  **Expected:** `BUILD SUCCESSFUL`. All of `AgentServiceTest`, `LangChainAgentRunnerTest`, `LangChainAgentToolLoopTest` pass; credential-dependent e2e specs (`AgentEndToEndTest`, `AgentToolEndToEndTest`) are skipped under smoke.

- [ ] **Gated E2E (only if `OPENAI_API_KEY` is set).** Command:
  ```bash
  if [ -n "$OPENAI_API_KEY" ]; then ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.AgentToolEndToEndTest" --tests "nextflow.agent.AgentEndToEndTest"; else echo "SKIP: OPENAI_API_KEY not set — e2e not run"; fi
  ```
  **Expected:** either `BUILD SUCCESSFUL` (both e2e specs pass against the live model, proving the AiServices path produces the same observable tool-call + final-text behaviour as before), or `SKIP: OPENAI_API_KEY not set — e2e not run`.

- [ ] **Confirm a clean worktree (verification only — no code commit here).** If the worktree is dirty from accidental edits, run `git status` and clean it before declaring Milestone 1 complete. Command:
  ```bash
  git status --porcelain || true
  ```
  **Expected:** empty output (`git status` shows `nothing to commit, working tree clean`).

### Milestone 1 exit criteria (all must hold)

1. `langchain4j-core` resolves to exactly one version, `1.16.3`.
2. `AgentService` interface exists with the single `String chat(String)` method, no method annotations.
3. `runWithTools()` uses `AiServices` with a `Map<ToolSpecification,ToolExecutor>`, seeds only the `SystemMessage`, passes full user text to `chat()`, and the conversation holds exactly one `UserMessage`.
4. The iteration cap is `maxIterations + 1`; the langchain4j `RuntimeException` on exceed is re-thrown as `IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})")`.
5. The single-shot structured-output path (`runSingleShot`, `composeMessages`) is behaviourally unchanged; `LangChainAgentRunnerTest` passes untouched.
6. Core (`modules/nextflow`) imports no `dev.langchain4j.*`; the plugin touches no `ProcessDef`/`Channel`/`Path`; `executeToolsConcurrently()` is never called.
7. The full nf-agent unit suite is green; e2e is green when credentials are available.
