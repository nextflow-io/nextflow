# Nextflow Agent Refactor — Milestone 2 (Goal-Directed Loop) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one optional `goal` directive that states a high-level objective, folded into the agent's system message to steer the existing AiServices multi-turn loop — without any new runtime primitive, grammar change, or completion-evaluator.

**Architecture:** `goal` is a `String`-valued agent directive shaped exactly like `instruction`. It is recognised at resolve time (`AgentDsl`), captured at runtime (`AgentBuilder.DIRECTIVES` → `AgentDef.directives` map), carried on `AgentRunnerRequest`, and combined with `instruction` into the single system message by a new `composeSystemMessage()` helper used by both the tool path and the single-shot path. `maxIterations` remains the hard brake; completion stays "model stops calling tools".

**Tech Stack:** Groovy 4.0.x (`@CompileStatic`/`@CompileDynamic`), Java (nf-lang ANTLR/AST/DSL), Spock, Gradle (`./gradlew`), langchain4j 1.16.3 (from M1).

## Global Constraints

- `goal` is **advisory**: it never raises `maxIterations` and adds no new loop/exit primitive. `maxIterations` stays the hard cap.
- `goal` is a `String` directive identical in shape to `instruction`; **no grammar/lexer/`AgentNode`/visitor/lowering change** is permitted or needed.
- Core module `modules/nextflow` must **never** import `dev.langchain4j.*`; the `nf-agent` plugin must **never** touch `ProcessDef`/`Channel`/`Path`.
- `AgentRunnerRequest` is `@Canonical`: **append** `goal` as the LAST field (after `requestTimeoutSeconds`); never insert it mid-list. Update the positional-order javadoc in lockstep.
- `AgentDef.run` is `@CompileDynamic` — construct `AgentRunnerRequest` with the Groovy **named-args** map form so future field additions (M3) are non-breaking.
- The system message must be **exactly one** `SystemMessage` (or none): `composeSystemMessage` returns a single combined string, or `null` when both `instruction` and `goal` are absent.
- All plugin/core production code keeps its existing `@CompileStatic`/`@CompileDynamic` annotations. All commits use `git commit -s` (DCO).

---

## File Structure

| File | Status | Responsibility |
|------|--------|----------------|
| `modules/nf-lang/src/main/java/nextflow/script/dsl/AgentDsl.java` | Modify (after the `instruction` method, ~line 38) | Declare the `goal(String)` directive so `goal '...'` resolves at compile time. |
| `modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy` | Modify (line 35) | Add `'goal'` to the `DIRECTIVES` allowlist so the runtime delegate captures it. |
| `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` | Modify (getter ~line 96; `run()` ~lines 170 & 195) | Add `getGoal()`; capture `agentGoal`; pass it via named-args construction. |
| `modules/nextflow/src/main/groovy/nextflow/agent/AgentRunnerRequest.groovy` | Modify (field list + javadoc) | Append `String goal` as the last `@Canonical` field. |
| `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy` | Modify (`runWithTools` seed ~lines 135-138; `composeMessages` ~lines 189-194; add `composeSystemMessage`) | Fold `instruction` + `goal` into the single system message in both paths. |
| `docs/agent.mdx` | Modify (Directives table + a short paragraph) | Document `goal` as an advisory directive. |
| Tests (see tasks) | Modify/Create | Resolution+capture, SPI field order, runner composition, gate. |

---

## Task 1: Recognise and capture the `goal` directive

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/script/dsl/AgentDsl.java`
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy:35`
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` (add getter after line 96)
- Test: `modules/nextflow/src/test/groovy/nextflow/script/AgentDefTest.groovy` (getGoal), and `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy` (a script with `goal` loads and exposes it)

**Interfaces:**
- Consumes: the existing directive machinery (`AgentBuilder.methodMissing` → `directives` map; `AgentDef.directives`).
- Produces: `AgentDef.getGoal() : String` returning the declared goal or `null`.

### Steps

- [ ] **Write the failing test.** Add to `AgentDefTest.groovy` a spec that builds an `AgentDef` whose directives include a goal and asserts `getGoal()`. Follow the existing `AgentDefTest` construction pattern; the direct-construction form is:
  ```groovy
  def 'should expose the goal directive'() {
      given:
      def directives = [model: 'openai/gpt-5-mini', instruction: 'be careful', goal: 'assemble then QC'] as Map<String,Object>
      def prompt = new PromptDef({ -> 'hi' }, 'hi')
      def agent = new AgentDef(Mock(BaseScript), 'a', directives, [], [], prompt)

      expect:
      agent.goal == 'assemble then QC'
      agent.instruction == 'be careful'
  }

  def 'should return null goal when not declared'() {
      given:
      def agent = new AgentDef(Mock(BaseScript), 'a', [model: 'openai/gpt-5-mini'] as Map<String,Object>, [], [],
          new PromptDef({ -> 'hi' }, 'hi'))
      expect:
      agent.goal == null
  }
  ```
  (If the existing `AgentDefTest` constructs `AgentDef`/`PromptDef` differently, match that file's existing pattern exactly — read it first.)

- [ ] **Run it — expect FAIL.** Command:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.AgentDefTest"
  ```
  Expected: FAIL — `agent.goal` does not resolve (`getGoal()` does not exist) → compile/MissingMethod failure.

- [ ] **Add the `goal` directive declaration** in `AgentDsl.java`, immediately after the `instruction(String)` method (after line 38):
  ```java
          @Description("""
              The `goal` directive states a high-level objective that steers the agent's multi-turn loop. It is advisory: the model is encouraged to keep working until the goal is met, while `maxIterations` remains the hard cap.
          """)
          void goal(String value);
  ```

- [ ] **Add `'goal'` to the runtime directive allowlist** in `AgentBuilder.groovy:35`. Replace:
  ```groovy
      static final List<String> DIRECTIVES = ['model', 'instruction', 'tools', 'maxIterations']
  ```
  with:
  ```groovy
      static final List<String> DIRECTIVES = ['model', 'instruction', 'goal', 'tools', 'maxIterations']
  ```

- [ ] **Add the `getGoal()` accessor** in `AgentDef.groovy`, immediately after `getInstruction()` (line 96):
  ```groovy
      String getGoal() { directives.get('goal') as String }
  ```

- [ ] **Run the unit test — expect PASS.** Command:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.AgentDefTest"
  ```
  Expected: PASS.

- [ ] **Add a loading test** in `AgentScriptLoadingTest.groovy` that a script declaring `goal` resolves and round-trips. Follow that file's existing pattern (it loads an agent script via `ScriptLoaderV2` and inspects the registered `AgentDef`). The agent body to use:
  ```
  agent qa {
      model 'openai/gpt-5-mini'
      instruction 'be concise'
      goal 'answer the question accurately'
      input:
          question: String
      output:
          answer: String
      prompt:
      """
      ${question}
      """
  }
  ```
  Assert the loaded `AgentDef.getGoal() == 'answer the question accurately'`. (Match the file's existing load/assert helpers exactly — read it first.)

- [ ] **Run the loading test — expect PASS.** Command:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
  ```
  Expected: PASS — confirms `goal` resolves through `AgentDsl` and is captured through `AgentBuilder`.

- [ ] **Commit.**
  ```bash
  git add modules/nf-lang/src/main/java/nextflow/script/dsl/AgentDsl.java \
          modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy \
          modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy \
          modules/nextflow/src/test/groovy/nextflow/script/AgentDefTest.groovy \
          modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy
  git commit -s -m "feat: add advisory \`goal\` agent directive

Declare goal(String) in AgentDsl (resolve time), add it to
AgentBuilder.DIRECTIVES (runtime capture), and expose AgentDef.getGoal().
A String directive shaped like instruction; no grammar/AST change."
  ```

---

## Task 2: Carry `goal` on `AgentRunnerRequest` and pass it from `AgentDef`

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/agent/AgentRunnerRequest.groovy`
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` (`run()`: capture ~line 170, construct ~line 195)
- Test: `modules/nextflow/src/test/groovy/nextflow/agent/AgentRunnerRequestTest.groovy` (create) and `modules/nextflow/src/test/groovy/nextflow/agent/AgentRunIntegrationTest.groovy` (extend — request carries goal)

**Interfaces:**
- Consumes: `AgentDef.getGoal()` (Task 1).
- Produces: `AgentRunnerRequest.goal : String` (last `@Canonical` field); `AgentDef.run` builds the request via named-args including `goal`.

### Steps

- [ ] **Write the failing field-order/named-args test.** Create `modules/nextflow/src/test/groovy/nextflow/agent/AgentRunnerRequestTest.groovy`:
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

  import spock.lang.Specification

  class AgentRunnerRequestTest extends Specification {

      def 'should build via named args with goal as the last field'() {
          when:
          def req = new AgentRunnerRequest(
              model: 'openai/gpt-5-mini',
              instruction: 'sys',
              prompt: 'p',
              maxIterations: 7,
              tools: [],
              outputSchema: null,
              inputJson: '{}',
              toolSpecs: null,
              dispatch: null,
              requestTimeoutSeconds: 30,
              goal: 'reach the objective')

          then:
          req.model == 'openai/gpt-5-mini'
          req.instruction == 'sys'
          req.prompt == 'p'
          req.maxIterations == 7
          req.inputJson == '{}'
          req.requestTimeoutSeconds == 30
          req.goal == 'reach the objective'
      }

      def 'should default goal to null when omitted'() {
          when:
          def req = new AgentRunnerRequest(model: 'm', prompt: 'p')
          then:
          req.goal == null
      }
  }
  ```

- [ ] **Run it — expect FAIL.** Command:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.AgentRunnerRequestTest"
  ```
  Expected: FAIL — `goal` is not a property of `AgentRunnerRequest` (`MissingPropertyException` / no such named arg).

- [ ] **Append the `goal` field** in `AgentRunnerRequest.groovy`. Change the field block (lines 44-53) by adding `goal` as the LAST field:
  ```groovy
  @Canonical
  @CompileStatic
  class AgentRunnerRequest {
      String model
      String instruction
      String prompt
      int maxIterations
      List tools
      Map outputSchema
      String inputJson
      List<ToolDescriptor> toolSpecs
      ToolDispatcher dispatch
      int requestTimeoutSeconds
      String goal
  }
  ```
  And update the positional-order javadoc note (lines 32-33) to:
  ```groovy
   * Being {@code @Canonical}, the positional constructor order is:
   * {@code (model, instruction, prompt, maxIterations, tools, outputSchema, inputJson, toolSpecs, dispatch, requestTimeoutSeconds, goal)}.
  ```
  Also add `goal` to the class-doc sentence (line 26-27 area) — append ", and the optional high-level goal" before the period.

- [ ] **Run the request test — expect PASS** (this also proves named-args construction works with `@Canonical`). Command:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.AgentRunnerRequestTest"
  ```
  Expected: PASS.
  **If the named-args form fails at runtime** (e.g. `@Canonical` shadows the map constructor), fall back to the positional constructor in Task-2's `AgentDef` edit below (append `agentGoal` as the 11th positional arg) and keep this test using positional construction; note the deviation in your report. (Confidence is high it works: `AgentDef.run` is `@CompileDynamic` and `model` is `String`, so a Map literal cannot bind to the tuple constructor's first param and Groovy uses the map-constructor path.)

- [ ] **Capture and pass `goal` in `AgentDef.run()`.** After `final agentInstruction = this.instruction` (line 170), add:
  ```groovy
          final agentGoal = this.goal
  ```
  Replace the positional request construction (line 195):
  ```groovy
          final req = new AgentRunnerRequest(agentModel, agentInstruction, promptText, agentMaxIter, agentTools, outputSchema, inputJson, toolSpecs, (bridge as ToolDispatcher), agentTimeoutSecs)
  ```
  with the named-args form (adds `goal`):
  ```groovy
          final req = new AgentRunnerRequest(
              model: agentModel,
              instruction: agentInstruction,
              prompt: promptText,
              maxIterations: agentMaxIter,
              tools: agentTools,
              outputSchema: outputSchema,
              inputJson: inputJson,
              toolSpecs: toolSpecs,
              dispatch: (bridge as ToolDispatcher),
              requestTimeoutSeconds: agentTimeoutSecs,
              goal: agentGoal)
  ```

- [ ] **Extend the integration test** so the request carries the goal. In `AgentRunIntegrationTest.groovy`, follow its existing pattern (it injects a capturing test runner via `AgentRunnerProvider` and runs an agent script). Add a case whose agent declares `goal 'do the thing'` and assert the captured `AgentRunnerRequest.goal == 'do the thing'`. (Read the file first; reuse its existing test-runner seam and script-run helper — do not invent a new harness.)

- [ ] **Run the integration test — expect PASS.** Command:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.AgentRunIntegrationTest"
  ```
  Expected: PASS — `AgentDef.run` propagates the declared goal into the request.

- [ ] **Commit.**
  ```bash
  git add modules/nextflow/src/main/groovy/nextflow/agent/AgentRunnerRequest.groovy \
          modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy \
          modules/nextflow/src/test/groovy/nextflow/agent/AgentRunnerRequestTest.groovy \
          modules/nextflow/src/test/groovy/nextflow/agent/AgentRunIntegrationTest.groovy
  git commit -s -m "feat: carry agent goal on AgentRunnerRequest

Append goal as the last @Canonical field; AgentDef.run captures the
directive and builds the request via named-args (non-breaking for future
fields). SPI stays langchain4j-free (plain String)."
  ```

---

## Task 3: Fold `goal` into the system message in the runner; document it

**Files:**
- Modify: `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy` (add `composeSystemMessage`; use it in `runWithTools` seed and in `composeMessages`)
- Modify: `docs/agent.mdx` (Directives table row + short paragraph)
- Test: `plugins/nf-agent/src/test/nextflow/agent/LangChainAgentToolLoopTest.groovy` (goal in system message; goal-only; neither) and `plugins/nf-agent/src/test/nextflow/agent/LangChainAgentRunnerTest.groovy` (single-shot path carries goal)

**Interfaces:**
- Consumes: `AgentRunnerRequest.goal` (Task 2).
- Produces: `LangChainAgentRunner.composeSystemMessage(AgentRunnerRequest) : String` (private static) — combined instruction+goal or `null`.

### Steps

- [ ] **Write the failing test (tool path).** Add to `LangChainAgentToolLoopTest.groovy` a spec that, with `instruction` AND `goal` set, asserts the first request's single `SystemMessage` text contains BOTH. Follow the file's existing mock pattern (closure-coerced `ChatModel` capturing `capturedRequests`, a `ToolDispatcher` stub, a `Stub(ChatModelFactory)`). Use a model that answers immediately (no tool call) to keep it focused:
  ```groovy
  def 'should fold goal into the single system message (tool path)'() {
      given:
      List<dev.langchain4j.model.chat.request.ChatRequest> capturedRequests = []
      dev.langchain4j.model.chat.ChatModel model = [
          chat: { dev.langchain4j.model.chat.request.ChatRequest req ->
              capturedRequests << req
              dev.langchain4j.model.chat.response.ChatResponse.builder()
                  .aiMessage(dev.langchain4j.data.message.AiMessage.from('done')).build()
          }
      ] as dev.langchain4j.model.chat.ChatModel
      ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
      def factory = Stub(ChatModelFactory) { createModel(_, _, _) >> model }
      def runner = new LangChainAgentRunner(modelFactory: factory)
      def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
      def req = new AgentRunnerRequest(
          model: 'openai/gpt-5-mini', instruction: 'You are careful.', prompt: 'go',
          maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
          toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0,
          goal: 'assemble then QC')

      when:
      def answer = runner.run(req)

      then:
      answer == 'done'
      def sys = capturedRequests[0].messages().find { it instanceof dev.langchain4j.data.message.SystemMessage } as dev.langchain4j.data.message.SystemMessage
      capturedRequests[0].messages().count { it instanceof dev.langchain4j.data.message.SystemMessage } == 1
      sys.text().contains('You are careful.')
      sys.text().contains('assemble then QC')
  }

  def 'should produce a system message from goal alone (no instruction)'() {
      given:
      List<dev.langchain4j.model.chat.request.ChatRequest> capturedRequests = []
      dev.langchain4j.model.chat.ChatModel model = [
          chat: { dev.langchain4j.model.chat.request.ChatRequest req ->
              capturedRequests << req
              dev.langchain4j.model.chat.response.ChatResponse.builder()
                  .aiMessage(dev.langchain4j.data.message.AiMessage.from('ok')).build()
          }
      ] as dev.langchain4j.model.chat.ChatModel
      ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
      def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) { createModel(_, _, _) >> model })
      def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
      def req = new AgentRunnerRequest(
          model: 'openai/gpt-5-mini', instruction: null, prompt: 'go',
          maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
          toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0,
          goal: 'reach the objective')

      when:
      runner.run(req)

      then:
      def msgs = capturedRequests[0].messages()
      msgs.count { it instanceof dev.langchain4j.data.message.SystemMessage } == 1
      (msgs.find { it instanceof dev.langchain4j.data.message.SystemMessage } as dev.langchain4j.data.message.SystemMessage).text().contains('reach the objective')
  }

  def 'should seed no system message when neither instruction nor goal is set'() {
      given:
      List<dev.langchain4j.model.chat.request.ChatRequest> capturedRequests = []
      dev.langchain4j.model.chat.ChatModel model = [
          chat: { dev.langchain4j.model.chat.request.ChatRequest req ->
              capturedRequests << req
              dev.langchain4j.model.chat.response.ChatResponse.builder()
                  .aiMessage(dev.langchain4j.data.message.AiMessage.from('x')).build()
          }
      ] as dev.langchain4j.model.chat.ChatModel
      ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
      def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) { createModel(_, _, _) >> model })
      def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
      def req = new AgentRunnerRequest(
          model: 'openai/gpt-5-mini', instruction: null, prompt: 'go',
          maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
          toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0, goal: null)

      when:
      runner.run(req)

      then:
      capturedRequests[0].messages().count { it instanceof dev.langchain4j.data.message.SystemMessage } == 0
  }
  ```
  (If the existing tests already import these langchain4j classes at the top of the file, use the short names instead of the fully-qualified ones — match the file's import style.)

- [ ] **Run it — expect FAIL.** Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.LangChainAgentToolLoopTest"
  ```
  Expected: FAIL — the goal text is not present in the SystemMessage (the runner currently seeds only `request.instruction`); the goal-only case currently seeds no SystemMessage.

- [ ] **Add the `composeSystemMessage` helper** in `LangChainAgentRunner.groovy`, immediately before `composeUserText` (line 177):
  ```groovy
      /**
       * Compose the system message: the {@code instruction} (role/persona),
       * followed by an optional {@code goal} section that steers the multi-turn
       * loop toward an objective. Returns {@code null} when neither is set, so the
       * caller seeds no {@link SystemMessage}.
       */
      private static String composeSystemMessage(AgentRunnerRequest request) {
          final String instruction = request.instruction
          final String goal = request.goal
          if( !instruction && !goal )
              return null
          final StringBuilder sb = new StringBuilder()
          if( instruction )
              sb.append(instruction)
          if( goal ) {
              if( sb.length() > 0 )
                  sb.append('\n\n')
              sb.append('Goal:\n').append(goal)
                .append('\nYou are done when the goal is met; produce your final answer as plain text.')
          }
          return sb.toString()
      }
  ```

- [ ] **Use it in `runWithTools`.** Replace the memory seed (lines 135-138):
  ```groovy
          // Seed memory with ONLY the system message; the user text is passed to chat().
          final ChatMemory memory = MessageWindowChatMemory.withMaxMessages(Integer.MAX_VALUE)
          if( request.instruction )
              memory.add(SystemMessage.from(request.instruction))
  ```
  with:
  ```groovy
          // Seed memory with ONLY the system message (instruction + optional goal);
          // the user text is passed to chat().
          final ChatMemory memory = MessageWindowChatMemory.withMaxMessages(Integer.MAX_VALUE)
          final String systemMessage = composeSystemMessage(request)
          if( systemMessage )
              memory.add(SystemMessage.from(systemMessage))
  ```

- [ ] **Use it in `composeMessages`** (single-shot path). Replace the body (lines 189-195):
  ```groovy
      private static List<ChatMessage> composeMessages(AgentRunnerRequest request) {
          final List<ChatMessage> messages = new ArrayList<ChatMessage>()
          if( request.instruction )
              messages.add(SystemMessage.from(request.instruction))
          messages.add(UserMessage.from(composeUserText(request)))
          return messages
      }
  ```
  with:
  ```groovy
      private static List<ChatMessage> composeMessages(AgentRunnerRequest request) {
          final List<ChatMessage> messages = new ArrayList<ChatMessage>()
          final String systemMessage = composeSystemMessage(request)
          if( systemMessage )
              messages.add(SystemMessage.from(systemMessage))
          messages.add(UserMessage.from(composeUserText(request)))
          return messages
      }
  ```

- [ ] **Run the tool-path tests — expect PASS.** Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.LangChainAgentToolLoopTest"
  ```
  Expected: PASS (all prior tests plus the three new goal cases).

- [ ] **Add a single-shot goal test** in `LangChainAgentRunnerTest.groovy` (no tools → `runSingleShot`/`composeMessages`). Follow that file's existing mock pattern (it mocks `ChatModel.chat(List)` for the single-shot path) and assert the goal text appears in the seeded `SystemMessage`. Read the file first and mirror its existing single-shot test exactly, adding `goal: '...'` to the request and asserting the captured messages' `SystemMessage` contains the goal text.

- [ ] **Run the single-shot test — expect PASS.** Command:
  ```bash
  ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.LangChainAgentRunnerTest"
  ```
  Expected: PASS.

- [ ] **Document `goal` in `docs/agent.mdx`.** In the Directives table, add a row after the `instruction` row:
  ```markdown
  | `goal` | no | A high-level objective that steers the agent's multi-turn loop. Advisory: folded into the system message; the model is encouraged to keep working (calling tools) until the goal is met, while `maxIterations` remains the hard cap. |
  ```
  And add a short paragraph under the Directives section:
  ```markdown
  The optional `goal` directive states what the agent should accomplish, distinct from the per-input `prompt`. It is appended to the system message as a labelled `Goal:` section. It does not change the loop mechanics or raise `maxIterations`; it only guides the model toward an objective across turns.
  ```

- [ ] **Commit.**
  ```bash
  git add plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy \
          plugins/nf-agent/src/test/nextflow/agent/LangChainAgentToolLoopTest.groovy \
          plugins/nf-agent/src/test/nextflow/agent/LangChainAgentRunnerTest.groovy \
          docs/agent.mdx
  git commit -s -m "feat: fold agent goal into the system message

Add composeSystemMessage (instruction + optional labelled Goal section),
used by both the tool path (memory seed) and the single-shot path
(composeMessages). Exactly one SystemMessage, or none when both absent.
Document the advisory goal directive."
  ```

---

## Task 4: Milestone 2 gate

**Files:** none modified — verification only.

### Exit-criteria steps

- [ ] **Boundary grep: core has no langchain4j.** Command:
  ```bash
  grep -RIl 'dev\.langchain4j' modules/nextflow/src/main || echo "CLEAN: no langchain4j imports in core"
  ```
  Expected: `CLEAN: no langchain4j imports in core`

- [ ] **Boundary grep: plugin does not touch ProcessDef/Channel/Path.** Command:
  ```bash
  grep -RIn -E 'import nextflow\.(script\.ProcessDef|extension\.|Channel)|import java\.nio\.file\.Path' plugins/nf-agent/src/main || echo "CLEAN: plugin does not touch ProcessDef/Channel/Path"
  ```
  Expected: `CLEAN: plugin does not touch ProcessDef/Channel/Path`

- [ ] **Guard: executeToolsConcurrently never called.** Command:
  ```bash
  grep -RIn 'executeToolsConcurrently' plugins/nf-agent/src/main || echo "CLEAN: executeToolsConcurrently never called"
  ```
  Expected: `CLEAN: executeToolsConcurrently never called`

- [ ] **nf-lang + nextflow agent unit tests.** Command:
  ```bash
  ./gradlew :nf-lang:test :nextflow:test --tests "nextflow.script.AgentDefTest" --tests "nextflow.script.parser.v2.AgentScriptLoadingTest" --tests "nextflow.agent.AgentRunnerRequestTest" --tests "nextflow.agent.AgentRunIntegrationTest"
  ```
  Expected: `BUILD SUCCESSFUL`, all green.

- [ ] **Full nf-agent suite under smoke.** Command:
  ```bash
  NXF_SMOKE=1 ./gradlew :plugins:nf-agent:test
  ```
  Expected: `BUILD SUCCESSFUL`, 0 failures; credential-gated e2e skipped under smoke.

- [ ] **Gated E2E (only if `OPENAI_API_KEY` is set).** Run a goal-directed agent end-to-end. Command:
  ```bash
  if [ -n "$OPENAI_API_KEY" ]; then ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.AgentToolEndToEndTest" --tests "nextflow.agent.AgentEndToEndTest"; else echo "SKIP: OPENAI_API_KEY not set"; fi
  ```
  Expected: `BUILD SUCCESSFUL` (existing e2e still pass with the goal plumbing inert when no `goal` is declared), or `SKIP`.

- [ ] **Clean worktree** (`.superpowers/` is git-ignored scratch; `examples/agents/two-agents/` is pre-existing untracked). Command:
  ```bash
  git status --porcelain
  ```
  Expected: no modifications to tracked source files from this milestone.

### Milestone 2 exit criteria (all must hold)

1. `goal '...'` resolves (AgentDsl) and is captured (AgentBuilder.DIRECTIVES) → `AgentDef.getGoal()`.
2. `AgentRunnerRequest` carries `goal` as its last `@Canonical` field; `AgentDef.run` builds the request via named-args; the request carries the declared goal.
3. `composeSystemMessage` folds instruction + goal into exactly one `SystemMessage` (or none) for BOTH the tool path and the single-shot path.
4. `goal` never raises `maxIterations`; no grammar/AST change; core stays langchain4j-free; plugin untouched boundary.
5. Docs document the advisory `goal` directive.
6. Full nf-agent suite green; e2e green when credentials present.
