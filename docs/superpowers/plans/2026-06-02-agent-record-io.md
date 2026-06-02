# Record-Typed Structured Agent I/O Implementation Plan (Plan F)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Replace the agent's free-style scalar I/O with **record-typed** I/O. The input record is bound into the prompt and serialized as JSON for the model; the output **record type's structure is reflected into a JSON schema** used as the LLM structured-output contract (langchain4j `responseFormat`), and the returned JSON is parsed and bound to an output **record** emitted on the channel.

**Architecture:** Named record types only (v1). No nf-lang *lowering* change is needed — a named record type already lowers to a concrete, field-introspectable JVM class at runtime (verified by probe). Work is: (1) resolution **requires** record-type I/O (rejects scalars + destructured `record(...)`); (2) a core **runtime schema deriver** (reflection on the output record `Class`) + extended `AgentRunnerRequest`; (3) `AgentDef.run` serializes the input record to JSON, derives the output schema, calls the runner, and binds the returned JSON to an output record via `TypeHelper.asRecordType`; (4) the plugin's runner sets `responseFormat` from the schema and returns raw JSON text (SPI unchanged); (5) migrate tests/example/docs to records.

**Tech Stack:** Groovy/Java (nf-lang resolution, core runtime, type reflection), langchain4j 1.15.1 structured outputs, Spock, Gradle.

---

## Verified facts this plan relies on (from research workflow + empirical probe, 2026-06-02)

- A named record type (`record Question { text: String; context: String? }`) used as agent I/O lowers to a **concrete class** at runtime: `agent.inputs[0].type == class Main.Question`, `agent.outputs[0].type == class Main.Answer`, with `getDeclaredFields()` returning the real fields (`[text, context]`, `[answer, confidence]`). **So schema can be derived at runtime by reflection; no `AgentToGroovyVisitorV2` change is needed for named types.**
- Destructured `record(...)` agent I/O is currently broken: input → name empty + bare `nextflow.script.types.Record` (no fields); output → dropped (`outputs == []`). v1 **rejects** it at resolution.
- Resolution does NOT currently constrain agent I/O types (scalars pass). `ScriptResolveVisitor.visitAgent` (`ScriptResolveVisitor.java:148-156`), `VariableScopeVisitor.visitAgent` (`VariableScopeVisitor.java:340-364`). Record-type detection: `cn.redirect() instanceof nextflow.script.ast.RecordNode` (named); the destructured form is a `TupleParameter` whose type equals `RECORD_TYPE` (`ClassHelper.makeCached(Record.class)`, constant already in `ScriptResolveVisitor`).
- Runtime record = `nextflow.util.RecordMap extends LinkedHashMap` (immutable), built via `nextflow.Nextflow.record(Map)`. `${q.text}` resolves because the bound delegate value is a Map.
- `nextflow.util.TypeHelper.asRecordType(Map, Class)` (nf-commons) validates required fields (via `@Nullable` reflection) and coerces+wraps into a `RecordMap`. Reusable for binding **provided** no `Path` output fields (its `Path` branch existence-checks and throws).
- langchain4j 1.15.1 (verified against jars in `plugins/nf-agent/build/target/libs/`): `ChatRequest.builder().messages(...).responseFormat(ResponseFormat.builder().type(ResponseFormatType.JSON).jsonSchema(JsonSchema).build()).build()`; `ChatModel.chat(ChatRequest)` is a `default` method → returns `ChatResponse`, `.aiMessage().text()`. `OpenAiChatModel.Builder` has `.responseFormat(ResponseFormat)`, `.strictJsonSchema(Boolean)`, `.supportedCapabilities(...)`. `JsonObjectSchema.builder().addStringProperty/addIntegerProperty/addNumberProperty/addBooleanProperty/.../required(List).build()`, `JsonSchema.builder().name(String).rootElement(JsonObjectSchema).build()`.
- The `@Nullable` annotation on optional record fields: confirm its FQN before use (`grep -rn "class Nullable" modules/`; expected `nextflow.script.types.Nullable` per research) — `TypeHelper` already imports it.

---

## File structure

### Create (core — `modules/nextflow/src/main/groovy/nextflow/agent/`)
- `RecordSchema.groovy` — reflects an output record `Class` → portable JSON-schema `Map`; rejects `Path`/unsupported field types.

### Modify (core)
- `script/AgentDef.groovy` — `run()`: serialize input record → JSON; derive output schema; build extended request; parse + bind runner JSON → output `RecordMap`; emit it.
- `agent/AgentRunnerRequest.groovy` — add `Map outputSchema` and `String inputJson`.

### Modify (nf-lang)
- `control/ScriptResolveVisitor.java` — `visitAgent`: require each input/output be a named record type; reject scalars and destructured `record(...)` with clear errors.

### Modify (plugin — `plugins/nf-agent/`)
- `src/main/nextflow/agent/JsonSchemaMapper.groovy` (new) — portable schema `Map` → langchain4j `JsonObjectSchema`/`JsonSchema`.
- `src/main/nextflow/agent/ChatModelFactory.groovy` — `createModel` accepts an optional output schema → set `responseFormat`+`strictJsonSchema` on the OpenAI model.
- `src/main/nextflow/agent/LangChainAgentRunner.groovy` — compose user message (`prompt` + input JSON), pass schema to the model, return JSON text.

### Migrate (tests/docs)
- `modules/nf-lang/.../parser/AgentParserTest.groovy`, `modules/nextflow/.../parser/v2/AgentScriptLoadingTest.groovy`, `modules/nextflow/.../script/{AgentDefTest,AgentBuilderTest}.groovy`, `modules/nextflow/.../agent/AgentRunIntegrationTest.groovy`, `plugins/nf-agent/.../{ChatModelFactoryTest,LangChainAgentRunnerTest,AgentEndToEndTest}.groovy`, `docs/agent.md`.

---

## Task 1: Resolution requires record-typed agent I/O (nf-lang)

**Files:** modify `ScriptResolveVisitor.java`; add tests to `AgentParserTest.groovy`; migrate existing `AgentParserTest` fixtures.

- [ ] **Step 1.1: Add failing tests** to `AgentParserTest.groovy` using `check(...)`:
  - `'should accept record-typed agent I/O'`: a script with `record Question { text: String }`, `record Answer { answer: String }`, agent `input: q: Question` / `output: a: Answer` / prompt `${q.text}` → `errors.isEmpty()`.
  - `'should reject scalar agent input/output'`: agent with `input: question: String` / `output: answer: String` → errors non-empty, message mentions record type.
  - `'should reject destructured record agent I/O'`: agent with `input: record(text: String)` → errors non-empty, message mentions named record type.
  - Migrate the existing fixtures in this file (which use `String`) to record types so they still pass.

- [ ] **Step 1.2: Run → expect FAIL** (`./gradlew :nf-lang:test --tests "nextflow.script.parser.AgentParserTest"`).

- [ ] **Step 1.3: Implement the constraint** in `ScriptResolveVisitor.visitAgent`. After resolving each input/output type, require it be a named record type. Sketch (adapt to real helpers/imports — `RecordNode`, `RECORD_TYPE` already referenced in this file):

```java
@Override
public void visitAgent(AgentNode node) {
    for( var input : asFlatParams(node.inputs) ) {
        resolver.resolveOrFail(input.getType(), input);
        requireRecordType(input.getType(), input, "Agent input `" + input.getName() + "`");
    }
    resolver.visit(node.directives);
    resolveTypedOutputs(node.outputs);
    requireRecordOutputs(node.outputs);
    resolver.visit(node.outputs);
    resolver.visit(node.prompt);
}

private void requireRecordType(ClassNode type, ASTNode ctx, String label) {
    if( RECORD_TYPE.equals(type) )   // destructured record(...) form -> TupleParameter type is bare Record
        addError(label + " must use a named record type; destructured `record(...)` is not yet supported for agents", ctx);
    else if( !(type.redirect() instanceof RecordNode) )
        addError(label + " must be a record type (got `" + type.getNameWithoutPackage() + "`)", ctx);
}
```
For outputs, iterate `asBlockStatements(node.outputs)`, get each target's type the same way `resolveTypedOutputs` does, and apply `requireRecordType` (label "Agent output `<name>`"). Detect the destructured `record(...)` output (a `MethodCallExpression`/`record` call rather than a typed `VariableExpression`) and reject it too. Use the visitor's existing error-reporting (`addError`) and helpers; confirm exact method names against the class.

> The destructured detection for outputs: a named output is an `ExpressionStatement` whose expression is a `VariableExpression`/`AssignmentExpression` with a record-type target; a destructured output is a `record(...)` call. Reject the latter for v1.

- [ ] **Step 1.4: Run → expect PASS** (4+ tests). Then `./gradlew :nf-lang:test --tests "nextflow.script.control.*"` (no regression).

- [ ] **Step 1.5: Commit** (`git commit -s -m "feat: require record-typed agent I/O at resolution"`).

---

## Task 2: Core — schema deriver, extended request, bind in AgentDef.run

**Files:** create `RecordSchema.groovy`; modify `AgentRunnerRequest.groovy`, `AgentDef.groovy`; tests in `nextflow/agent/`.

- [ ] **Step 2.1: Failing unit test for `RecordSchema`** (`modules/nextflow/src/test/groovy/nextflow/agent/RecordSchemaTest.groovy`): declare a Groovy class with `@Nullable`-style fields OR load a record type via the script loader, and assert `RecordSchema.of(AnswerClass)` returns `[type:'object', properties:[answer:[type:'string'], confidence:[type:'number']], required:['answer','confidence'], additionalProperties:false]`. Include a case with an optional field (not in `required`) and a case asserting a `Path` output field throws a clear error.

> Prefer building the test record class via the script loader (as the probe did) so the real compiled record class + `@Nullable` annotations are exercised. Confirm the `@Nullable` FQN first.

- [ ] **Step 2.2: Implement `RecordSchema.groovy`** — `static Map of(Class recordType)` reflecting `recordType.getDeclaredFields()` (skip synthetic), mapping each field type → JSON-schema fragment, optionality via `field.isAnnotationPresent(<Nullable>)`. Type map: `String`→`{type:'string'}`; `Integer`/`int`/`Long`/`long`→`{type:'integer'}`; `Double`/`Float`/`BigDecimal`/`double`/`float`→`{type:'number'}`; `Boolean`/`boolean`→`{type:'boolean'}`; a nested record (`Record.isAssignableFrom(fieldType)`)→recurse; `List`/`Collection` (element type via generic signature, fallback string)→`{type:'array', items:<elem>}`. **Reject** `Path` and any unmapped type with `IllegalArgumentException` naming the field. Set `additionalProperties:false` and `required:[...non-optional field names...]`.

- [ ] **Step 2.3: Run → PASS** (`./gradlew :nextflow:test --tests "nextflow.agent.RecordSchemaTest"`).

- [ ] **Step 2.4: Extend `AgentRunnerRequest`** — add `Map outputSchema` and `String inputJson` to the `@Canonical` class (after `tools`). Update the ctor-order note (callers must pass them). 

- [ ] **Step 2.5: Update `AgentRunIntegrationTest`** to the record-typed world (failing): script declares record types, agent `input: q: Question`/`output: a: Answer`; the **mock runner returns a JSON string** `{"answer":"...","confidence":0.9}`; assert the emitted value is a record (`result.val instanceof Map`, `result.val.answer == '...'`) AND the captured request has a non-null `outputSchema` (with `properties.answer`) and an `inputJson` containing the input field value.

- [ ] **Step 2.6: Implement `AgentDef.run`** changes (it's already `@CompileDynamic`):
  - Resolve `outputClass = outputs[0].type`.
  - Derive `outputSchema = RecordSchema.of(outputClass)`.
  - In the mapper, per input item: bind the record under the input name (existing delegate trick) to render `promptText`; serialize the input record to JSON (`inputJson = toJson(item)` with `Path`→abs-string handling); build `new AgentRunnerRequest(model, instruction, promptText, maxIter, tools, outputSchema, inputJson)`; `final json = runner.run(req)`; parse `final map = new groovy.json.JsonSlurper().parseText(stripFences(json)) as Map`; bind `final rec = nextflow.util.TypeHelper.asRecordType(map, outputClass)`; return `rec`.
  - Emit via the existing `MapOp` → `ChannelOut`. Add a small `stripFences` helper (remove ```json fences if present) and a JSON serializer that renders `Path` as `toString()`.

> Confirm `TypeHelper.asRecordType(Map, Class)` is accessible from core (nf-commons is a core dep) and returns a `RecordMap`. If it can't coerce nested records / lists as needed, narrow v1 output types to flat scalars and note it.

- [ ] **Step 2.7: Run → PASS** (`AgentRunIntegrationTest`). Reconcile `AgentDefTest`/`AgentBuilderTest` if they construct `AgentRunnerRequest` or assume String output.

- [ ] **Step 2.8: Commit** (`git commit -s -m "feat: derive output schema and bind structured agent output to a record"`).

---

## Task 3: Plugin — structured output via responseFormat

**Files:** create `JsonSchemaMapper.groovy`; modify `ChatModelFactory.groovy`, `LangChainAgentRunner.groovy`; tests.

- [ ] **Step 3.1: Failing test for `JsonSchemaMapper`** (`plugins/nf-agent/src/test/nextflow/agent/JsonSchemaMapperTest.groovy`): assert `JsonSchemaMapper.toJsonSchema('Answer', [type:'object', properties:[answer:[type:'string'], confidence:[type:'number']], required:['answer']])` returns a langchain4j `JsonSchema` whose root is a `JsonObjectSchema` with the expected properties/required (introspect what the API exposes; at minimum assert non-null and `name == 'Answer'`).

- [ ] **Step 3.2: Implement `JsonSchemaMapper.groovy`** — portable `Map` → `dev.langchain4j.model.chat.request.json.JsonSchema`. Build a `JsonObjectSchema.builder()`, add a typed property per entry (`string`→`addStringProperty`, `integer`→`addIntegerProperty`, `number`→`addNumberProperty`, `boolean`→`addBooleanProperty`, `array`→`JsonArraySchema`, nested `object`→recurse), set `.required(requiredList)`, wrap in `JsonSchema.builder().name(name).rootElement(obj).build()`. Verify exact builder method names against the jar (`javap`).

- [ ] **Step 3.3: Run → PASS** (`./gradlew :plugins:nf-agent:test --tests "nextflow.agent.JsonSchemaMapperTest"`).

- [ ] **Step 3.4: Update `LangChainAgentRunnerTest`** (failing): the mock `ChatModel` must implement `chat(ChatRequest)` (the request-object overload), assert the request's `responseFormat` is JSON with a non-null `jsonSchema`, that the user message contains the prompt AND the input JSON, and return a `ChatResponse` whose `aiMessage().text()` is a canned JSON string; assert `run(req)` returns that JSON string.

- [ ] **Step 3.5: Implement plugin changes:**
  - `ChatModelFactory.createModel(String modelId, int timeoutSeconds, JsonSchema schema)` (overload or extra param): when `schema != null`, set `.responseFormat(ResponseFormat.builder().type(JSON).jsonSchema(schema).build())` and `.strictJsonSchema(true)` on the `OpenAiChatModel.Builder`. Keep the no-schema overload working.
  - `LangChainAgentRunner.run(AgentRunnerRequest req)`: if `req.outputSchema`, `schema = JsonSchemaMapper.toJsonSchema('Output', req.outputSchema)`; build model via the schema overload; compose the user message text = `req.prompt + (req.inputJson ? "\n\nInput (JSON):\n" + req.inputJson : "")`; messages = `[SystemMessage(instruction)?, UserMessage(userText)]`; `final resp = model.chat(ChatRequest.builder().messages(messages).build())` (responseFormat is on the model) — OR put responseFormat on the `ChatRequest` if that's the verified-cleaner path; return `resp.aiMessage().text()`.

> Decide responseFormat-on-model vs responseFormat-on-request based on which the jar makes work for OpenAI strict mode; the research notes the model-builder path is the reliable one. Keep the SPI return type `String`.

- [ ] **Step 3.6: Run → PASS** (`ChatModelFactoryTest`, `LangChainAgentRunnerTest`, `JsonSchemaMapperTest`).

- [ ] **Step 3.7: Commit** (`git commit -s -m "feat: nf-agent structured output via langchain4j responseFormat"`).

---

## Task 4: Migrate remaining tests, docs, and the gated E2E

- [ ] **Step 4.1: `AgentScriptLoadingTest`** — change the directive-rich fixture to record-typed I/O; assert the loaded `AgentDef.outputs[0].type` is the record class.
- [ ] **Step 4.2: `docs/agent.md`** — update the example + directives/IO sections to record types (`record` decls, `input: q: Question`, `output: a: Answer`), document the structured-output behavior and input-JSON presentation, and the v1 limits (named records only, no `Path` outputs, scalars-not-allowed).
- [ ] **Step 4.3: `AgentEndToEndTest`** (gated, real OpenAI) — build an `AgentRunnerRequest` with a real `outputSchema` (e.g. `{answer:string, capital:string}`) and assert the returned JSON parses to those fields (e.g. contains `paris`). Keep `@Requires({System.getenv('OPENAI_API_KEY')})`.
- [ ] **Step 4.4: Run** the full agent suite (keyless): `env -u OPENAI_API_KEY ./gradlew :nf-lang:test --tests "*Agent*" :nextflow:test --tests "nextflow.script.Agent*" --tests "nextflow.agent.*" --tests "*AgentScriptLoading*" :plugins:nf-agent:test --rerun-tasks` → all pass (E2E skipped).
- [ ] **Step 4.5: Commit** (`git commit -s -m "test+docs: migrate agent examples to record-typed structured I/O"`).

---

## Task 5: Verification + ADR status

- [ ] **Step 5.1: Full regression** — the command in 4.4 plus `:nf-lang:test --tests "nextflow.script.control.*"` and a real end-to-end `./launch.sh` run (with key) of a record-typed agent producing a structured record (capstone; report the JSON output).
- [ ] **Step 5.2: Update the ADR** `20260505-llm-agent-primitive.md` — mark the record-typed structured-I/O direction delivered; note remaining deferrals (destructured form, `Path`/enum/`Map` outputs, tool dispatch).
- [ ] **Step 5.3: Commit.**

---

## What this plan does NOT deliver
- Destructured `record(...)` agent I/O (rejected with a clear error; needs lowering work).
- `Path`/enum/`Map`/`Set` fields in agent **outputs**; multiple inputs/outputs.
- Tool dispatch (the separate ToolDispatcher bridge) and `agent {}` config scope.
- Input-side JSON *schema* injection (only the input JSON *value* is sent, per the 2026-06-02 decision).

## Self-Review
- **Spec coverage:** Implements the 2026-06-02 ADR direction: record-required I/O, output-record→JSON-schema→`responseFormat`→parsed→bound record; input record bound into prompt + serialized as JSON value (the chosen presentation). Each deferral is explicit and fails loudly.
- **Placeholder scan:** Code sketches are concrete; the ">" notes are verification reminders for the two uncertain surfaces (langchain4j `JsonObjectSchema` builder names; `TypeHelper` nested/list coercion + `@Nullable` FQN).
- **Type consistency:** `RecordSchema.of(Class)→Map`; `AgentRunnerRequest(model, instruction, prompt, maxIterations, tools, outputSchema, inputJson)`; `JsonSchemaMapper.toJsonSchema(String, Map)→JsonSchema`; `ChatModelFactory.createModel(String, int, JsonSchema)`; `AgentRunner.run(AgentRunnerRequest)→String` (unchanged); core binds via `TypeHelper.asRecordType(Map, Class)→RecordMap`. Used consistently across resolution, core, plugin, and tests.
