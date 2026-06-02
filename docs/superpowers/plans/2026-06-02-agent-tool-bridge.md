# Agent Tool Bridge — Module-as-LLM-Tool Implementation Plan (Plan G)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Let an agent invoke Nextflow modules/processes as **LLM tools that run as real dataflow nodes**. The LLM calls a tool; the harness marshals the call's JSON args into channel values, runs the module through the standard executor/container/retry/cache machinery, and serializes its outputs back to the LLM as JSON. This is the ADR's headline goal.

**Architecture (feasibility validated by two empirical spikes — see ADR §"Module-as-tool execution"):** the agent is an ordinary dataflow node. For each declared tool, the harness — at agent construction (inside the workflow body, **before** ignition) — pre-wires the tool process into the network: `toolOut = toolProcess(toolIn)`. Post-ignition, the agent operator drives the LLM loop and, per tool call, emits the marshalled args onto `toolIn` and **blocks** on `toolOut.val` (serialized → correct correlation; pool-safe). On completion the harness poisons `toolIn`. **No core runtime change is required.** Core owns resolution + wiring + dispatch + serialization (langchain4j-free, exposed via a `ToolDescriptor` Map DTO + a `ToolDispatcher` callback); the plugin owns the langchain4j tool-call loop.

**Tech Stack:** Groovy (core dataflow/runtime), langchain4j 1.x tool-calling (core jar 1.14.1, openai 1.15.1), GPars dataflow, Spock, Gradle.

---

## Validated facts this plan relies on (spikes + research, 2026-06-02)

- **Pre-wire + blocking-pull works under the real executor.** Recipe: `def toolIn = CH.queue(); def toolOut = CH.getReadChannel(toolDef.run([toolIn] as Object[])[0])` BEFORE `fireDataflowNetwork()`; post-fire, `toolIn.bind(arg); def r = toolOut.val`. The tool task runs in a real work dir; `r` is the genuine transformed output. (Spike #2.)
- **Termination obligation:** the harness MUST `toolIn.bind(groovyx.gpars.dataflow.operator.PoisonPill.instance)` once the agent finishes, or `session.await()` hangs. (Spike #2.)
- **Correlation:** serialize tool calls (one `toolIn.bind` → one `toolOut.val` at a time). Outputs reorder under `maxForks>1`; serialized blocking-pull is correct and simplest. Pool-safe (task bodies run on `session.execService`, separate from the GPars operator pool). (Spike #2.)
- **Post-fire node *creation* deadlocks** — do NOT call `ProcessDef.run()` from inside the agent operator at tool-call time; all wiring is pre-fire. (Spike #1.)
- **langchain4j tool API (verified in jars):** `ChatRequest.builder().messages(...).toolSpecifications(List<ToolSpecification>).build()`; `resp.aiMessage().hasToolExecutionRequests()` / `.toolExecutionRequests()` → each `ToolExecutionRequest` has `id()`, `name()`, `arguments()` (JSON string); feed results back via `ToolExecutionResultMessage.from(req, resultString)`; loop until no tool requests, then `aiMessage().text()`. `ToolSpecification.builder().name(..).description(..).parameters(JsonObjectSchema).build()`.
- **Current agent runtime:** `AgentDef.run` (`modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy`) runs a `MapOp` over the single input; `AgentRunnerRequest` carries model/instruction/prompt/maxIterations/tools/outputSchema/inputJson; `AgentRunner.run`→String; `LangChainAgentRunner` does single-shot chat. The `tools(Object...)` directive (`AgentDsl.java`) is captured but ignored. `JsonSchemaMapper` (plugin) already maps a portable schema Map → langchain4j `JsonSchema`.
- **Module resolution:** `nextflow.module.{ModuleReference,ModuleResolver,ModuleSpec,ModuleSpecFactory}`; `ModuleResolver.resolve(ref)` → `Path` to main.nf (NOT a ProcessDef — compilation to a ProcessDef is its own step, deferred to Phase 3). An in-scope process is reachable via `ScriptMeta` (confirm the exact lookup, e.g. `ScriptMeta.get(owner).getProcess(name)` / `getDefinitions()`).

---

## Phase 1 — Relax agent I/O to process-style (prerequisite, independent, low-risk)

Removes Plan-F's record-only restriction so an agent can be `input: val question / output: val answer`; keeps the record structured-output path as opt-in. Mirrors a typed process's `val`/`path` I/O.

### Task 1.1 — Resolution: allow `val`/`path` agent I/O again
**Files:** `modules/nf-lang/src/main/java/nextflow/script/control/ScriptResolveVisitor.java`; tests `AgentParserTest.groovy`.
- [ ] Failing tests: `'should accept val agent I/O'` (`input: question: String` / `output: answer: String`, prompt `${question}`) → `check(...).isEmpty()`; keep a record-typed accept test; keep the destructured-`record(...)` rejection. (Scalars are now ALLOWED again; only destructured `record(...)` stays rejected.)
- [ ] Implement: relax `requireRecordType`/`requireRecordOutputs` in `visitAgent` — permit any resolvable scalar/`path`/record type; reject ONLY the destructured `record(...)` form (the `RECORD_TYPE` bare-tuple input and the `record(...)` output method-call) with the existing message. Migrate the Plan-F fixtures that asserted scalar rejection.
- [ ] Run `:nf-lang:test --tests "*AgentParserTest" --tests "nextflow.script.control.*"` → green. Commit `feat: allow process-style (val/path) agent I/O; keep record output opt-in`.

### Task 1.2 — Runtime: optional structured output
**Files:** `AgentDef.groovy`, `AgentBuilder.groovy`; tests `AgentRunIntegrationTest.groovy`, `AgentScriptLoadingTest.groovy`.
- [ ] Failing test: a `val`-output agent whose mock runner returns plain text emits that text verbatim (no JSON parse / no `asRecordType`); a record-output agent still binds structured JSON (keep one such test).
- [ ] Implement in `AgentDef.run`: derive `outputSchema` and parse+bind via `TypeHelper.asRecordType` **only when** `outputs[0].type` is a record type; otherwise pass `outputSchema=null` and emit the runner's raw `String` (current behavior pre-Plan-F). Broaden `AgentBuilder._output_`/`AgentInput`/`AgentOutput` to carry non-record types. Keep `normalizeForJson` for input serialization.
- [ ] Plugin: `LangChainAgentRunner` already sends `responseFormat` only when `outputSchema != null` — verify the null path is single-shot text. 
- [ ] Run the full agent suite keyless → green. Migrate docs example to `val` I/O. Commit `feat: make agent structured-output binding opt-in (record output only)`.

---

## Phase 2 — Headline slice: invoke an IN-SCOPE process as a tool, end-to-end

Proves the headline with a **local process tool** (no registry yet). Example target:
```groovy
process greet {
    input:  name: String
    output: String
    exec:   "Hello ${name}!"          // or the codebase's typed exec form
}

agent assistant {
    model 'openai/gpt-5-mini'
    instruction 'Use the greet tool to greet the user by name, then report what it returned.'
    tools 'greet'                      // names an in-scope process
    input:  request: String
    output: String
    prompt: "${request}"
}

workflow {
    assistant(channel.of('greet Ada'))
        .view { it }
}
```

### Task 2.1 — Core SPI: `ToolDescriptor` + `ToolDispatcher` + request fields
**Files:** new `modules/nextflow/src/main/groovy/nextflow/agent/ToolDescriptor.groovy`, `ToolDispatcher.groovy`; modify `AgentRunnerRequest.groovy`; tests.
- [ ] `ToolDescriptor` — `@Canonical` DTO: `String name`, `String description`, `Map inputSchema`, `Map outputSchema` (plain Maps; langchain4j-free).
- [ ] `ToolDispatcher` — `@FunctionalInterface interface ToolDispatcher { String call(String toolName, String argsJson) }`.
- [ ] `AgentRunnerRequest` — add `List<ToolDescriptor> toolSpecs` and `ToolDispatcher dispatch` (after existing fields). Note: these are in-JVM only (not serialized).
- [ ] Unit test the DTOs; commit `feat: add ToolDescriptor and ToolDispatcher agent SPI`.

### Task 2.2 — Core: tool input/output schema from a process's typed I/O
**Files:** new `modules/nextflow/src/main/groovy/nextflow/agent/ProcessToolSchema.groovy` (or extend `RecordSchema`); tests.
- [ ] `static Map inputSchema(ProcessDef)` and `static Map outputSchema(ProcessDef)` deriving portable JSON-schema Maps from the process's declared typed inputs/outputs. For Phase 2 support scalar `val`/`String`/number/boolean inputs+outputs; map each input param `name:Type`. (Reuse the `RecordSchema` type→fragment table.) Defer `tuple`/`path`/`meta` to Phase 3 (throw a clear "not yet supported as a tool input" error for those, so failures are loud).
- [ ] Find the ProcessDef typed-I/O accessors (the declared inputs/outputs the lowering produced — `ProcessDef`/`ProcessConfig`/`ProcessInputsDef`). Confirm names by reading the process runtime.
- [ ] Tests: load a script with `process greet { input: name: String; output: String; ... }`, assert `inputSchema` == `{type:object, properties:{name:{type:string}}, required:[name]}`. Commit `feat: derive tool JSON-schema from a process typed I/O`.

### Task 2.3 — Core: `ModuleToolBridge` (pre-wire + dispatch + poison)
**Files:** new `modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy`; modify `AgentDef.groovy`; integration test.
- [ ] `ModuleToolBridge`: constructed in `AgentDef.run` (pre-fire) with the resolved tool `ProcessDef`s. For each tool: build a `ToolDescriptor` (name, description from the process, `inputSchema`/`outputSchema` via Task 2.2); create `toolIn = CH.queue()`, `toolOut = CH.getReadChannel(toolDef.cloneWithName(name).run([toolIn] as Object[])[0])` (clone to get a stable invocable; confirm single-output assumption for Phase 2). Hold `{name → (toolIn, toolOut)}`.
  - `dispatch(name, argsJson)`: parse JSON → marshal to the channel value matching the process input (Phase 2: single scalar → the value itself) → `toolIn.bind(value)` → `result = toolOut.val` → serialize output to a JSON object (Phase 2: `{<outName>: value}`) → return JSON string. Serialized (synchronized) so calls don't interleave.
  - `close()`: `toolIn.bind(PoisonPill.instance)` for every tool.
- [ ] `AgentDef.run`: resolve `tools` entries to in-scope `ProcessDef`s (Phase 2: each tools entry is a String naming a process; look it up via `ScriptMeta.get(owner).getProcess(name)` — confirm API). If `tools` non-empty: build the `ModuleToolBridge` (pre-fire), set `request.toolSpecs` + `request.dispatch`, and ensure `bridge.close()` is called when the agent's `MapOp` source completes (hook the operator `onComplete`/a terminal subscriber). If `tools` empty, behave as today.
- [ ] Integration test (no real LLM): use the `AgentRunnerProvider.testRunner` seam with a mock runner that, given `request.dispatch`, calls `dispatch.call('greet', '{"name":"Ada"}')` and asserts it returns `{"...":"Hello Ada!"}` after the **real** `greet` process executed (assert via the transformed value). This proves the full pre-wire+dispatch+execute path end-to-end WITHOUT the LLM. Use `ScriptHelper.runScript` + `@Timeout`. Commit `feat: ModuleToolBridge — run in-scope process as a tool via pre-wired dataflow`.

### Task 2.4 — Plugin: langchain4j tool-call loop
**Files:** `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy`, extend `JsonSchemaMapper.groovy`, new `ModuleToolAdapter.groovy`; tests.
- [ ] `ModuleToolAdapter.toToolSpecification(ToolDescriptor)`: `inputSchema` Map → `JsonObjectSchema` (reuse/extend `JsonSchemaMapper`) → `ToolSpecification`.
- [ ] `LangChainAgentRunner.run`: when `request.toolSpecs` non-empty, run the manual loop (≤ maxIterations): `chat(ChatRequest.builder().messages(messages).toolSpecifications(specs).build())`; if `aiMessage().hasToolExecutionRequests()`, append the Aی message, then for each request call `request.dispatch.call(req.name(), req.arguments())`, append `ToolExecutionResultMessage.from(req, result)`, and loop; else return `aiMessage().text()`. When empty, keep single-shot.
- [ ] Test with a mock `ChatModel` that emits one tool-execution request then a final text answer; a stub `ToolDispatcher` returns a canned JSON; assert the loop calls the dispatcher with the right args and returns the final text. Commit `feat: langchain4j tool-call loop in nf-agent runner`.

### Task 2.5 — Verification + real end-to-end run
- [ ] Full agent suite keyless → green. 
- [ ] Real `./launch.sh` run (dev plugin + `OPENAI_API_KEY`) of the `assistant`/`greet` example; capture the output proving the LLM called `greet`, the process ran (work dir), and the result returned to the LLM. Commit a docs example + ADR status update.

---

## Phase 3 — External/registry modules + `meta.yml`-driven schema

**Spike (2026-06-02) confirmed feasibility:** a module FILE compiles to a runnable `ProcessDef` at agent-run time via `new ScriptLoaderV2(session).setModule(true); loader.parse(path); loader.runScript(); ScriptMeta.get(loader.script).getProcess(name)` — on the OWNER's live `Session`, pre-ignition (where `createToolBridge` runs). `setModule(true)` is mandatory (else the standalone process auto-builds an entry workflow and fails on `--params`). The resulting `ProcessDef` pre-wires through the existing `ModuleToolBridge` unchanged. No core change needed.

### Task 3.1 — External module file as a tool (compile → ProcessDef)
- In `AgentDef.createToolBridge`, a `tools` entry that is NOT an in-scope process name and looks like a path (`./x.nf`, `/x.nf`, `x.nf`) is compiled via the spike recipe → `ProcessDef`, added to the bridge map. In-scope names still resolve via `ScriptMeta.getProcess`. Tool name = sanitized basename/process name.
- Test: a LOCAL module file with a typed scalar process used as `tools './mod.nf'`, end-to-end via the no-LLM dispatch test (reusing `AgentToolBridgeIntegrationTest` style); assert the external process ran + result returned + session terminates.

### Task 3.2 — `ModuleSpec`/`meta.yml`-driven tool schema + tuple/path/map marshalling
- Derive the tool input/output schema from a module's `ModuleSpec` (`nextflow.module.ModuleSpecFactory.fromYaml(sibling meta.yml)`) when present — this is the rich "free" descriptor (names, types, **descriptions**). New `ModuleSpecToolSchema` (or extend `ProcessToolSchema`): flatten the input tuple `{meta:map, reads:file}` to top-level properties (`meta`→object, `file`/`path`→string path-handle, `val`→typed, with the item `description`); `required` = non-optional items.
- Extend `ModuleToolBridge` marshalling: LLM args (flattened JSON) → the tuple channel value the `ProcessDef` expects (`[meta, file(readsPath)]`); `file`/`path` args → a `Path` (via `Nextflow.file`/`FileHelper`). Serialize tuple/file outputs back to JSON (files → absolute path-handle strings; `meta`→object; `eval`/topic → string). Decide `meta`-map provenance: LLM-supplied (documented in the system prompt) for v1.
- Handle optional inputs (`path(..., optional:true)` / `?` → not in `required`).
- Test: a LOCAL module file + sibling `meta.yml` mimicking the fastqc tuple shape (trivial `exec` body, no real container); assert the tuple-input tool schema (flattened + descriptions), arg marshalling, and tuple/file output serialization. The opaque-path contract is stated in the agent system prompt.

### Task 3.3 — Registry reference resolution (`nf-core/fastqc`)
- A `tools 'nf-core/fastqc'` entry → `ModuleReference.parse` + `ModuleResolver.resolve(ref, autoInstall=true)` → `Path` (downloaded) → the 3.1 compile path; fetch the registry `ModuleMetadata` (or the installed `meta.yml`) for the 3.2 schema.
- **External dependency:** registry download needs network + (likely) registry auth, and running a real nf-core tool needs a container runtime — so the REAL end-to-end (download + run `fastqc`) is a validation-environment concern. Build + unit/integration-test the resolution+compile CODE PATH with a LOCAL stand-in; gate any real-registry test (skip without registry access) and clearly LOG/flag the requirement. Do not fake success.

## Phase 4 — Polish (OUTLINE)

- Multiple tools per agent; correlation IDs to allow non-serialized/parallel tool calls; error-as-tool-result (return an error JSON so the LLM can recover) counting toward `maxIterations`.
- Concurrency when the agent input is a queue (multiple in-flight LLM loops): per-loop correlation, thread-safety of `ScriptMeta`/`Global.session` access, and per-tool `toolIn` sharing.
- Termination edge cases (poison on error paths); `agent {}` config scope; docs + gated real-LLM E2E with a real module.

### Task 4.1 — Multi-tool support + dispatch-level error-as-tool-result (DELIVERED)
- **Multi-tool** already worked end-to-end (no fix needed): `AgentDef.createToolBridge` resolves every `tools` entry into the `ModuleToolBridge` tool map; `descriptors()` returns all; `LangChainAgentRunner.runWithTools` advertises every `ToolSpecification` and dispatches each `ToolExecutionRequest` by name. Covered by a no-LLM integration test (`AgentMultiToolBridgeIntegrationTest`): two in-scope processes (`shout`/`whisper`) wired as tools, both descriptors advertised, each dispatched to its correct distinct process.
- **Dispatch-level error-as-tool-result**: `ModuleToolBridge.call` returns a well-formed `{"error": "<message>"}` JSON string (naming the tool + the problem) instead of throwing for: an unknown tool name, an unparseable/blank/non-object `argsJson`, and any argument-marshalling failure. The error is logged at warn and counts as a normal tool result so the LLM can retry rather than the agent loop aborting.
- **Limitation (documented, future hardening):** a failure of the underlying tool *process task* is NOT recovered as a tool result — it surfaces through the process's own dataflow operator and aborts the session via the standard error model. Intercepting it cleanly fights the runtime and is deferred. Only dispatch-level errors (unknown tool / arg parsing / arg marshalling) are turned into a tool-result error.

---

## What this plan does NOT deliver (v1 / Phases 1–2)
Registry module resolution + `meta.yml` schemas (Phase 3); tuple/`path`/`meta`/`eval` tool I/O (Phase 3); multi-tool, parallel/correlated calls, queue-input concurrency, error recovery (Phase 4). Phases 1–2 prove the headline with a single in-scope scalar-I/O process tool, end to end.

## Self-Review
- **Spec coverage:** Implements the ADR headline (tool call = real dataflow node; module-as-tool) via the spike-validated pre-wire mechanism, and the agent-I/O relaxation the user requested. Registry/`meta.yml` schemas and richer I/O are staged to Phase 3 with loud "not yet supported" failures meanwhile.
- **Placeholder scan:** Code recipes are concrete; the ">"/"confirm" notes are verification points (ProcessDef typed-I/O accessors; `ScriptMeta` process lookup; the operator onComplete hook) the implementer must verify against the runtime, not deferred work.
- **Type consistency:** `ToolDescriptor(name, description, Map inputSchema, Map outputSchema)`; `ToolDispatcher.call(String,String)→String`; `AgentRunnerRequest(+toolSpecs:List<ToolDescriptor>, +dispatch:ToolDispatcher)`; `ModuleToolBridge.dispatch/close`; `ProcessToolSchema.inputSchema/outputSchema(ProcessDef)→Map`; `ModuleToolAdapter.toToolSpecification(ToolDescriptor)→ToolSpecification`. Used consistently across core, plugin, and tests. The pre-wire recipe matches the spike-verified `CH.queue()` / `toolDef.run([toolIn])[0]` / `toolOut.val` / `PoisonPill` sequence.
