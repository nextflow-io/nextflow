# Nextflow Agent — langchain4j AiServices Refactor

- Author: Paolo Di Tommaso
- Status: draft (design approved in principle; pending written-spec review)
- Date: 2026-06-22
- Branch: `agent-syntax-and-model`
- Supersedes (engine layer only): the hand-rolled `ChatModel` tool loop in `nf-agent`
- Related: [`adr/20260505-llm-agent-primitive.md`](../../../adr/20260505-llm-agent-primitive.md), [`docs/agent.mdx`](../../agent.mdx)

## Summary

Refactor the agent runtime engine onto langchain4j's **`AiServices`** API (core
artifact), in three independently-shippable milestones:

1. **Engine swap** — replace the hand-rolled chat→tool→chat loop with an
   `AiServices` proxy. Plugin-only; behaviour-preserving.
2. **Goal-directed loop** — one optional `goal` directive folded into the system
   message to steer the existing multi-turn loop toward a declared objective.
3. **Target refactor** — a per-invocation **work-dir sandbox**, a generic
   **`filesystem`** tool (read/write/list/exists), and a single generic
   **`module_run`** tool whose available modules are inferred from the script's
   `include` statements. The `tools` directive lists **generic capabilities only**.

The agent stays a process-shaped primitive (one record in → one record out);
Nextflow channels/workflows remain the multi-agent composition layer. We do **not**
adopt `langchain4j-agentic` / `AgenticServices` — its orchestration overlaps
Nextflow's own composition and is a non-goal per the ADR.

## Motivation

The current `nf-agent` plugin (`LangChainAgentRunner.runWithTools`) hand-rolls the
tool-call loop with the low-level `ChatModel` API: it advertises tool specs on every
request, inspects `AiMessage.toolExecutionRequests()`, dispatches each, appends
`ToolExecutionResultMessage`, and loops to an iteration cap. This is exactly what
`AiServices` does for free, with memory and structured-output support on top.

Moving to `AiServices`:

- deletes the manual loop, message bookkeeping, and cap handling;
- gives per-invocation chat memory and (later) structured-output-with-tools;
- keeps the proven core/plugin split and the module-as-tool binding rules untouched;
- positions the feature on the framework's supported, evolving API.

This refactor is an **engine swap plus two additive capabilities** — not a redesign
of the agent model, the dataflow tool bridge, or the I/O binding rules.

## Locked decisions

1. **`AiServices` (core `langchain4j`), single agent.** No `AgenticServices`
   orchestration; channels stay the composition layer.
2. **Per-agent work-dir sandbox**, like any Nextflow task. The `filesystem` tool may
   read/write that work dir and read the output work dirs of modules the agent ran.
3. **`tools` lists generic capabilities only** (`'filesystem'`, `'module_run'`;
   `'bash'` deferred). The **specific modules are inferred solely from the script's
   `include` statements** — never enumerated in `tools`. Enabling `'module_run'`
   exposes **one** generic `module_run(module, args)` tool whose `module` enum is the
   set of `include`d modules.
4. **Module execution stays native**, via the existing `ModuleToolBridge` (run a
   process as a real dataflow node) and its binding rules — now invoked behind the
   single `module_run` tool. `bash` and a `nextflow run` CLI tool are **deferred**.

## Architecture after the refactor

### Two-layer split (hard invariant — verified clean today)

`grep -r dev.langchain4j modules/` is empty and must stay empty.

- **Core (`modules/nextflow`, `modules/nf-lang`)** — owns the DSL, AST, dataflow
  wiring, and tool *execution*. Never imports `dev.langchain4j.*`. Exposes the SPI:
  - `AgentRunnerRequest` — plain fields (model, instruction, prompt, goal,
    maxIterations, outputSchema, inputJson, toolSpecs, dispatch, timeout, work dir).
  - `ToolDescriptor` — a Map-carrying POJO (name, description, input schema, `type`).
  - `ToolDispatcher` — the SAM `String call(String name, String argsJson)`.
  `ModuleToolBridge implements ToolDispatcher` is the single execution chokepoint;
  its `call()` is `synchronized` (keeps module tool input/output correlated).
- **Plugin (`nf-agent`)** — owns everything langchain4j. Never touches
  `ProcessDef`/`Channel`/`Path`. Maps each `ToolDescriptor` → `ToolSpecification`
  (`ModuleToolAdapter`/`JsonSchemaMapper`), wraps the single `dispatch` callback in
  `ToolExecutor`s, and runs the `AiServices` proxy.

The plugin sees a **uniform descriptor list and one dispatch callback** and cannot
tell `module_run` from `filesystem` — all `type` discrimination lives in core's
`ModuleToolBridge.call()`.

### Where AiServices sits

Inside `LangChainAgentRunner.runWithTools()` (plugin), constructed **fresh per input
record** (no caching → no cross-record state leakage). It receives a `ChatModel`, a
`Map<ToolSpecification, ToolExecutor>`, and a `MessageWindowChatMemory`. Every tool
call routes through the one `request.dispatch.call(...)` callback.

**Critical concurrency invariant (verified against langchain4j 1.16.3 bytecode):**
the design never calls `executeToolsConcurrently()`, so AiServices' `executor` field
stays `null` and **all tool calls run sequentially on the calling GPars operator
thread** — identical to today's manual loop. This is what preserves the pre-wired,
serialized-dispatch contract. Treat "never call `executeToolsConcurrently()`" as a
documented hard invariant with a test asserting tools execute on the calling thread.

### Pre-wiring constraint

Module tools are wired into the dataflow graph as operators **before network
ignition** (`AgentDef.createToolBridge()` runs in the workflow body; post-fire node
creation deadlocks). The **full module set must be known at agent-construction time**
— so when `module_run` is enabled, **every `include`d module is pre-wired** as a
candidate (the LLM may select any of them; `module_run` dispatch routes to the chosen
one by name). No mid-loop langchain4j `ToolProvider`. The `filesystem` capability adds
no dataflow node.

### Tool taxonomy

The LLM is shown at most two tools, gated by the `tools` directive:

| | `module_run` capability | `filesystem` capability |
|---|---|---|
| DSL | `tools 'module_run'` | `tools 'filesystem'` |
| LLM-facing tool | one `module_run(module, args)`; `module` is an **enum of the script's `include`d module names**; `args` is a generic object | one `filesystem(path, operation, content?)` |
| Per-module hints | the tool **description** aggregates each available module's one-line summary + input fields (reusing `ModuleMetadataToolSchema`/`ModuleSpecToolSchema`) | n/a |
| `ToolDescriptor.type` | `'module_run'` | `'filesystem'` |
| Execution | parse `{module, args}` → look up the pre-wired bridge for `module` → marshal `args` via the **existing binding rules** (`ProcessEntryHandler.getProcessArguments`) → run as a dataflow node → serialize outputs to JSON | read/write/list/exists confined by `SandboxGuard` to the agent work dir (+ module output dirs for reads) |
| Lives in | core `ModuleToolBridge.callModuleRun` (routes by module name) | core `ModuleToolBridge.callFilesystem` |

Bad `module` names or malformed `args` are returned to the LLM as `{"error":…}` tool
results (the existing dispatch-level-error mechanism), so it can recover and retry.

---

## Milestone 1 — AiServices engine swap (behaviour-preserving)

**Effort: M.** Plugin-only. No DSL/AST/grammar/config/core change.

### Scope
Replace the hand-rolled loop in `LangChainAgentRunner.runWithTools()` with an
`AiServices` proxy. Both the tool path (now AiServices-driven) and the
structured-output single-shot path produce byte-identical observable behaviour. The
tools-XOR-structured-output limitation stays enforced upstream (not lifted here).

### Tasks (file-level)
1. **`plugins/nf-agent/build.gradle`** — the AiServices classes ship in the
   aggregator `langchain4j` artifact, not `langchain4j-open-ai`.
   **Must-fix (version):** the aggregator `langchain4j:1.15.1` is **not published**.
   Bump `langchain4j-open-ai` to `1.16.3` **and** add `dev.langchain4j:langchain4j:1.16.3`
   so `langchain4j-core` resolves to a single version. Verify with
   `./gradlew :plugins:nf-agent:dependencies` that exactly one `langchain4j-core` (1.16.3)
   appears.
2. **New `plugins/nf-agent/src/main/nextflow/agent/AgentService.groovy`** — minimal
   proxy interface `interface AgentService { String chat(String userMessage) }`.
   Named `AgentService` (not `Agent`) to avoid clashing with langchain4j-agentic's
   `@Agent`. No annotations.
3. **`plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy`** —
   rewrite **only** `runWithTools()`:
   - Delete the manual loop (the `AiMessage` accumulator, `for(maxIterations)`,
     `ChatRequest.builder()`, `model.chat(req)`, `hasToolExecutionRequests()`, the
     inner dispatch loop, `ToolExecutionResultMessage.from(...)`, trailing `throw`).
   - Keep `run()`, `runSingleShot()`, `composeMessages()`, `timeoutSeconds()`, the
     `DEFAULT_*` constants, `ModuleToolAdapter.toToolSpecification(...)`, and the
     `ChatModelFactory.createModel(model, timeout, null)` null-schema call (XOR holds).
   - New body: build `Map<ToolSpecification, ToolExecutor>` (one entry per spec),
     each executor `{ ToolExecutionRequest ter, Object memId -> request.dispatch.call(ter.name(), ter.arguments()) } as ToolExecutor`;
     seed memory; build the proxy with `systemMessageProvider`; invoke.
   - Verified-correct 1.16.3 signatures: `ToolExecutor.execute(ToolExecutionRequest, Object)`
     (2-arg, in `dev.langchain4j.service.tool`), `tools(Map<ToolSpecification, ToolExecutor>)`,
     `systemMessageProvider(Function<Object,String>)`, `MessageWindowChatMemory.withMaxMessages(int)`.
   - **Prompt-duplication wiring:** seed **only** the `SystemMessage` into memory and
     pass the full composed user text (prompt + `inputJson`) to `chat()`. Pin with a
     test asserting **exactly one `UserMessage`**.
   - **Must-fix (cap off-by-one + exception type):** `maxSequentialToolsInvocations(N)`
     permits only **N−1** model round-trips and throws a plain `RuntimeException`, not
     the existing `IllegalStateException`. Fix: pass `maxIterations + 1` (reconcile
     against the pinned dispatch-count test) **and** wrap `agent.chat(...)` in
     try/catch to re-throw
     `IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})")`.
4. **`plugins/nf-agent/src/test/.../LangChainAgentToolLoopTest.groovy`** — rewrite.
   **Must-fix (mock strategy):** AiServices routes through `ChatExecutor.execute()` →
   `ChatModel.chat(ChatRequest)`, so the mocked-`ChatModel` injection via
   `AiServices.builder().chatModel(mock)` still works — but **prototype it against the
   real 1.16.3 artifact before finalizing the test list**. Assert on **memory
   contents** (AiServices reuses the seeded `MessageWindowChatMemory`), not a
   hand-built message list.

**Unchanged (verify):** all `modules/nextflow/.../agent/*`, all `nf-lang`, and
`ChatModelFactory`/`ModuleToolAdapter`/`JsonSchemaMapper`.

### Check-in / exit criteria
- `./gradlew :plugins:nf-agent:dependencies` shows exactly one `langchain4j-core` (1.16.3).
- Unit suite green, asserting: final text == model's last `AiMessage.text()`;
  `dispatch` called exactly once for a single-tool scenario with the right args;
  factory called with **null** schema on the tool path; model saw the tool result fed
  back; **exactly one `UserMessage`**; cap-exceeded throws `IllegalStateException` with
  the existing message.
- Boundary: `grep -r dev.langchain4j modules/nextflow/src/main` returns empty.
- **Gated real-LLM run** (`[e2e prod]`, `OPENAI_API_KEY`): `AgentToolEndToEndTest`
  (real multi-step tool loop) and `AgentEndToEndTest` (structured output) both pass.

### Risks
Wrong artifact version → resolved per task 1. Cap semantics → test-pinned + re-wrapped
exception. Mock injection → prototype against the real artifact early.

---

## Milestone 2 — Goal-directed multi-turn loop

**Effort: S.** Additive; orthogonal to Milestone 3.

### Scope
Add one optional `goal` directive folded into the system message to *steer* the
existing AiServices loop toward a declared objective. Deliberately thin: the
multi-turn loop already exists; completion stays "model stops calling tools";
`maxIterations` stays the hard brake. No new runtime primitive, no AST/grammar change,
no `done()` tool, no per-turn completion evaluator (explicitly deferred). Works with
one tool, many tools, or zero tools.

### DSL
```groovy
include { skesa }        from 'nf-core/skesa'
include { assemblyscan } from 'nf-core/assemblyscan'
include { prokka }       from 'nf-core/prokka'

agent triage {
    model 'openai/gpt-5-mini'
    instruction 'You are a careful bioinformatics assistant.'
    goal 'Assemble the reads, QC the assembly, and only annotate if N50 >= 50kb.'
    tools 'module_run'

    input:  isolate: Isolate
    output: report: String
    prompt: "Process isolate ${isolate.sample_id}."
}
```

### Tasks (file-level)
1. **`modules/nf-lang/.../dsl/AgentDsl.java`** — add `void goal(String value);` to the
   directive DSL, mirroring `instruction(String)` exactly. A `String`-valued directive
   identical in shape to `instruction` needs no grammar/`AgentNode`/visitor change — it
   lands in the generic directives map.
2. **`modules/nextflow/.../script/AgentDef.groovy`** — add
   `String getGoal() { directives.get('goal') as String }` next to `getInstruction()`;
   capture it in `run()` and pass it to the request.
3. **`modules/nextflow/.../agent/AgentRunnerRequest.groovy`** — add a `String goal`
   field. **Must-fix (positional churn):** this class is `@Canonical`. **Append `goal`
   at the END** (after `requestTimeoutSeconds`); do not insert it adjacent to
   `instruction`. **Switch the single call site (`AgentDef.run`) to the Groovy
   named-args map constructor** (`new AgentRunnerRequest(model:…, goal:…)`) so
   Milestone 3's further field addition is non-breaking. Update the positional-order
   javadoc in lockstep. `AgentConfig` — no change (`goal` is per-agent intent).
4. **`plugins/nf-agent/.../LangChainAgentRunner.groovy`** — add
   `composeSystemMessage(request)` = instruction + a labelled goal suffix
   (`\n\nGoal:\n${goal}\nYou are done when the goal is met; produce your final answer
   as plain text.`), returning `null` when both absent. Use it for the system text in
   both `runWithTools` (the seeded `SystemMessage`) and `runSingleShot` so
   single-tool/zero-tool agents benefit. No change to loop mechanics, `maxIterations`,
   dispatch, factory, or adapter.

### Check-in / exit criteria
- Unit: `composeSystemMessage` covers instruction-only / goal-only / both (ordering) /
  neither (null); mock model receives goal text in the `SystemMessage`; a multi-turn
  loop preserves dispatch order and returns final text; `maxIterations` brake still
  trips → `IllegalStateException`.
- Core: `getGoal()` returns value/null; a test pins `AgentRunnerRequest` field-by-field
  (guards constructor order); boundary grep stays empty.
- **Gated real-LLM run**: a goal-directed agent (goal requiring 2–3 turns) reasons
  across turns and terminates on goal satisfaction within the cap.

### Risks
Weak models ignore goal guidance → `maxIterations` is the hard brake; document `goal`
as advisory. The agent operator thread is held for the whole (now-longer) conversation
→ pre-existing behaviour; **do not let `goal` silently raise the default
`maxIterations`**.

---

## Milestone 3 — `filesystem` tool, work-dir sandbox & include-driven `module_run`

**Effort: L.** `bash` and `nextflow run` deferred; module execution stays native.

### Scope
- A **per-invocation agent work-dir sandbox**, allocated like a task work dir.
- A generic **`filesystem`** tool — read/write/list/exists — confined to the sandbox
  (writes) and the sandbox + module output dirs (reads), enforced by a new
  `SandboxGuard`.
- A single generic **`module_run`** tool: `module_run(module, args)`, where `module`
  is an enum of the script's `include`d modules and `args` is marshalled per the
  chosen module's existing binding rules.
- `tools` becomes a **capability list** (`'filesystem'`, `'module_run'`). The module
  set is inferred from `include` statements — no module enumeration in `tools`, no
  symbol resolution.

### DSL
```groovy
include { skesa }        from 'nf-core/skesa'
include { assemblyscan } from 'nf-core/assemblyscan'

agent triage {
    model 'openai/gpt-5-mini'
    goal  'Assemble and QC the isolate; write a one-line summary file.'
    tools 'module_run', 'filesystem'   // generic capabilities only

    input:  isolate: Isolate
    output: report: String
    prompt: "Process isolate ${isolate.sample_id}."
}
// The LLM sees: module_run(module: enum[skesa, assemblyscan], args), filesystem(...)
```

### Must-fix items (gating this milestone)
- **Per-record state must not live on the shared bridge.** `ModuleToolBridge` is
  constructed **once pre-ignition** and reused for every record. The work dir and the
  module-output read-whitelist are **per-record** and must be threaded through a
  **per-call dispatch context** (resolved from the in-flight request — the plugin
  already carries `AgentRunnerRequest.agentWorkDir`), keeping the bridge stateless
  across records. Putting them on the shared bridge is a latent isolation leak that
  silently breaks under any future parallelism (`maxForks`, multiple agents,
  `executeToolsConcurrently`).
- **`tools` capability validation.** `tools` entries are now reserved capability
  strings (`'filesystem'`, `'module_run'`). Validate each entry against the known set
  and reject an unknown capability with a clear error (prefer compile time in the
  nf-lang resolver; runtime is acceptable). No symbol/`VariableScopeVisitor` work is
  needed — the earlier symbol-resolution must-fix no longer applies.

### Tasks (file-level)
1. **`modules/nf-lang/.../dsl/AgentDsl.java`** — keep `tools(String... values)` as a
   capability list; update `@Description` to document the reserved capabilities. Add
   validation that each entry is a known capability. No grammar / visitor / lowering
   change.
2. **`modules/nextflow/.../agent/ToolDescriptor.groovy`** — add `String type`
   (default `'module_run'` for module tools / set per generic tool).
3. **`modules/nextflow/.../agent/AgentRunnerRequest.groovy`** — append
   `Path agentWorkDir` (in-JVM only, never serialized; like `dispatch`). Construction
   is already named-args after Milestone 2.
4. **New `modules/nextflow/.../agent/SandboxGuard.groovy`** —
   `isPathInSandbox(candidate, agentWorkDir, readableDirs)`: realpath/normalize, reject
   `..` traversal, reject absolute-outside-whitelist, **and reject symlink-target
   escape** (test the symlink case, not just literal `..`). Module-output reads are
   **read-only** and canonicalized to the specific output file parents actually
   returned — never a broad dir grant.
5. **New `modules/nextflow/.../agent/FilesystemToolSchema.groovy`** — static factory
   returning a `ToolDescriptor` of `type='filesystem'` with schema
   `{path:string, operation:enum[read,write,list,exists], content:string?}`. The work
   dir is bound by the dispatcher per call, never supplied by the LLM.
6. **New `modules/nextflow/.../agent/ModuleRunToolSchema.groovy`** — static factory
   building the single `module_run` `ToolDescriptor` (`type='module_run'`): `module`
   is a string **enum of the included module names**; `args` is a generic object. The
   description aggregates, per available module, its one-line summary + input fields
   (reuse `ModuleMetadataToolSchema`/`ModuleSpecToolSchema` to render compact text, not
   per-tool JSON schemas).
7. **`modules/nextflow/.../script/AgentDef.groovy`** —
   - **Module discovery:** when `'module_run'` is enabled, collect the candidate
     modules from the script's `include` statements via `ScriptMeta` (the resolved
     included `ProcessDef`s / module references), instead of from `tools` strings.
   - **Pre-wire** each candidate module's bridge before ignition (as today), keyed by
     module name; build the `module_run` descriptor with the name enum.
   - When `'filesystem'` is enabled, register the `FilesystemToolSchema` descriptor.
   - In `run()` allocate the per-invocation sandbox under an `agent/`-scoped bucket of
     `session.workDir` (reuse the hashed `work/xx/yy` scheme, keyed by `[name, idx,
     inputJson]`), `mkdirs()`, and pass it into the **per-call dispatch context**.
   - `afterStop`: delete the sandbox only when `session.config.cleanup` is set.
8. **`modules/nextflow/.../agent/ModuleToolBridge.groovy`** — keep `call()`
   `synchronized`. Route on tool `type`:
   - `module_run` → new `callModuleRun(argsJson)`: parse `{module, args}`, resolve the
     pre-wired bridge for `module` (return `{"error":…}` if unknown), marshal `args`
     via `ProcessEntryHandler.getProcessArguments` (the existing rules), run, serialize
     outputs to JSON. Collect that call's serialized output file parents into the
     **per-request** read whitelist.
   - `filesystem` → new `callFilesystem`: read the **per-call context** (sandbox +
     module-output read whitelist), return the same `{...}`/`{"error":…}` JSON contract.
9. **Plugin — no change.** `ModuleToolAdapter` maps any `ToolDescriptor` Map (including
   the `module_run` enum and `filesystem` schema); `LangChainAgentRunner` builds one
   tool map with both capabilities side by side.

### Check-in / exit criteria
- Unit (no network): `SandboxGuardTest` (accept sandbox + module-output paths; reject
  `..`, absolute-outside, **symlink escape**); `ModuleToolBridgeGenericTest`
  (`filesystem` read/write/list/exists in an injected temp dir; `{"error":…}` on
  out-of-sandbox; `module_run` routes to the right module, marshals args via the
  existing rules, serializes outputs; `{"error":…}` on unknown module / malformed
  args); `AgentDefToolResolutionTest` (`tools 'module_run','filesystem'` registers both
  capabilities; `module` enum == the included module names; unknown capability errors).
- Regression: existing examples/tests are **migrated** from `tools 'nf-core/skesa'`
  (module strings) to `include { skesa } …` + `tools 'module_run'`; the migrated
  examples pass.
- **Gated real-LLM run**: the `triage` agent calls `module_run(module:'skesa', …)`,
  reads a module output via `filesystem`, writes a summary file in its sandbox, the
  sandbox exists then is cleaned per `cleanup`, and **no record sees another record's
  sandbox or module-output dir** (isolation assertion).

### Backward compatibility / migration
The legacy module reference forms in `tools` — registry strings (`tools 'nf-core/…'`),
local-path strings, and in-scope process names (`tools 'uppercase'`) — are **removed**
in favour of `include` + `'module_run'`. Existing examples (`isolate-triage`, `skesa`,
`tool`, `two-agents`) and the docs are migrated. (The feature is experimental; this is
an accepted breaking change to its surface.)

### Deferred (explicit, out of Milestone 3)
- **`bash`/script tool** — ungated host arbitrary command execution. A path sandbox is
  not a command sandbox. When taken up, ship behind **double opt-in**
  (`agent.allowShellTools` config flag, default false, **plus** `tools 'bash'`) and
  **containerized execution** (run as an executor task in the agent work dir).
- **`nextflow run` (CLI / nested pipeline launch)** — nested-session lifecycle hazard.
  Native `module_run` already covers "run a module"; a CLI-style pipeline-launch tool
  is a separate design.

### Risks
Serialized `synchronized` dispatch also serializes `filesystem` ops → document them as
quick ops; if a generic tool can block indefinitely, give it a separate non-serialized
path. Remote work dir → `filesystem` is `file://`-only; detect and return
`{"error":…}`. Pre-wiring **every** included module even when unused → acceptable for a
curated include set; revisit if include counts grow large.

---

## Independence & sequencing

- **Milestone 1** is self-contained (one plugin, no DSL/core change), lowest risk, and
  ships/releases on its own — every existing `.nf` script behaves identically.
- **Milestone 2** depends on M1's AiServices builder for the tool path (it edits the
  `systemMessageProvider`/memory seeding). The `goal` directive, `getGoal()`, and the
  `AgentRunnerRequest.goal` field are purely additive.
- **Milestone 3** is orthogonal to M2 (`goal` does not touch tool wiring) and only
  lightly coupled to M1 (it adds generic descriptors to the same tool map the
  AiServices builder consumes). Within M3, the include-discovery + `module_run` work
  and the sandbox/`SandboxGuard`/`filesystem` work are largely independent.

**Shared seam:** `AgentRunnerRequest` grows a field in both M2 (`goal`) and M3
(`agentWorkDir`). Adopt named-args construction at the single call site in whichever
lands first; append new fields; pin field order with a unit test.

**Parallelizable:** M2 and M3 can be developed concurrently after M1 merges.

## Resolved defaults

- `tools` = generic capabilities only (`'filesystem'`, `'module_run'`); unknown
  capability → error.
- Modules inferred from `include` statements; `module_run` is a single tool with a
  `module` enum + generic `args`; per-module hints in the tool description.
- The legacy module-reference string forms in `tools` are removed; examples migrated.
- `maxIterations` default **unchanged** under `goal`.
- `filesystem` is read **and** write (write confined to the agent work dir).
  `bash`, `nextflow run` deferred.

## Open questions (none blocking)

- Whether in-scope processes (defined in the same script, not `include`d) should also
  be `module_run` candidates, or strictly includes-only (current default:
  includes-only; the `tool`/`uppercase` toy example migrates to a module file).
- Final `filesystem` schema details (`write` content size cap; whether `list` returns
  names or stat records).
- How richly to render per-module hints in the `module_run` description (one-line vs
  full input-field listing) — a prompt-quality tuning question.

## Testing strategy (cross-cutting)

- **Boundary test in CI**: `grep -r dev.langchain4j modules/` must be empty.
- **Concurrency invariant test**: assert tool calls execute on the calling thread
  (AiServices `executor` stays null).
- **Constructor-order test**: pin `AgentRunnerRequest` field-by-field.
- **Gated real-LLM E2E** per milestone (`@Requires(OPENAI_API_KEY)`, `[e2e prod]`),
  skipped keyless.
