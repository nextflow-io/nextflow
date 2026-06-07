# Nextflow Agent — Design & Implementation Guide

**Branch:** `agent-syntax-and-model`  ·  **Status:** v1 (OpenAI-only)
**Companion design spec:** `docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md`

This is the engineering reference for the Nextflow **agent** primitive: how an `agent`
flows from DSL syntax to a running dataflow operator, how it exposes modules/processes
as LLM **tools** that execute as real dataflow nodes, and the critical implementation
choices that make synchronous tool calls work inside Nextflow's asynchronous runtime.
It is written to serve as a manual for review and future decisions — every non-obvious
choice carries its *why*.

> **Scope.** The whole agent feature, with emphasis on the **module-as-tool** execution
> path. File/line references are approximate — anchor on method names, not line numbers.
> Read alongside the design spec: the spec states intent and the v1 boundary; this guide
> states mechanism and rationale.

---

## 1. What an agent is

```groovy
nextflow.enable.types = true

agent triage {
    model 'openai/gpt-5-mini'
    instruction '''high-level goal; the LLM decides which tools to call and when'''
    tools 'nf-core/skesa', 'nf-core/assemblyscan', 'nf-core/prokka'   // module tools
    input:  isolate: Isolate          // a typed record (the agent's single input)
    output: verdict: String           // a plain value when tools are used
    prompt: """Triage '${isolate.sample_id}'. Reads: ${isolate.reads}"""
}
```

An `agent` is a **dataflow operator like a process**: it consumes one input channel and
emits one output channel. For each input item it runs **one LLM conversation**; inside that
conversation the LLM may call the declared `tools` any number of times. Each tool call runs
the corresponding **module as a real Nextflow task** (its container, work dir, caching) and
returns a result the LLM reads. Composition with the rest of a pipeline is via channels —
agents do not call each other directly (v1).

**Two output modes** (mutually exclusive in v1):

- **Tools mode** — when `tools` are declared, the runner drives a tool-call loop and the
  final answer is **free text** bound to the output channel verbatim.
- **Structured mode** — when the output is a **record type** *and no tools are declared*,
  the model is constrained to a JSON schema (`RecordSchema`), and the JSON is parsed into a
  typed record. Declaring tools *and* a record output is rejected up front.

---

## 2. Architecture at a glance

```
DSL ─ parse ─► AgentNode ─ lower ─► agent('triage', { … }) ─ run ─► AgentBuilder.build() ─► AgentDef
                                                                                              │
                                       AgentDef.run()  (BUILD TIME, before ignition)          │
                                       ├─ createToolBridge(): resolve+compile tools           │
                                       ├─ ModuleToolBridge: PRE-WIRE queues, capture outputs  │
                                       └─ applyAgentOperator(): one MapOp over the input      │
                                                                                              ▼
input record ──► mapper ──► prompt + inputJson ──► AgentRunnerRequest ──► AgentRunner (plugin: LLM loop)
                                                                                   │
                          ┌────────────────────────────────────────────────────────┘
                          ▼ (per tool call, on the agent operator thread)
        LLM emits tool name + flattened JSON args
                          │  ModuleToolBridge.call(name, argsJson)   [SYNCHRONIZED]
                          ▼
        ProcessEntryHandler.getProcessArguments ──► [ [meta], file(reads) ]  (channel value)
                          │  bind onto pre-wired input queue(s)
                          ▼
        REAL module task runs on the executor (container / work dir / cache)
                          │  read pre-captured output channels (.val); skip eval/topic
                          ▼
        serialize outputs ──► JSON result   (small structured file ⇒ INLINED contents;
                          │                   bulk/binary/oversized file ⇒ absolute-path HANDLE)
                          ▼
        result fed back to the LLM as a ToolExecutionResultMessage
                          │  (LLM chains a handle into the next tool, or reasons over inlined data)
                          ▼
        LLM returns a final text answer ──► AgentDef emits it on the output channel
                          │
        agent input exhausted ──► bridge.close() poisons input queues ──► operators terminate
```

**Layering.** The DSL keyword/AST and *all* dataflow + tool-bridge logic live in **core**
(`modules/nextflow`, `modules/nf-lang`). Only the langchain4j client — building the
`ChatModel` and running the chat/tool loop — lives in the **`nf-agent` plugin**, reached
through the `AgentRunner` SPI. Core never depends on langchain4j; it passes a portable
`ToolDescriptor` (plain `Map` schemas) and a `ToolDispatcher` callback across the boundary.

---

## 3. Component & file map

### DSL front-end (`modules/nf-lang`)

| Concern | File |
|---|---|
| `AGENT` token | `src/main/antlr/ScriptLexer.g4` |
| `agentDef` / `agentBody` / `agentDirectives` / `agentInputs` / `agentOutputs` / `agentPrompt` rules (I/O reuse `processInput`/`processOutput`) | `src/main/antlr/ScriptParser.g4` |
| AST node | `src/main/java/nextflow/script/ast/AgentNode.java` |
| AST build | `src/main/java/nextflow/script/parser/ScriptAstBuilder.java` (`agentDef`, `agentInputs`, …) |
| Lowering to `agent('name', { … })` | `src/main/java/nextflow/script/control/AgentToGroovyVisitorV2.java` |
| Resolve/scope/visitor support | `ScriptResolveVisitor`, `VariableScope*`, `ScriptVisitor*` |

### Core runtime (`modules/nextflow/src/main/groovy/nextflow`)

| Concern | File |
|---|---|
| Runtime delegate that captures the lowered agent body, builds `AgentDef` | `script/AgentBuilder.groovy` |
| Agent runtime model + operator + tool resolution/compilation | `script/AgentDef.groovy` |
| LLM-args → channel-values marshalling (shared with `module run`) | `script/ProcessEntryHandler.groovy` |
| **Tool execution bridge**: pre-wire, marshal, dispatch, serialize, poison | `agent/ModuleToolBridge.groovy` |
| Dispatch callback interface (SAM) | `agent/ToolDispatcher.groovy` |
| Portable, langchain4j-free tool descriptor | `agent/ToolDescriptor.groovy` |
| Output handle-vs-inline decision | `agent/ToolOutputReader.groovy` |
| `agent` config scope | `agent/AgentConfig.groovy` |
| Runner SPI + plugin lookup (test seam) | `agent/AgentRunner.groovy`, `agent/AgentRunnerProvider.groovy` |
| Immutable request DTO to the runner | `agent/AgentRunnerRequest.groovy` |
| Tool descriptor/schema from registry metadata | `agent/ModuleMetadataToolSchema.groovy` |
| Tool descriptor/schema from sibling `meta.yml` `ModuleSpec` | `agent/ModuleSpecToolSchema.groovy` |
| Tool schema from a typed in-scope process | `agent/ProcessToolSchema.groovy` |
| Structured-output schema by reflecting a record type | `agent/RecordSchema.groovy` |

### `nf-agent` plugin (`plugins/nf-agent/src/main/nextflow/agent`)

| Concern | File |
|---|---|
| pf4j plugin entry | `AgentPlugin.groovy` |
| `AgentRunner` impl: single-shot / structured / tool-loop | `LangChainAgentRunner.groovy` |
| `provider/model` → langchain4j `ChatModel` (**OpenAI only, v1**) | `ChatModelFactory.groovy` |
| portable `ToolDescriptor` → langchain4j `ToolSpecification` | `ModuleToolAdapter.groovy` |
| portable schema `Map` → langchain4j `JsonSchema` | `JsonSchemaMapper.groovy` |
| Dep: `dev.langchain4j:langchain4j-open-ai:1.15.1` | `build.gradle` |

---

## 4. From syntax to a running operator

1. **Parse / lower.** The grammar recognizes `agent name { … }`; `ScriptAstBuilder` builds an
   `AgentNode` (directives statement block, typed input/output params, a prompt statement);
   `AgentToGroovyVisitorV2` lowers it to a runtime call `agent('triage', { … })` whose closure
   carries the directives, `_input_`/`_output_` declarations and a `PromptDef`.
2. **Capture.** Running that closure drives an `AgentBuilder`: directives are captured via
   `methodMissing` against the fixed set `[model, instruction, tools, maxIterations]` (unknown
   ⇒ error); `_input_`/`_output_` record typed declarations; `withPrompt` stores the `PromptDef`.
   `build()` produces the `AgentDef`.
3. **Wire (build time).** `AgentDef.run(args)` runs in the workflow body **before ignition**.
   It validates exactly-one-input / exactly-one-output, resolves the `agent` config defaults,
   decides tools-vs-structured, builds the tool bridge (§6), and applies the operator (§5).

`AgentDef.run` effective-value precedence: the agent's own directive first, then the `agent`
config scope, then a built-in default (`maxIterations` 20, `requestTimeout` 120 s,
`maxToolOutputInlineSize` 32 KB). The model is `directive ?: agent.defaultModel` (a null model
surfaces downstream in the runner).

---

## 5. The agent operator

`applyAgentOperator` builds a single `MapOp`-style operator over the input channel (mirrors
`nextflow.extension.MapOp`) with one addition: a `DataflowEventAdapter` whose
- `afterStop` calls `bridge.close()` (poison the tool queues) once the input is exhausted, and
- `onException` *also* poisons the queues and aborts the session, so the network unwinds on failure.

The per-item **mapper** closure:
```groovy
cl.setDelegate([(inputName): item]); final promptText = cl.call()?.toString()
final inputJson = toJson(item)                       // input record as JSON (Path ⇒ abs string)
final req = new AgentRunnerRequest(model, instruction, promptText, maxIter,
                                   tools, outputSchema, inputJson, toolSpecs, bridge, timeoutSecs)
final result = runner.run(req)
// structured mode only: stripFences(result) → JsonSlurper → TypeHelper.asRecordType(map, outputClass)
```
The input record is rendered into the `prompt` **and** appended verbatim as JSON (the runner adds
`Input (JSON): …` to the user message). The agent's input shape and a tool's input shape are
**independent** — the LLM bridges them.

`AgentRunnerProvider.get()` resolves the active `AgentRunner` from loaded plugins (priority
extensions), with a package-scope `testRunner` seam for unit tests; absent ⇒ abort with a hint to
enable `nf-agent`.

---

## 6. Wiring the tools (build time, **before** the dataflow network ignites)

Tools are resolved and pre-wired in the workflow body, before ignition. **This is the central
design decision.** Invoking a process synchronously *after* ignition deadlocks GPars; pre-wiring
plus a blocking pull on the captured outputs is what makes synchronous tool calls work.

### 6a. Resolution — `AgentDef.createToolBridge()`

Each `tools` entry (a String) is resolved in order:

1. **In-scope process name** (`ScriptMeta.getProcess`) → used directly → **scalar mode**.
2. **Local module path** (`looksLikeModulePath`: ends `.nf`, or starts `./` `../` `/`) →
   `resolveModuleTool` → `compileModuleProcess` compiles the module's single process to a
   `ProcessDef`. A sibling `meta.yml`/`meta.yaml` (if present) → `ModuleSpec` → **spec mode**.
3. **Registry ref** (`looksLikeRegistryRef`: parses as `ModuleReference`, optional `@version`) →
   `resolveRegistryTool`: `ModuleResolver.resolve(ref, version, autoInstall=true)` (local install
   first, else download) → `compileModuleProcess`. The **public registry `ModuleMetadata`** is
   fetched via the same `RegistryClient` (`getModule(name).latest.metadata`, or
   `getModuleRelease(name, version).metadata` when pinned) for the descriptor; the sibling
   `meta.yml` spec still drives marshalling. Tool name is sanitized (`/`,`-`,`@`,`.` → `_`).
4. Otherwise → error (`tool … is not a process in scope`).

> **Critical (per-tool module isolation).** `compileModuleProcess` gives each tool module its
> **own** `ScriptBinding` (seeded with the owner's params) via `ScriptLoaderV2.setMainBinding`,
> mirroring `IncludeDef.loadModuleV2`, and sets `setModule(true)` so the standalone process does
> not auto-build an entry workflow. Without the isolated binding, all tool modules shared
> `session.binding` and `BaseScript.setup()`'s `moduleDir` write was last-wins — so every tool's
> `conda "${moduleDir}/environment.yml"` / `container` resolved to the **last-compiled** module
> → one shared container for all tools. The resolved `modPath` is remembered (keyed by the
> `ProcessDef`) so the sibling `meta.yml` can be located.

### 6b. Pre-wire — `ModuleToolBridge.wireSpec()` (once per spec-driven tool)

```groovy
// one persistent queue per declared input CHANNEL (== ProcessEntryHandler's arg count)
final int nInputs = declaredInputChannelCount(proc)      // V2 input params, or V1 inputs
final toolIns = (0..<nInputs).collect { CH.queue() }
final cloned  = proc.cloneWithName(name)
final ChannelOut out = cloned.run(toolIns as Object[])    // pre-wire into the network

// CRITICAL: capture the output READ channels NOW, at wiring time
final outReadChannels = (0..<out.size()).collect { CH.getReadChannel(out[it]) }
```

The number of queues is the **process's** declared input-channel count (what
`getProcessArguments` returns), not the spec's input count; a mismatch logs a drift warning and
uses the process count. The **descriptor** (name + description + flattened input JSON-schema) is
built from the registry `ModuleMetadata` when present (`ModuleMetadataToolSchema`), else the
sibling `meta.yml` `ModuleSpec` (`ModuleSpecToolSchema`). When both exist, `warnOnInputDrift`
warns if their flattened input names differ.

`wireScalar()` is the simpler in-scope-process path: one queue per typed scalar input, a single
scalar output, schema from `ProcessToolSchema` (scalar types only; tuple/path/map/record fail loud).

> **Critical (read-channel timing).** Output read channels are captured at **wiring**, not lazily
> at dispatch. These process outputs are broadcast-style — an adapter created *after* the task
> binds its value never receives it. Subscribing at wiring (before any value binds) lets `callSpec`
> read every output regardless of order; lazily subscribing per-read deadlocked the 2nd+ output of
> multi-output tools (e.g. prokka).

---

## 7. Agent input → specific module input (marshalling)

The LLM emits a tool call with a **flattened JSON** argument object shaped by the tool's input
schema, e.g. `nf_core_skesa({"meta":{"id":"isolate_001"}, "fastq":"…/sample.fastq"})`. It maps
`isolate.sample_id → meta.id` and `isolate.reads → fastq` itself; no field-by-field correspondence
between the agent input and a tool input is required.

`ModuleToolBridge.callSpec()` reuses the **exact same binding logic as `nextflow module run`**:
```groovy
final List args = ProcessEntryHandler.getProcessArguments(tool.processDef, parsed, spec)
for( int i=0; i<tool.toolIns.size(); i++ )
    tool.toolIns[i].bind(args[i])           // bind onto the pre-wired queues → task runs
```
`getProcessArguments` (→ `getProcessArgumentsV1`) treats the flattened arg map as the params map
(dot-notation folded into nested maps), looks up each declared input by name, **type-coerces** from
the sibling `meta.yml` (`file`→`Nextflow.file`, `map`→Map, `integer`→Integer, …), and assembles a
tuple input in order into a channel value such as `[[id:'s1'], file(reads)]`.

> **Critical (optional path inputs).** `getValueForInputV1` treats a **null OR blank-string** file
> argument as *not provided* → empty list `[]` (the nf-core convention for optional `path` inputs).
> This avoids `file("")` throwing when the LLM passes `""` for an optional input it doesn't have
> (e.g. prokka's `proteins`/`prodigal_tf`).

A dispatch-level failure (unknown tool, unparseable JSON, bad argument) is **returned as a
`{"error": "..."}` tool result** rather than thrown (`ModuleToolBridge.call`), so the LLM can
correct itself and retry instead of aborting the run.

---

## 8. The LLM tool-calling loop (`nf-agent` plugin)

`LangChainAgentRunner.run` branches on `request.toolSpecs`:

- **No tools** → `runSingleShot`: a single `model.chat(messages)`. When an output record schema is
  present it is compiled to a langchain4j `JsonSchema` (`JsonSchemaMapper`) and the OpenAI model is
  built with a **strict JSON-schema response format** (`ChatModelFactory` sets `responseFormat` +
  `strictJsonSchema(true)`).
- **Tools** → `runWithTools`, the manual loop (tools and structured output are mutually exclusive —
  no `responseFormat` is forced):

```groovy
final specs = request.toolSpecs.collect { ModuleToolAdapter.toToolSpecification(it) }
for( int i = 0; i < maxIterations; i++ ) {
    final response = model.chat(ChatRequest.builder().messages(messages)
                                   .toolSpecifications(specs).build())
    final ai = response.aiMessage()
    if( !ai.hasToolExecutionRequests() )
        return ai.text()                                   // final answer
    messages.add(ai)                                       // assistant turn (the tool requests)
    for( ToolExecutionRequest ter : ai.toolExecutionRequests() ) {
        final String result = request.dispatch.call(ter.name(), ter.arguments())   // ← the bridge
        messages.add(ToolExecutionResultMessage.from(ter, result))                 // feed result back
    }
}
throw new IllegalStateException("Agent exceeded the maximum number of tool-call iterations (...)")
```

- Tool specs are advertised on **every** chat request.
- `request.dispatch` is the `ModuleToolBridge` (a `ToolDispatcher`); `call(...)` is
  **`synchronized`**: the bridge submits and awaits **one tool call at a time** so input/output
  stay correlated over each tool's single shared pre-wired instance. This serializes *dispatch*,
  **not** module execution — the module runs as a normal Nextflow task on its configured executor
  (its own thread/process, or a remote grid/cloud backend) while the dispatcher thread blocks on
  `.val` (see §6 and §12). Several tool requests in one assistant turn are therefore submitted
  sequentially, not overlapped.
- The loop ends on a text-only assistant message, or throws on hitting `maxIterations` (the operator
  then aborts the session).
- `composeMessages` builds: optional `SystemMessage(instruction)` + `UserMessage(prompt + "\n\nInput (JSON):\n" + inputJson)`.

`ChatModelFactory` (v1) supports **only** `openai`; any other provider, or a missing
`OPENAI_API_KEY`, throws here (at run time, on the first input record).

---

## 9. Tool output → result, and chaining into downstream tools

After the task completes, `callSpec` reads outputs in declared order and serializes them:
```groovy
final outParams = declaredOutputParams(tool.processDef)   // authoritative param list
for( int i=0; i<outputs.size(); i++ ) {
    final param = outputs[i]
    if( isEvalOutput(param) || isTopicOutput(outParams, i) ) continue   // skip versions/topic
    if( i >= nOut ) continue
    final ch = tool.outReadChannels[i]                    // read channel captured at WIRING
    final value = ch.val                                  // blocking pull
    result.put(key, serializeOutput(param, value))
}
return JsonOutput.toJson(result)
```

> **Critical (skip non-binding bookkeeping outputs).** An nf-core `versions` output is routed to a
> `topic:` (a topic-source channel that never binds a readable per-task value). It is detected
> **authoritatively from the ProcessDef** — `isTopicOutput` checks `OutParam.getChannelTopicName()`,
> guarded to `BaseOutParam` since typed/V2 `ProcessOutput.getChannelTopicName()` throws — *not* from
> the `meta.yml` `type`, which the registry types inconsistently (`eval` vs `string`, e.g.
> assemblyscan vs skesa). Reading such a channel's `.val` would block the dispatch forever.

### Handle vs inlined — `ToolOutputReader.readOrHandle(path, maxBytes)`

`serializeOutput` turns a tuple value into a record keyed by component name; each file-typed value
goes through the reader:

- **Bulk / binary / oversized** → **absolute path string** (an opaque handle). Decision: extension
  *not* in `TEXT_EXTENSIONS` (`json,tsv,csv,txt,tab,yaml,yml,log,md`), or size `> maxBytes`, or a NUL
  byte in the first 8 KB. The LLM never sees the bytes; it forwards the handle to the next tool.
  (Oversized inline candidates return `[path: <abs>, note: "content not inlined: … exceeds … limit"]`.)
- **Small structured text** (text extension, under the cap, non-binary) → **contents inlined** as a
  UTF-8 string so the LLM can reason over them.

The cap is `agent.maxToolOutputInlineSize` (default 32 KB), wired into the bridge as `maxInlineBytes`.

This split enables both **downstream chaining** and **data-driven control flow** — e.g. in
`examples/agents/isolate-triage`:
1. `nf_core_skesa` → `fasta` is a `.fa` → **path handle**.
2. The LLM passes that handle as `nf_core_assemblyscan`'s `assembly` input (chaining).
3. `nf_core_assemblyscan` → stats `.tsv`/`.json` small/structured → **inlined**; the LLM reads N50 /
   #contigs and applies the QC gate.
4. On PASS it calls `nf_core_prokka` (forwarding the contigs handle) and reports the annotation path;
   on FAIL it stops — prokka is never invoked.

The serialized JSON becomes the `ToolExecutionResultMessage` content the LLM reads next iteration.

---

## 10. Completion & termination

The operator runs one LLM conversation per input item. When the input source is exhausted, the
operator's `afterStop` calls `ModuleToolBridge.close()`, which binds a `PoisonPill` to every
pre-wired tool input queue (idempotent). The tool process operators receive the poison, terminate,
and propagate poison to their outputs (`TaskProcessor.sendPoisonPill`), so the dataflow network
unwinds and the session completes cleanly.

---

## 11. Configuration reference

```groovy
agent {
    defaultModel            = 'openai/gpt-5-mini'   // when an agent declares no `model`
    maxIterationsDefault    = 20                    // when an agent declares no `maxIterations`
    requestTimeout          = '120s'                // per-request LLM chat timeout
    maxToolOutputInlineSize = '32 KB'               // inline cap for structured tool outputs
}
```
All values are **defaults**: an agent's own directive always wins (`AgentConfig` doc: "the scope
only fills in a value the agent did not declare"). Credentials: `OPENAI_API_KEY` (only provider in
v1). Enable the runner with the `nf-agent` plugin.

---

## 12. Design decisions & rationale (quick index)

| Decision | Where | Why |
|---|---|---|
| **Pre-wire tools before ignition**; bind+`.val` at dispatch | `ModuleToolBridge.wireSpec/wireScalar` | synchronous process invocation *after* ignition deadlocks GPars |
| **Capture output read channels at wiring** | `wireSpec` | broadcast outputs: a subscriber created after a value binds misses it → deadlock on multi-output tools |
| **Serialized dispatch** (`synchronized`) — submit-and-await one call at a time | `ModuleToolBridge.call` | single pre-wired instance per tool — keep input/output correlated. Serializes *dispatch*, not the module task (which runs async on its executor) |
| **Registry `ModuleMetadata` first, `meta.yml` fallback** | `createToolBridge`, `ModuleMetadataToolSchema` | richer self-describing tools (descriptions/patterns/enums/`meta.id`); works offline via spec |
| **Marshalling shared with `module run`** | `ProcessEntryHandler.getProcessArguments` | one code path for flattened-JSON → channel values; blank/optional path → `[]` |
| **Inline small structured outputs; handle bulk** | `ToolOutputReader.readOrHandle` | LLM must reason over stats (gate) but only chain bulk artifacts (handles) |
| **Per-tool `ScriptBinding` isolation** | `AgentDef.compileModuleProcess` + `ScriptLoaderV2.setMainBinding` | otherwise `moduleDir` is last-wins → all tools share one container |
| **Skip `eval`/`topic` outputs from the ProcessDef** | `isTopicOutput`/`isEvalOutput` (guarded to `BaseOutParam`) | topic-source channels never bind → `.val` hangs; `meta.yml` type is unreliable |
| **Dispatch errors as `{"error":…}`** | `ModuleToolBridge.call` | LLM can recover from bad args without aborting the run |
| **langchain4j only in the plugin; portable `Map` schemas in core** | `ToolDescriptor`, `ModuleToolAdapter`, `JsonSchemaMapper` | keep core free of any LLM-client dependency |
| **Tools XOR structured output** | `AgentDef.run` guard | the tool loop produces free text; a forced response schema would fail to parse |

---

## 13. Bugs fixed in this work

1. **Inline small structured outputs** (`ToolOutputReader` + `agent.maxToolOutputInlineSize`) —
   surface small `.json`/`.tsv` tool outputs to the LLM; bulk data stays a handle.
2. **Per-tool `moduleDir` isolation** — each tool module gets its own `ScriptBinding` (else all
   tools share one container).
3. **Skip topic-routed `versions` outputs via the ProcessDef** — the `meta.yml` `type` is unreliable;
   reading a topic-source `.val` hangs.
4. **V2 guard** — `isTopicOutput` guards `BaseOutParam` (typed `ProcessOutput.getChannelTopicName()`
   throws).
5. **Blank/optional file inputs** — a `""` path arg → not provided (`[]`), never `file("")`.
6. **Capture output read channels at wiring time** — fixes the multi-output dispatch deadlock.

---

## 14. Known limitations (open)

1. **Genuinely-absent optional output.** A tool that produces *no* file for an `optional: true`
   output is still unhandled: `TaskProcessor.bindOutputsV1` skips binding it, so the channel never
   gets a value for that call and the per-call `.val` blocks until process termination (which the
   blocked dispatch prevents). The bundled example tools always produce their declared outputs, so
   this does not arise there. A robust fix needs either a non-terminating "empty for this task"
   sentinel on the output channel in core, or a bridge-local guard — a core dataflow-semantics
   decision deferred pending direction.
2. **Tool-task failure is not recoverable.** Only *dispatch-level* errors become `{"error":…}` tool
   results. A failure of the underlying module *task* surfaces through the process's own dataflow
   operator and aborts the session; intercepting it as a recoverable tool result fights the runtime
   and is deferred.
3. **OpenAI only.** `ChatModelFactory` rejects non-`openai` providers; only `langchain4j-open-ai` is
   on the classpath.
4. **One input, one output.** Multiple/zero inputs or outputs are rejected; `path` declarations at
   the agent boundary are out of scope (record fields carry path *strings*).
5. **No agent operator directives.** The agent operator is not a `TaskProcessor`, so
   `errorStrategy`/`maxRetries`/`cache` apply to the *tools* (real tasks), not the agent loop.
