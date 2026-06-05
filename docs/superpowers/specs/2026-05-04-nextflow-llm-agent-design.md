# Nextflow LLM Agent — Design

**Date:** 2026-05-04 (revised to track the implementation on branch `agent-syntax-and-model`)
**Status:** Implemented (v1) — see `docs/superpowers/agent-module-tools-implementation.md` for the implementation guide
**Author:** Paolo Di Tommaso

> This document is the original design intent, revised so that every claim matches the
> code as built. Where the implementation deliberately departed from the first draft, the
> change is called out inline and summarized in §10 *Divergences from the original design*.
> The authoritative, code-level reference is the companion **agent-module-tools-implementation.md**.

## 1. Goal & non-goals

### Goal

Add a new top-level `agent` construct to the Nextflow DSL — a process-shaped primitive that wraps an LLM-driven tool-calling loop. Each invocation of an agent receives one input record, runs an LLM with access to a declared set of Nextflow modules as tools, and emits one output record. Agents compose with processes and other agents through the existing channel/workflow model.

The defining property of the design: **a tool call from the LLM is not a function call — it is a real Nextflow module invocation.** The args are wrapped as channel values, the module runs through the standard dataflow runtime (executor, container, retries, work dir, caching), and the resulting output channel record is serialized back to the LLM. The harness is the adapter between the LLM tool-call protocol and Nextflow's dataflow model.

### Non-goals (v1)

- Long-lived conversational agents (one input → one output is the contract)
- Channel-aware orchestrator agents (agents do not subscribe to channels mid-loop)
- Multi-provider LLM support — **v1 ships OpenAI only** (langchain4j is the abstraction, but only `langchain4j-open-ai` is on the classpath and `ChatModelFactory` rejects any other provider)
- Agent-to-agent invocation as tools (composition is via channels only)
- Cost tracking, prompt analytics, fine-grained token budgets
- RAG, persistent memory, conversation state across invocations
- **Concurrent** tool dispatch — multiple tool calls in one assistant turn run **serially** (see §3)
- Streaming partial outputs (langchain4j supports it; surfaced as a follow-up)

> Note: *typed/structured LLM outputs* were listed as a non-goal in the first draft but are
> **implemented in v1** for the no-tools path (see §2 *I/O contract* and §9). They are mutually
> exclusive with tools.

## 2. DSL surface

A new top-level `agent` definition, parsed by `nf-lang` as a sibling of `processDef` and `workflowDef`. Tools are referenced as **modules** (registry coordinate, local include, or an in-scope process), so they bring their own description and I/O metadata.

Agent inputs and outputs use the **typed DSL** (`processInput` / `processOutput` grammar rules), i.e. `name: Type` declarations — `nextflow.enable.types` must be on. Exactly one input and exactly one output are supported in v1.

```groovy
nextflow.enable.types = true

record Isolate {
    sample_id: String
    organism:  String
    reads:     String   // path string — a handle the tools chain
}

agent triage {
    model 'openai/gpt-5-mini'
    instruction 'You triage bacterial isolate assemblies step by step.'
    tools 'nf-core/skesa', 'nf-core/assemblyscan', 'nf-core/prokka'
    maxIterations 20

    input:
        isolate: Isolate

    output:
        verdict: String         // tools present ⇒ plain (non-record) output: the LLM's final text

    prompt:
    """
    Triage isolate '${isolate.sample_id}' (${isolate.organism}).
    Short reads: ${isolate.reads}
    """
}

workflow {
    triage(channel.of(
        record(sample_id: 'isolate_001', organism: 'Escherichia coli', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { v -> "TRIAGE: ${v}" }
}
```

### Directives

| Directive | Required | Meaning |
|---|---|---|
| `model` | yes (or `agent.defaultModel`) | Model identifier in `provider/model` form. **v1: provider must be `openai`** (e.g. `openai/gpt-5-mini`). The provider prefix selects the langchain4j chat-model factory. |
| `instruction` | no | System prompt (agent role/persona) |
| `tools` | no | Varargs list of module references: an in-scope process name, a local `.nf` path, or a registry coordinate (`nf-core/fastqc`, optionally `@version`) |
| `maxIterations` | no | LLM-loop iteration cap (default `agent.maxIterationsDefault`, else 20) |
| `input:` / `output:` | yes | Typed (`name: Type`) channel I/O — exactly one of each in v1 |
| `prompt:` | yes | Templated user prompt with `${var}` interpolation over the single input binding |

The recognized directive set is fixed (`model`, `instruction`, `tools`, `maxIterations`) — an unknown directive fails at build time (`AgentBuilder.checkName`).

### I/O contract (v1)

- **Input:** exactly one typed declaration bound into the `prompt:` template via `${name}`. The whole input record is **also appended to the user message as JSON** (`Input (JSON): …`) so the model sees the structured value even if the prompt omits a field. `Path` values are normalized to absolute strings in that JSON.
- **Output:** exactly one typed declaration.
  - When the output type is a **record type** *and no tools are declared*, structured output is enabled: the runner constrains the model to a JSON schema derived from the record (`RecordSchema`), parses the JSON and binds a typed record. (Code fences around the JSON are stripped first.)
  - Otherwise (a plain type such as `String`, or whenever tools are declared) the LLM's final assistant text is emitted verbatim.
- **Tools XOR structured output:** declaring `tools` together with a record output is rejected at agent-run time with a clear message — the tool loop drives the conversation to free text and cannot also honor a forced response schema.

Standard process directives (`errorStrategy`, `maxRetries`, `time`, `cache`) do **not** apply to the agent operator itself in v1 — an agent is a `MapOp`-style dataflow operator, not a `TaskProcessor`. The *tools* it invokes are real tasks and carry their own directives/retries.

## 3. Architecture

The runtime is **in-JVM**. The DSL keyword/AST and all dataflow/tool-bridge logic live in **core**; only the langchain4j client (model construction + the chat/tool loop) lives in the **`nf-agent` plugin**, reached through the `AgentRunner` SPI. No subprocesses, no external binaries, no IPC.

```
┌─────────────────────────────── Nextflow JVM ─────────────────────────────────┐
│                                                                              │
│  workflow ──► AgentDef.run()  (build time, BEFORE ignition)                  │
│        │                                                                     │
│        ├── createToolBridge(): resolve each `tools` entry                    │
│        │     • in-scope process            → scalar mode                     │
│        │     • local .nf path / registry   → compile to ProcessDef           │
│        │       ref (nf-core/…)               + load sibling meta.yml spec    │
│        │       (+ fetch public registry ModuleMetadata for the descriptor)   │
│        │                                                                     │
│        ├── ModuleToolBridge PRE-WIRES every tool into the dataflow network:  │
│        │     • one persistent CH.queue() per declared input channel          │
│        │     • clone+run the ProcessDef over those queues                    │
│        │     • CAPTURE the output read channels NOW (before any value binds) │
│        │                                                                     │
│        └── applyAgentOperator(source, mapper): one MapOp over the input      │
│              channel; afterStop → bridge.close() poisons the tool queues     │
│                                                                              │
│  per input record (agent operator thread):                                   │
│        render prompt + input-as-JSON → AgentRunnerRequest → runner.run()     │
│                                                                              │
│  runner (nf-agent plugin, LangChainAgentRunner):                             │
│        build OpenAI ChatModel; advertise tool specs on EVERY chat request;   │
│        loop ≤ maxIterations: chat → if tool_calls, dispatch each → append    │
│        results → repeat; plain assistant message ⇒ final answer              │
│                                                                              │
│  on each tool call (back in core, ModuleToolBridge.call — SYNCHRONIZED):     │
│        1. parse JSON args (errors-as-result, not thrown)                     │
│        2. marshal flattened args → channel values                            │
│           (ProcessEntryHandler.getProcessArguments — same path as            │
│            `nextflow module run`)                                            │
│        3. BIND each value onto the pre-wired input queue(s)                  │
│        4. BLOCK on the pre-captured output read channels (.val) — the task   │
│           runs on the real executor/container/work dir, caching applies      │
│        5. serialize outputs → JSON (path ⇒ absolute-path handle, OR small    │
│           structured text inlined; skip topic/eval bookkeeping outputs)      │
│        6. return JSON tool result to the loop                                │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Why tools are pre-wired (the central decision)

A tool call must run *synchronously* from inside the agent's map closure — the LLM is blocked waiting for the result. Invoking a `ProcessDef` synchronously *after* the dataflow network has ignited deadlocks GPars. The bridge therefore **pre-wires** each tool in the workflow body (before ignition): it creates one persistent queue per input channel, runs the cloned process over those queues once, and **captures the output read channels at wiring time**. At dispatch it merely `bind()`s args onto the queues and blocks on `.val`. Capturing the read channels at wiring (not lazily at dispatch) is required because process outputs are broadcast-style — a subscriber created after a value binds never sees it, which would deadlock the 2nd+ output of a multi-output tool.

### Why dispatch is serialized (not concurrent)

The original draft proposed running concurrent tool calls in parallel. The implementation does the opposite: `ModuleToolBridge.call` is **`synchronized`**. The pre-wired queues are shared per tool, so binding args and pulling `.val` must stay paired; serializing dispatch keeps each call's input correlated with its output. Multiple tool requests in one assistant turn are executed one at a time.

### Why in-JVM rather than subprocess

- No external binary or container daemon dependency — the agent feature is available wherever Nextflow is.
- Per-invocation overhead is the chat-model HTTP call, not a subprocess startup.
- Tool callbacks dispatch directly into the dataflow runtime; no MCP serialization round-trip.
- Errors propagate as Java exceptions; observability is direct; tests are plain Spock — no subprocess mocking.

## 4. Tool bridge: module → tool descriptor & dispatch

A tool descriptor (name, description, input JSON-schema) is built **once per tool at wiring time**. The descriptor source is chosen by availability:

| Source | When | Class |
|---|---|---|
| Public registry `ModuleMetadata` (`GET /api/v1/modules/{name}`) | registry tool resolves and metadata is reachable — **the canonical, richer source** (per-field descriptions, patterns, enums, tool homepages, nf-core `meta.id` convention) | `ModuleMetadataToolSchema` |
| Sibling `meta.yml` `ModuleSpec` | local `.nf` tool, or registry tool offline / metadata fetch failed (graceful fallback) | `ModuleSpecToolSchema` |
| Declared typed process I/O | in-scope process tool (scalar mode, no `meta.yml`) | `ProcessToolSchema` |

Schemas are **flattened**: a tuple input channel (`tuple(meta, reads)`) contributes one top-level property per component, so the LLM passes `{"meta": {...}, "reads": "/abs/path"}` rather than a nested array. For an nf-core `map` input named `meta`, a nested `{id: string}` sub-schema is emitted (the `meta.id` convention) so the model can build `meta.id` from the schema alone. The portable schema is a plain `Map` in core (`ToolDescriptor`); the plugin maps it onto langchain4j (`ModuleToolAdapter` → `JsonSchemaMapper`), keeping core free of any LLM-client dependency.

When both the registry metadata (descriptor) and the installed `meta.yml` spec (marshalling) are present, a **drift guard** warns if their flattened input names differ.

### Per tool-call flow

1. **Arg parsing** — the LLM-supplied JSON is parsed to a `Map`. Empty payload ⇒ no args. Unparseable / non-object payloads, an unknown tool name, or marshalling failures are returned as a `{"error": "..."}` tool result (not thrown), so the LLM can retry (counts toward `maxIterations`).
2. **Marshalling** — `ProcessEntryHandler.getProcessArguments` (the *same* path as `nextflow module run`) treats the flattened args as the params map, looks up each declared input by name, type-coerces per the spec (`file`→`Nextflow.file`, `map`→`Map`, `integer`→`Integer`, …) and assembles tuple inputs in order (`[[id:'s1'], file(reads)]`). A null **or blank-string** path arg is treated as *not provided* (`[]`) — the nf-core optional-path convention — so `file("")` never throws.
3. **Module invocation** — the value is bound onto the pre-wired queue; the cloned `ProcessDef` runs on whatever executor/container the user configured. Nothing in the agent path overrides that.
4. **Output serialization** — read each output channel in declared order from the channels captured at wiring; map back per output param. `eval`/`topic`-routed bookkeeping outputs (nf-core `versions`) are **skipped authoritatively from the ProcessDef** (`isTopicOutput`, guarded to `BaseOutParam`) because their channel never binds a per-call value and `.val` would block forever.
5. **File outputs** — `ToolOutputReader.readOrHandle`: bulk/binary/oversized → **absolute path string** (opaque handle); small text/structured files (`.json/.tsv/.csv/.txt/.yaml/…`, under `agent.maxToolOutputInlineSize`, default 32 KB, non-binary) → **contents inlined** so the LLM can reason over them. This handle-vs-inline split is what enables both downstream chaining (handles) and data-driven control flow (inlined stats gates).
6. **Return** — a single JSON string fed back to the loop as a `ToolExecutionResultMessage`.

### Path handling (v1 contract)

Bulk/binary paths are **opaque handles** — the LLM sees absolute path strings, chains them between tools, and never reads/writes contents. Small structured text outputs are the deliberate exception (inlined). The output description prepended to each tool makes the opaque-path contract explicit.

## 5. Configuration & credentials

`nextflow.config` gains an `agent` scope (`AgentConfig`). All values are **defaults** — an agent's own directives take precedence.

```groovy
agent {
    defaultModel            = 'openai/gpt-5-mini'   // used when an agent declares no `model`
    maxIterationsDefault    = 20                    // used when an agent declares no `maxIterations`
    requestTimeout          = '120s'                // per-request LLM chat timeout
    maxToolOutputInlineSize = '32 KB'               // cap for inlining structured tool outputs
}
```

The API key is read from `OPENAI_API_KEY` (the only provider in v1). No new secret-management surface. The `nf-agent` plugin (see §7) brings in langchain4j; the user opts in by declaring the plugin.

## 6. Error handling & limits

| Failure | Behavior |
|---|---|
| `nf-agent` plugin not enabled (no `AgentRunner`) | `AgentRunnerProvider.get()` aborts with a hint to enable the plugin |
| Provider in `model` not `openai` | `ChatModelFactory.createModel` throws at **run time** (first input record) — *not* at agent construction |
| `OPENAI_API_KEY` missing | throws at **run time** in `ChatModelFactory` with a hint naming the env var |
| Dispatch-level tool error (unknown tool / unparseable args / bad marshalling) | returned to the LLM as `{"error": "..."}`; the model may correct and retry (counts toward `maxIterations`) |
| Underlying tool **task** fails | **not** recovered as a tool result — surfaces through the process's own dataflow operator and aborts the session (standard error model). Documented limitation, future hardening. |
| LLM run exceeds `maxIterations` | the runner throws `IllegalStateException`; the agent operator's `onException` poisons the tool queues and aborts the session |
| LLM provider error (rate limit, network, auth) | surfaces as the operator failing → session abort |

## 7. Module layout

The DSL keyword and AST stay in **core** (the parser must recognize `agent`); the langchain4j client lives in the **`nf-agent` plugin**. **All** dataflow and tool-bridge logic is in core — the plugin is intentionally thin (model construction + the chat/tool loop).

### Core (always loaded)

- `modules/nf-lang/src/main/antlr/ScriptLexer.g4` / `ScriptParser.g4` — `AGENT` token + `agentDef`/`agentBody`/`agentDirectives`/`agentInputs`/`agentOutputs`/`agentPrompt` rules (inputs/outputs reuse `processInput`/`processOutput`).
- `modules/nf-lang/.../ast/AgentNode.java` + `control/AgentToGroovyVisitorV2.java` — AST node and lowering to a runtime `agent('name', { … })` call.
- `nextflow.script.AgentBuilder` — runtime delegate that captures directives/inputs/outputs/prompt and builds the `AgentDef` (**not** an `AgentFactory`).
- `nextflow.script.AgentDef` — the runtime model + operator (`run`, `createToolBridge`, tool resolution/compilation, `applyAgentOperator`).
- `nextflow.script.ProcessEntryHandler` — LLM-args → channel-values marshalling (shared with `nextflow module run`).
- `nextflow.agent.ModuleToolBridge` — the dispatcher (`ToolDispatcher`): pre-wire, marshal, dispatch, serialize, poison-on-close.
- `nextflow.agent.{AgentConfig, AgentRunner, AgentRunnerProvider, AgentRunnerRequest, ToolDispatcher, ToolDescriptor, ToolOutputReader}` — config scope, SPI, plugin lookup, request DTO, dispatch interface, portable descriptor, output reader.
- `nextflow.agent.{ModuleMetadataToolSchema, ModuleSpecToolSchema, ProcessToolSchema, RecordSchema}` — descriptor/schema derivation from the three sources + record reflection.

### `plugins/nf-agent` (opt-in)

- `AgentPlugin` — pf4j plugin entry.
- `LangChainAgentRunner` — implements `AgentRunner`: single-shot (optionally structured) chat, or the manual tool-call loop.
- `ChatModelFactory` — `provider/model` → langchain4j `ChatModel`. **v1: OpenAI only.**
- `ModuleToolAdapter` — portable `ToolDescriptor` → langchain4j `ToolSpecification`.
- `JsonSchemaMapper` — portable schema `Map` → langchain4j `JsonSchema` (structured output + tool parameters).
- `build.gradle` dep — **`dev.langchain4j:langchain4j-open-ai:1.15.1` only** (no anthropic/gemini modules in v1).

### Existing infrastructure reused

- `nextflow.module.{ModuleSpec, ModuleSpecFactory, ModuleReference, ModuleResolver, RegistryClientFactory}` — spec loading, registry resolution/auto-install, metadata fetch.
- Standard include/registry resolution + `ScriptLoaderV2` (per-module `ScriptBinding`) — for compiling tool modules in isolation.
- `ProcessDef` / `ChannelOut` / dataflow runtime — for actually executing tool invocations.
- Plugin loader — for discovering and activating `nf-agent`.

## 8. Open questions (resolved)

| # | Question | Resolution |
|---|---|---|
| 1 | Path representation | String handle (absolute path), **except** small structured text outputs which are inlined |
| 2 | Tool source | In-scope process, local `.nf`, or registry ref; descriptor from registry `ModuleMetadata` (fallback `meta.yml`) |
| 3 | LLM↔dataflow mapping | **Pre-wired** tool queues; dispatch binds args + blocks on captured output channels — see §3 |
| 4 | Agent-as-tool | Not in v1 |
| 5 | Runtime engine | langchain4j (in-JVM), **OpenAI only** in v1 |
| 6 | Provider/model identifier | `provider/model`, parsed by `ChatModelFactory` |
| 7 | Concurrent tool calls | **Serialized** (`synchronized` dispatch over shared pre-wired queues) — reversed from the draft |

## 9. Future extensions

- **Additional providers** — add langchain4j modules (anthropic, gemini, bedrock, vertex) and extend `ChatModelFactory` beyond OpenAI.
- **Structured outputs *with* tools** — today they are mutually exclusive; a future phase could request a JSON-schema response after the tool loop converges.
- **Multiple inputs/outputs & `path` declarations at the agent boundary** — currently exactly one input and one output.
- **External MCP toolsets** — treat any MCP server as another `tools` source via `langchain4j-mcp`.
- **Agent-as-tool** — allow one `agent` in another agent's `tools` list.
- **Streaming partial outputs** — surface intermediate tokens / tool-call traces to channels.
- **Tool-task failure recovery** — turn an underlying tool task failure into a recoverable tool result instead of a session abort.
- **Genuinely-absent optional outputs** — a tool that produces no file for an `optional: true` output currently blocks the per-call `.val` (see the implementation guide §10). Needs a core dataflow "empty for this task" sentinel or a bridge-local guard.

## 10. Divergences from the original design

The first draft (2026-05-04) and the shipped v1 differ in these material ways:

1. **Concurrency reversed** — draft said concurrent tool dispatch; v1 serializes (`synchronized`) because tool queues are shared and pre-wired.
2. **Pre-wiring** — draft created channels per call (`Channel.value`/`Channel.fromPath`); v1 pre-wires persistent queues before ignition and captures output read channels at wiring (a GPars-deadlock requirement the draft didn't anticipate).
3. **Single provider** — draft implied multi-provider via langchain4j; v1 ships OpenAI only.
4. **Tool-descriptor source** — draft sourced everything from `meta.yml`; v1 prefers the public registry `ModuleMetadata` and falls back to `meta.yml`.
5. **Output inlining** — draft treated all paths as opaque handles; v1 inlines small structured text outputs so the LLM can reason over them (new `agent.maxToolOutputInlineSize`).
6. **Structured outputs** — listed as a non-goal in the draft; implemented in v1 for the no-tools path (tools XOR structured).
7. **Typed DSL I/O** — draft used `val name` process-style inputs; v1 uses typed `name: Type` declarations (one input, one output).
8. **Error timing & layout** — provider/key checks fire at run time (not construction); the dispatcher (`ModuleToolBridge`) and `ToolDispatcher` live in **core**, not the plugin; the builder is `AgentBuilder`, not `AgentFactory`.
