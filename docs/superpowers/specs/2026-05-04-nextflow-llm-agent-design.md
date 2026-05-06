# Nextflow LLM Agent — Design

**Date:** 2026-05-04
**Status:** Draft (pending implementation plan)
**Author:** Paolo Di Tommaso

## 1. Goal & non-goals

### Goal

Add a new top-level `agent` construct to the Nextflow DSL — a process-shaped primitive that wraps an LLM-driven tool-calling loop. Each invocation of an agent receives one input record, runs an LLM with access to a declared set of Nextflow modules as tools, and emits one output record. Agents compose with processes and other agents through the existing channel/workflow model.

The defining property of the design: **a tool call from the LLM is not a function call — it is a real Nextflow module invocation.** The args are wrapped as channel values, the module runs through the standard dataflow runtime (executor, container, retries, work dir, caching), and the resulting output channel record is serialized back to the LLM. The harness is the adapter between the LLM tool-call protocol and Nextflow's dataflow model.

### Non-goals (v1)

- Long-lived conversational agents (one input → one output is the contract)
- Channel-aware orchestrator agents (agents do not subscribe to channels mid-loop)
- Built-in multi-provider LLM abstraction (delegated to langchain4j)
- Agent-to-agent invocation as tools (composition is via channels only)
- Cost tracking, prompt analytics, fine-grained token budgets
- RAG, persistent memory, conversation state across invocations
- Typed structured outputs from the LLM (deferred until record-types lands)
- Streaming partial outputs (langchain4j supports it; surfaced as a follow-up)

## 2. DSL surface

A new top-level `agent` definition, parsed by `nf-lang` as a sibling of `processDef` and `workflowDef`. Tools are referenced as **modules** (registry coordinate or local include), so they bring their own description and I/O metadata via `meta.yml`.

```groovy
include { fastqc }     from 'nf-core/fastqc'
include { multiqc }    from 'nf-core/multiqc'
include { fetch_sra }  from './local/fetch_sra'

agent eval_agent {
    model 'openai/gpt-5-mini'
    instruction 'You break biological questions into structured plans.'
    tools fastqc, multiqc, fetch_sra
    maxIterations 20

    input:
        val question

    output:
        val plan

    prompt:
    """
    Question: ${question}
    Plan the steps to answer this using the available tools.
    """
}

workflow {
    Channel.of(params.question)
        | eval_agent
        | view
}
```

### Directives

| Directive | Required | Meaning |
|---|---|---|
| `model` | yes | Model identifier in `provider/model` form (e.g. `openai/gpt-5-mini`, `anthropic/claude-opus-4-7`). Provider prefix selects the langchain4j chat-model factory. |
| `instruction` | yes | System prompt (agent role/persona) |
| `tools` | yes | Comma-separated module references resolved through standard include/registry machinery |
| `maxIterations` | no | LLM-loop iteration cap (default from `agent.maxIterationsDefault` config) |
| `input:` / `output:` | yes | Standard process-style channel I/O |
| `prompt:` | yes | Templated user prompt with `${var}` interpolation over input bindings |

### I/O contract (v1)

- Input: any number of `val` declarations bound into the `prompt:` template via `${name}` interpolation. `path` inputs are not supported in v1 (out of scope).
- Output: a single `val` declaration that receives the LLM's final assistant message as a string. Multiple outputs and typed/structured outputs are deferred — see §9.

Standard process directives (`errorStrategy`, `maxRetries`, `time`, `cache`) apply to agents as they do to processes; agent failures use existing retry/error machinery without a parallel control path.

## 3. Architecture

The runtime is **in-JVM**. Each agent invocation builds a langchain4j chat model, registers the agent's tools as Java callbacks, and drives the tool-call loop directly. No subprocesses, no external binaries, no IPC.

```
┌─────────────────────────── Nextflow JVM ─────────────────────────────┐
│                                                                      │
│  workflow ──► AgentTaskRun  (one per input record)                   │
│                    │                                                 │
│                    ├── resolves tool list                            │
│                    │     • each ref → ModuleSpec (from meta.yml)     │
│                    │     • spec → langchain4j ToolSpecification      │
│                    │                                                 │
│                    ├── builds ChatModel (per `model` directive)      │
│                    │     • provider prefix → factory                 │
│                    │     • API key from env (provider-standard var)  │
│                    │                                                 │
│                    ├── renders `prompt:` template                    │
│                    │                                                 │
│                    └── runs agent loop (≤ maxIterations):            │
│                          1. send messages to ChatModel               │
│                          2. if assistant returns tool_calls,         │
│                             dispatch each → ToolDispatcher           │
│                             (concurrent if multiple)                 │
│                          3. append tool results, repeat              │
│                          4. on plain assistant message → done        │
│                                                                      │
│  on tool call from LLM:                                              │
│      ToolDispatcher.invoke(toolName, jsonArgs)                       │
│         1. validate jsonArgs against ModuleSpec.input schema         │
│         2. wrap each arg as a one-shot channel value                 │
│              (val→Channel.value, path→Channel.fromPath, etc.)        │
│         3. invoke ModuleDef.apply(channels) — full dataflow path:    │
│              executor, container, retries, work dir, caching all     │
│              behave exactly as a normal pipeline invocation          │
│         4. await ChannelOut, serialize each output per .output       │
│              schema (path → absolute string handle, val → JSON)      │
│         5. return JSON tool result to the loop                       │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘
```

### Per-invocation lifecycle

1. `AgentTaskRun` materializes the tool list, resolving each reference to a `ModuleSpec` and a callable `ModuleDef`.
2. Builds a langchain4j `ChatModel` from the `model` directive (provider prefix selects the factory; api key from env).
3. Builds a `ToolSpecification` for each tool from its `ModuleSpec`.
4. Renders the `prompt:` template against input bindings.
5. Runs the agent loop: send messages → handle `tool_calls` via `ToolDispatcher` → repeat until the LLM emits a plain assistant message or `maxIterations` is hit.
6. Final assistant message text is bound to the agent's `output:` channel.

### Why in-JVM rather than subprocess

- No external binary or container daemon dependency — the agent feature is available wherever Nextflow is.
- Per-invocation overhead is the chat-model HTTP call, not a subprocess startup; cleanly supports workflows that fan out hundreds of agent calls.
- Tool callbacks dispatch directly into the dataflow runtime; no MCP serialization round-trip.
- Errors propagate as Java exceptions, observability is direct (token counts, tool traces), and tests are plain JUnit/Spock — no subprocess mocking.

## 4. Tool bridge: ModuleSpec → langchain4j ToolSpecification

The bridge consumes the already-validated `ModuleSpec` (loaded by `ModuleSchemaValidator`, introduced in 26.04 via commit `571274552`) and produces a langchain4j `ToolSpecification`:

| `ModuleSpec` field | `ToolSpecification` field |
|---|---|
| `name` | tool name (registry-qualified, e.g. `nf-core/fastqc`) |
| `description` | tool description |
| `tools[].description` (joined) | appended to tool description for richer LLM context |
| `input[]` (paramSpec list) | `JsonObjectSchema` properties — each param's `name`, `type`, `description` |
| `output[]` (paramSpec list) | included in tool description as a hint (langchain4j tools don't carry an output schema; the LLM learns the shape from the returned JSON and the description) |

Each tool is registered with the chat model's tool list. The `ToolDispatcher` is a single Java callback (one method per tool) that receives `(toolName, ToolExecutionRequest)` and returns a JSON string.

### Per tool-call flow

The steps below describe how the harness handles a **single LLM tool-call** — i.e. one invocation of one of the modules listed in the agent's `tools`. They are distinct from the agent's own `input:` block (covered in §2), which is rendered into the prompt template once per agent invocation.

1. **Arg validation** — JSON Schema check of the LLM-provided args against `ModuleSpec.input`. On failure the dispatcher returns a tool-result string describing the error; the LLM may retry (counts toward `maxIterations`).
2. **Channel materialization** — for each tool-call arg, build a one-shot input channel matching the module's declared input shape:
   - `val` arg → `Channel.value(arg)`
   - `path` arg → `Channel.fromPath(arg)` (the LLM passes a path string previously emitted by another tool — a handle, not contents)
   - `tuple` arg → composed value channel matching the tuple shape
3. **Module invocation** — `ModuleDef.apply(channels)` returns a `ChannelOut`. The module runs on whatever executor and container the user has configured for it; nothing in the agent path overrides that.
4. **Output serialization** — collect the first record from each output channel; map back per `output[]` paramSpec:
   - `path` → absolute path string (handle)
   - `val` → JSON-serialized value
   - `tuple` → object with named fields
5. **Return** — single JSON string returned to langchain4j as the tool result; LLM continues its loop.

### Concurrent tool calls

If the LLM emits multiple `tool_calls` in one assistant turn, the dispatcher runs them concurrently. Each invocation is independent — different work dirs, no shared state — matching standard Nextflow process semantics.

### Path handling (v1 contract)

Paths are **opaque handles**. The LLM sees absolute path strings, passes them between tools, and never reads or writes the contents directly. The system prompt the harness prepends to `instruction` should make this contract explicit so the LLM does not attempt to inline file contents.

## 5. Configuration & credentials

`nextflow.config` gains an `agent` scope:

```groovy
agent {
    defaultModel         = 'openai/gpt-5-mini'
    maxIterationsDefault = 20
    requestTimeout       = '120s'
}
```

API keys are read from the launching environment using each provider's standard variable (`OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `GOOGLE_API_KEY`, etc.). No new secret-management surface in v1.

The `nf-agent` plugin (see §7) brings in the langchain4j core and per-provider modules. The user enables providers by declaring the plugin; no external binary or container daemon is required.

## 6. Error handling & limits

| Failure | Behavior |
|---|---|
| Provider for `model` not enabled (no factory loaded) | Fail fast at agent construction with hint to enable the matching langchain4j module |
| API key missing for selected provider | Fail fast at agent construction with hint naming the expected env var |
| Tool (module) invocation fails | Return error string as tool result; LLM may recover or give up |
| LLM run exceeds `maxIterations` | AgentTaskRun fails with iteration-cap error (standard Nextflow `errorStrategy` applies) |
| LLM provider error (rate limit, network, auth) | Surface as task failure; user opts in via `errorStrategy 'retry'` |
| Schema validation fails on tool args | Return error string as tool result; counts toward iteration budget |

Agents inherit the standard process directives `errorStrategy`, `maxRetries`, `time`, `cache`. No parallel retry logic.

## 7. Module layout

The DSL keyword and AST stay in **core** (the parser must recognize `agent`); the runtime, langchain4j dependencies, and provider modules live in a new **`nf-agent` plugin**. This keeps core's classpath clean and lets users opt in.

### Core (always loaded)

- `modules/nf-lang/src/main/antlr/ScriptParser.g4` — add `agentDef` rule, parallel to `processDef`. AST node + visitor support in the v2 parser package.
- `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` — runtime model for an agent definition.
- `modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy` — construction from script body.
- `modules/nextflow/src/main/groovy/nextflow/agent/AgentRunner.groovy` — SPI interface implemented by the plugin.
- `modules/nextflow/src/main/groovy/nextflow/agent/AgentConfig.groovy` — `agent` config scope binding.

### `plugins/nf-agent` (opt-in)

- `LangChainAgentRunner.groovy` — implements the `AgentRunner` SPI: builds `ChatModel`, registers tools, runs the loop.
- `ToolDispatcher.groovy` — single Java callback that receives `(toolName, jsonArgs)` and returns a JSON string.
- `ModuleToolAdapter.groovy` — one `ModuleSpec` → one `ToolSpecification` + dispatch handler.
- `ChatModelFactory.groovy` — `provider/model` → `ChatModel` (delegates to langchain4j's per-provider builders).
- `build.gradle` deps — `langchain4j-core`, plus `langchain4j-open-ai`, `langchain4j-anthropic`, `langchain4j-google-ai-gemini`, etc. (kept narrow for v1; expand on demand).

### Existing infrastructure reused

- `nextflow.module.ModuleSpec` and `nextflow.module.ModuleSchemaValidator` (26.04+) — the spec source.
- Standard include/registry resolution machinery — for resolving `tools` references.
- `ProcessDef` / `ChannelOut` / dataflow runtime — for actually executing tool invocations.
- Plugin loader — for discovering and activating `nf-agent`.

## 8. Open questions (resolved)

| # | Question | Resolution |
|---|---|---|
| 1 | Path representation | String handle (absolute path) — opaque to the LLM |
| 2 | Tool source | Modules from registry or local include — leverages `meta.yml` |
| 3 | LLM↔dataflow mapping | Tool callback wraps args as channel values, awaits outputs, serializes back — see §3 |
| 4 | Agent-as-tool | Not in v1 |
| 5 | Runtime engine | langchain4j (in-JVM) — see §3 rationale |
| 6 | Provider/model identifier | `provider/model` convention, parsed by `ChatModelFactory` |
| 7 | Concurrent tool calls | Run concurrently; one work dir per call (matches process semantics) |

## 9. Future extensions

- **Typed/structured LLM outputs** — once record-types lands, agents declare typed `output:` and the runtime requests a JSON-schema response from the LLM (langchain4j supports structured outputs natively), validated and bound to typed channel records.
- **Multiple outputs and `path` inputs/outputs at the agent boundary** — currently restricted to single `val` output for simplicity.
- **External MCP toolsets** — treat any MCP server (Seqera Cloud MCP, third-party) as another `tools` source via `langchain4j-mcp`. Aligns with the Foundry pattern, inverted (Nextflow as host).
- **Agent-as-tool** — allow one `agent` to be referenced in another agent's `tools` list, enabling hierarchical agent composition without channel plumbing.
- **Streaming partial outputs** — surface intermediate LLM tokens or tool-call traces to channels for live observability (langchain4j has a streaming chat API).
- **Alternative runners** — the `AgentRunner` SPI in core leaves room for a future `nf-agent-bedrock`, `nf-agent-vertex`, or `nf-agent-docker` (subprocess to docker-agent) without changes to the DSL or core.
