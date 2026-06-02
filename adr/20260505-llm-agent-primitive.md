# LLM Agent Primitive

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2026-05-05
- Tags: dsl, agent, llm, language

## Summary

Introduce a new top-level `agent` construct to the Nextflow DSL — a process-shaped primitive that wraps an LLM-driven tool-calling loop. Agents accept input records, run an LLM with access to a declared set of Nextflow modules as tools, and emit output records. Tool calls from the LLM are dispatched as real Nextflow module invocations, going through the standard dataflow runtime.

## Problem Statement

Bioinformatics pipelines are increasingly being orchestrated by LLM agents (Microsoft Foundry, OpenAI Assistants, Claude tool use, etc.) that interpret biological intent, parameterise pipelines, and post-process results. Today this orchestration lives **outside** Nextflow — agents in Foundry call into Seqera Cloud's MCP server to launch nf-core pipelines.

This design proposes the inverse: agents become a first-class Nextflow construct, defined inside `.nf` scripts alongside processes and workflows. The pipeline reasons about itself, with the LLM choosing which processes to invoke based on intermediate data — closing the loop between hypothesis and computation without an external orchestrator.

## Goals or Decision Drivers

- **First-class language construct**: agents compose with processes and other agents through the existing channel/workflow model.
- **Process-shaped primitive**: one input record → one output record. Multi-turn LLM↔tool reasoning happens inside the agent for each input; composition of multiple agents happens at the workflow level.
- **Tool calls are real dataflow nodes**: when the LLM calls a tool, it is **not** a function call — args are wrapped as channel values, the module runs through the standard executor/container/retry/cache machinery, and outputs flow back to the LLM as JSON.
- **Tools are registry modules**: each tool brings its own description and I/O metadata via `meta.yml`, so the LLM gets a JSON-schema tool descriptor for free.
- **In-JVM runtime, no external binary**: delegate multi-provider LLM client work to **langchain4j**; no Docker Desktop / Go subprocess prerequisite.
- **Backward compatibility**: zero impact on existing `process`/`workflow`/`function` declarations.

## Non-goals (v1)

- Long-lived conversational agents (one input → one output is the contract).
- Channel-aware orchestrator agents (agents do not subscribe to channels mid-loop).
- Agent-to-agent invocation as tools (composition is via channels only).
- Cost tracking, prompt analytics, fine-grained token budgets.
- RAG, persistent memory, conversation state across invocations.
- Typed structured outputs from the LLM (deferred until record-types lands).
- Streaming partial outputs (langchain4j supports it; surfaced as a follow-up).

## Solution or decision outcome

A new `agent` keyword is added to the Nextflow DSL with the same composition model as `process`. Tools are referenced as registry modules (the existing `include { fastqc } from 'nf-core/fastqc'` syntax). The runtime engine is **langchain4j** in a new `nf-agent` plugin; the DSL keyword and AST live in core (parser must recognise `agent`).

```groovy
include { fastqc } from 'nf-core/fastqc'

agent eval_agent {
    model 'openai/gpt-5-mini'
    instruction 'Plan FASTQ analysis steps.'
    tools fastqc

    input:
        question: String

    output:
        plan: String

    prompt:
    """
    Question: ${question}
    """
}

workflow {
    Channel.of(params.question)
        | eval_agent
        | view
}
```

## Rationale & discussion

The defining design choice is **A. agent-as-process** over the alternatives (B. channel-aware orchestrator, C. long-lived service). Each agent invocation is one black-box record-in / record-out node; multi-turn LLM reasoning happens internally, and multi-agent pipelines compose like multi-process pipelines. This reuses the existing channel composition story and only adds one new node type.

The runtime engine choice — **langchain4j** over docker-agent — eliminates the Docker Desktop prerequisite, removes subprocess startup overhead, simplifies the tool-call path (no MCP-over-stdio bridge), and supports the same multi-provider abstraction natively.

Tools-as-modules over a custom tool format reuses the `meta.yml` schema (validated by `ModuleSchemaValidator`, introduced in 26.04), so tool descriptors come for free with no parallel metadata layer.

Full design spec: [`docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md`](../docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md).

Implementation plan: [`docs/superpowers/plans/2026-05-04-agent-syntax-and-model.md`](../docs/superpowers/plans/2026-05-04-agent-syntax-and-model.md).

## Implementation status (as of 2026-05-05)

### Delivered — minimum working POC on branch `agent-syntax-and-model`

The `agent` keyword is recognised end-to-end. A script containing `agent foo { prompt: "..." }` parses, lowers to executable Groovy, registers an `AgentDef` on the script's `ScriptMeta`, and would invoke the agent body on workflow execution (where it currently throws `UnsupportedOperationException("agent execution not yet implemented")`).

| Layer | Change |
|---|---|
| Lexer (`ScriptLexer.g4`) | New `AGENT` and `PROMPT` tokens |
| Parser (`ScriptParser.g4`) | New `agentDef`/`agentBody`/`agentDirectives`/`agentInputs`/`agentOutputs`/`agentPrompt` rules; `AGENT` and `PROMPT` added to `identifier` and `keywords` rules |
| AST | New `AgentNode` class; `ScriptNode.agents` collection; `ScriptVisitor.visitAgent` + `ScriptVisitorSupport` default; agents in `getDeclarations()` |
| AST builder (`ScriptAstBuilder.java`) | `agentDef(...)` + helpers; agents added to `SCRIPT_DEF_NAMES`; `AgentDefAltContext` dispatch |
| Lowering (`ScriptToGroovyVisitor.java`) | `visitAgent` emits `agent('name', { /* empty closure */ })` |
| Resolution (`ScriptResolveVisitor`, `VariableScopeVisitor`, `VariableScopeChecker`, `ResolveIncludeVisitor`) | Agents iterate alongside processes; `isDataflowMethod` recognises `AgentNode` so pipe syntax works; agents are includable across scripts |
| Runtime (`AgentDef`, `AgentFactory`) | Mirror `ProcessDef`/`ProcessFactory`. `AgentDef extends BindableDef implements ChainableDef`; `run()` throws `UnsupportedOperationException` |
| Runtime DSL (`BaseScript.agent(name, body)`) | Constructs `AgentDef` via `AgentFactory`, registers via `meta.addDefinition(...)` |
| Tests | 8 new Spock tests across `AgentParserTest`, `AgentDefTest`, `AgentFactoryTest`, `AgentScriptLoadingTest` — all passing |

19 implementation commits; no regressions in `nf-lang`/`nextflow` test suites (2 unrelated `RegistryClientFactoryTest` failures pre-existed on `master`).

### Known limitations / what's missing

1. ~~**Agent body resolution**~~ — **DELIVERED (Plan B, 2026-06-01)**. The full directive-rich form (`model`, `instruction`, `tools()`, `maxIterations`, typed `input:`/`output:`, and a `${var}`-interpolated `prompt:`) now resolves under typed mode and loads via `ScriptLoaderV2`. Added `AgentDsl` (parallel to `ProcessDsl`, declaring the four directive method shapes) and a `VariableScopeVisitor.visitAgent` override that pushes the DSL scope, declares input parameters as locals, and visits directives/prompt/outputs. See [`docs/superpowers/plans/2026-06-01-agent-body-resolution.md`](../docs/superpowers/plans/2026-06-01-agent-body-resolution.md). Note: `tools` with module references (`tools fastqc`) is still deferred to the tool-bridge work — only empty `tools()` is exercised — and agent outputs are not wrapped in an output DSL scope (intentional, since agent outputs are plain `name: Type` declarations).
2. **Better missing-prompt error** — the ANTLR `agentBody` rule enforces `prompt:` non-optionally, so a missing prompt produces a generic "Invalid agent definition" error rather than the intended specific message. The AST-builder guard `if( ctx.body.agentPrompt() == null ) return invalidAgent("Missing prompt: section", ...)` is dead code. Make `agentPrompt` grammar-optional and let the AST-builder guard fire.
3. ~~**Lowered closure body is empty**~~ — **DELIVERED (Plan C, 2026-06-01)**. `AgentToGroovyVisitorV2` now lowers `AgentNode` into a populated `agent('name', { <directives>; _input_(...); _output_(...); new PromptDef(...) })` closure; `BaseScript.agent` runs it against an `AgentBuilder` delegate and builds an `AgentDef` exposing `model`/`instruction`/`tools`/`maxIterations`/`inputs`/`outputs`/`prompt`. `AgentFactory` was removed (superseded by `AgentBuilder`). See [`docs/superpowers/plans/2026-06-01-agent-closure-population.md`](../docs/superpowers/plans/2026-06-01-agent-closure-population.md). Note: unrecognized `output:` targets are silently dropped (no `$out` fallback like processes) — fine for the single-`val` v1 contract.
4. ~~**Agent execution is stubbed**~~ — **DELIVERED (Plan D, 2026-06-02)** for the no-tools path. `AgentDef.run` now executes as a dataflow `MapOp` operator: per input record it renders the prompt and delegates to an `AgentRunner` SPI (`nextflow.agent.AgentRunner`/`AgentRunnerRequest`/`AgentRunnerProvider`) discovered via `Plugins.getPriorityExtensions`. The new `plugins/nf-agent/` plugin provides `LangChainAgentRunner` (langchain4j 1.15.1, single-shot chat) + `ChatModelFactory` (`provider/model`, OpenAI). See [`docs/superpowers/plans/2026-06-02-agent-execution.md`](../docs/superpowers/plans/2026-06-02-agent-execution.md). **Scoped v1:** exactly one `val` input, OpenAI provider only, no tools (the request carries `tools` but the runner ignores them) — multiple/zero inputs and non-OpenAI providers fail with a clear error. Tested without API keys (mock runner in core, mock `ChatModel` in the plugin).
5. ~~**Tool bridge is unimplemented**~~ — **DELIVERED (Plan G Phase 2, 2026-06-02)** for in-scope processes with scalar I/O. Core `ModuleToolBridge`/`ToolDescriptor`/`ToolDispatcher` pre-wire each declared process, marshal the LLM's JSON args into channel values, run the process as a real dataflow node, and serialize the output back to JSON; the plugin's adapter builds langchain4j `ToolSpecification`s and runs the tool-call loop. Tools require a non-record output (tools XOR structured output). Registry-module (`nf-core/...`) tools and `meta.yml`-derived schemas are still future (staged-scope item 3). See the real-run evidence under "Staged scope" above.
6. **`agent` config scope** — `nextflow.config` `agent { defaultModel = ..., maxIterationsDefault = ... }` block is not parsed yet.

### Next steps

In rough order of dependency:

1. ~~**Plan B: Agent body resolution**~~ — **DONE (2026-06-01)**, see [`docs/superpowers/plans/2026-06-01-agent-body-resolution.md`](../docs/superpowers/plans/2026-06-01-agent-body-resolution.md). `AgentDsl` + `VariableScopeVisitor.visitAgent` landed; the directive-rich agent now resolves and loads. The next blocker is Plan C (closure population), since the lowered closure is still empty.
2. ~~**Plan C: Closure population**~~ — **DONE (2026-06-01)**, see [`docs/superpowers/plans/2026-06-01-agent-closure-population.md`](../docs/superpowers/plans/2026-06-01-agent-closure-population.md). `AgentToGroovyVisitorV2` + `AgentBuilder` + populated `AgentDef` landed. The runtime now has the full declarative content; the next blocker is Plan D (execution).
3. ~~**Plan D: nf-agent plugin**~~ — **DONE (2026-06-02)** for the no-tools, OpenAI, single-input path, see [`docs/superpowers/plans/2026-06-02-agent-execution.md`](../docs/superpowers/plans/2026-06-02-agent-execution.md). The `AgentRunner` SPI, langchain4j `ChatModelFactory`, and `LangChainAgentRunner` landed. Still future: the `ToolDispatcher` (tool→module bridge), additional provider modules (`langchain4j-anthropic`, `-google-ai-gemini`), and multi-input support.
4. ~~**Plan E: End-to-end demo**~~ — **DONE (2026-06-02)**, see [`docs/superpowers/plans/2026-06-02-agent-demo-and-docs.md`](../docs/superpowers/plans/2026-06-02-agent-demo-and-docs.md). Added the `docs/agent.md` user guide and a gated (`@Requires(OPENAI_API_KEY)`) end-to-end test. **Verified with a real run:** `./launch.sh run` of an agent pipeline with `plugins { id 'nf-agent@0.1.0' }` (dev mode) produced a live OpenAI `gpt-5-mini` answer (`ANSWER=Paris`), with PF4J discovering `LangChainAgentRunner` via the real extension point (not the test seam). Tool dispatch through a real module is still future (depends on the `ToolDispatcher`).
5. **Polish backlog** (any time): make `agentPrompt` grammar-optional for better errors; add `agent` config scope; add an `@author` tag style sweep across new files; revisit `simpleName` invariants when scoped sub-workflows invoke agents; add a one-line comment in `VariableScopeVisitor.visitAgent` documenting why agent outputs are intentionally not wrapped in an output DSL scope (divergence from `visitProcessV2`); harden the resolution test to assert no unused-variable warning is emitted for an input referenced only in the `prompt:` (the `findVariableDeclaration` suppression line is currently uncovered); clarify the Elvis-precedence idiom in `ChatModelFactory.providerOf/modelOf` (`modelId?.indexOf('/') ?: -1` masks a `0` index — benign today but fragile); investigate the typed-DSL operator nuance where `agent_name | view` (bare pipe to `view`) fails to resolve under `nextflow.enable.types = true` while the chained `agent_name(ch).view { ... }` form works — confirm whether this is general to typed operators or specific to the agent `ChannelOut`, and align the docs/tests accordingly.

## Record-typed structured I/O (DELIVERED 2026-06-02)

**Status: delivered** (Plan F, [`docs/superpowers/plans/2026-06-02-agent-record-io.md`](../docs/superpowers/plans/2026-06-02-agent-record-io.md)). Verified with a real run: `./launch.sh run` of a record-typed agent (`input: q: Question` / `output: a: Answer`) produced a structured record bound from the live OpenAI structured-output path — output `CAPITAL=Paris COUNTRY=France`, log `LangChainAgentRunner - Running agent model=openai/gpt-5-mini; messages=2; structured=true`, with the `Answer` fields read individually off the emitted record.

The free-style scalar I/O used in the original POC (`input: question: String` / `output: answer: String`, emitting the LLM's raw text) was inadequate for real usage. **Agent I/O now uses record types** ([`adr/20260306-record-types.md`](20260306-record-types.md)), and the record structure drives the LLM's structured I/O:

```groovy
record Question { text: String; context: String? }
record Answer   { answer: String; confidence: Double }

agent qa {
    model 'openai/gpt-5-mini'
    instruction 'You are a careful question-answering assistant.'
    input:  q: Question
    output: a: Answer
    prompt: "Answer the question: ${q.text}"
}
```

**Output side (the core capability):** the output record type's structure is reflected into a JSON schema and used as the LLM **structured-output contract** (langchain4j `responseFormat` + OpenAI strict mode). The returned JSON is parsed and bound to a record instance emitted on the agent's output channel, flowing into a downstream typed agent/process.

**Input side:** the input record (a `RecordMap`) is bound into the `prompt:` closure (so `${q.text}` resolves) **and** serialized as a JSON object appended to the user message, so the model receives the structured input value (decision 2026-06-02).

### Decisions (locked 2026-06-02, after research + adversarial review)

- **Named record types only** for v1 inputs and outputs. Destructured `record(...)` agent I/O currently loses field metadata at lowering (input → bare `Record`, output → dropped) and is rejected at resolution with a clear error; supporting it is deferred.
- **Records required** — scalar I/O (`String`, etc.) is rejected at resolution. This is a breaking change to the experimental feature; the existing scalar-`String` tests/example/docs are migrated to records.
- **Schema derived at runtime by reflection** on the output record `Class`. The probe confirmed a named record type lowers to a concrete, field-introspectable JVM class at runtime (the brief's feared compile-time "metadata-loss bug" does not apply to named types), so no nf-lang lowering change is needed — only a resolution constraint.
- **SPI stays text-based.** `AgentRunner.run` still returns the raw JSON `String`; the derived schema (a portable `Map`) and input JSON travel in `AgentRunnerRequest`; the plugin sets `responseFormat` from the schema; **core** parses + binds the JSON to the output record (reusing `TypeHelper.asRecordType`). This keeps record/type coupling in core and the plugin thin.
- **`Path` fields forbidden in agent *outputs*** for v1 (LLMs don't produce real files; `TypeHelper.asType`'s `Path` branch existence-checks and would throw). Allowed in inputs (serialized to the absolute path string). v1 output field types: `String`/`Integer`/`Long`/`Double`/`Boolean`, nested records, and `List` of those; enum/`Map`/`Set` deferred.
- **`prompt:` stays mandatory** (grammar-enforced) and **single input / single output** for v1 (a record already aggregates many fields).

### Known semantic caveat (discovered during implementation)

Under OpenAI **strict** structured output, langchain4j force-promotes *all* output schema properties to `required` before sending. So an optional (`?`) field in an agent **output** record is effectively required — the model must always emit it. Optional fields on **input** records are honored normally. This is documented in `docs/agent.md`.

### Remaining deferrals after Plan F

Destructured `record(...)` agent I/O (rejected with a clear error; needs lowering work); `Path`/enum/`Map`/`Set` fields in agent outputs; multiple inputs/outputs; non-OpenAI providers; the tool-dispatch bridge (limitation #5 below); and the `agent {}` config scope (limitation #6).

Follow-up plan: [`docs/superpowers/plans/2026-06-02-agent-record-io.md`](../docs/superpowers/plans/2026-06-02-agent-record-io.md).

## Module-as-tool execution — validated direction (2026-06-02)

The headline goal (tool calls are real dataflow nodes; tools are modules with `meta.yml`-derived schemas) is the real point of the feature. Re-framing accordingly:

- **Agent I/O is ordinary process-style I/O** (`val`/`path`/record), not a special structured mechanism. The Plan-F LLM structured-output binding becomes **opt-in** (used only when the user declares a record output); it is no longer the centerpiece. The schema mapping that matters lives at the **tool boundary**.
- **A tool is a Nextflow module/process.** Its `meta.yml`/`ModuleSpec` `input` metadata (e.g. `nextflow module view nf-core/fastqc -o json` → input tuple `{meta:map, reads:file}`, with per-item descriptions) becomes the LLM tool's parameter JSON-schema "for free". The LLM's tool-call args are marshalled into channel values, the module runs through the standard executor/container/retry/cache machinery, and its outputs are serialized back to the LLM as JSON (files as opaque absolute-path handles).

### Execution mechanism — feasibility settled by two empirical spikes

The crux was whether a module can run as a real dataflow node driven by the (post-ignition) agent loop. Two throwaway spikes settled it under the **real async executor** (`jstack`-verified):

- **Spike #1 — naive synchronous invocation does NOT work.** Calling `ProcessDef.run()` from inside the agent operator *after* `session.fireDataflowNetwork()` deadlocks: the new operator's start is deferred via `session.addIgniter{…}` and the igniter list is drained only once at ignition, so the late operator never starts; `.val` and `session.await()` hang. **Creating a process node post-fire is out.**
- **Spike #2 — pre-wiring WORKS, no core change required.** Wire `toolOut = tool(toolIn)` *before* ignition (the tool operator starts normally), then post-fire the agent operator emits onto `toolIn` and blocks on `toolOut.val`. The task genuinely executes on the real executor (real work dir, transformed result); correlation is correct when calls are serialized; the GPars operator pool can't starve because task bodies run on a separate `execService`.

**Validated recipe** (per declared tool, at agent construction inside the workflow body, i.e. pre-ignition):
```groovy
def toolIn  = CH.queue()
def toolOut = CH.getReadChannel( toolDef.run([toolIn] as Object[])[0] )   // == toolOut = tool(toolIn)
// agent operator (post-fire) drives it; on completion:
toolIn.bind(PoisonPill.instance)   // REQUIRED — otherwise session.await() hangs forever
```

**Hard runtime obligations:** (1) the harness MUST poison every pre-wired `toolIn` once the agent finishes issuing calls, or the network never terminates; (2) tool calls are serialized (or carry explicit correlation IDs — outputs reorder under `maxForks>1`). No new core runtime primitive is needed.

### Core/plugin split

Core owns module resolution, the pre-wiring (`toolIn`/`toolOut`), the dispatch (marshal JSON args → channel value → block on output → serialize to JSON), and poison-on-complete — all langchain4j-free, exposed to the plugin via a `ToolDescriptor` (plain Maps) + a `ToolDispatcher` (`String call(name, argsJson)`) callback. The plugin owns only langchain4j: building `ToolSpecification`s from the schema Maps and running the manual tool-call loop (`ChatRequest.toolSpecifications` → `AiMessage.toolExecutionRequests()` → `ToolExecutionResultMessage`). The plugin never touches a `ProcessDef`/`Channel`/`Path`; core never imports `dev.langchain4j.*`.

### Staged scope (this is multi-plan work)

1. **Headline slice** — **DELIVERED (Plan G Phase 2, 2026-06-02).** An agent invokes a **local/in-scope process** as a tool, end-to-end with a real LLM, via pre-wire + blocking-pull + poison. The core `ModuleToolBridge`/`ToolDispatcher` derives each tool's input schema from the process's typed inputs and serializes its output back to JSON; the plugin's `LangChainAgentRunner` runs the manual tool-call loop. **Tools XOR structured output:** when `tools` are declared the LLM's free-text answer is emitted, so a record (structured) output is rejected at agent run time (`AgentDef` guard + test). **Verified with a real run** (`NXF_PLUGINS_MODE=dev ./launch.sh run` of a `uppercase`-tool agent against live OpenAI `gpt-5-mini`): the tool loop fired the dispatcher (`LangChainAgentRunner - Agent tool call name=uppercase; args={"text":"hello"}`), the `uppercase` process then ran as a real dataflow node on the executor (`Submitted process > uppercase (1)` with work dir `work/63/7bef28...`, `Task completed`), and the result flowed back to the model to produce the final answer `ANSWER=HELLO`. Registry-module tools (`nf-core/...`) remain in stage 3.
2. **Agent I/O relaxation** — allow `val`/`path` I/O; make the record structured-output opt-in (independent; can land first).
3. **Registry modules** — resolve `nf-core/fastqc` to a `ProcessDef` + `ModuleSpec`; derive the tool schema from `meta.yml` (descriptions, tuple flattening, file handles, the `meta` map, `eval`/topic outputs).
4. **Polish** — multi-tool, parallel/correlated calls, error-as-tool-result, termination edge cases, docs, gated real-LLM E2E.

Known open sub-problems (flagged for the plan): runtime module→`ProcessDef` compilation outside the normal include flow; concurrency when the agent input is a queue (multiple in-flight loops mutating `Global.session`/`ScriptMeta`); optional module inputs and the `meta`-map provenance (LLM-supplied). The earlier Plan-F deferrals (record-only I/O, etc.) are superseded by the I/O relaxation above.

Follow-up plan: [`docs/superpowers/plans/2026-06-02-agent-tool-bridge.md`](../docs/superpowers/plans/2026-06-02-agent-tool-bridge.md).

## Links

- Design spec: [`docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md`](../docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md)
- Implementation plan: [`docs/superpowers/plans/2026-05-04-agent-syntax-and-model.md`](../docs/superpowers/plans/2026-05-04-agent-syntax-and-model.md)
- Branch: `agent-syntax-and-model`
- Inspiration / inverse pattern: [`colbyford/nf-foundry-workflow`](https://github.com/colbyford/nf-foundry-workflow) — Microsoft Foundry agents calling Seqera Cloud MCP to launch Nextflow pipelines (the external-orchestrator alternative this design inverts)
- Runtime engine: [`langchain4j`](https://docs.langchain4j.dev/tutorials/agents)
- Module spec source of truth: [`adr/module-spec-schema.json`](module-spec-schema.json) (introduced in commit `571274552`)
