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
3. **Lowered closure body is empty** — `ScriptToGroovyVisitor.visitAgent` emits `agent('name', { /* empty */ })` for the POC. The directive/I/O/prompt content lives on the `AgentNode` instance but is not propagated into the runtime closure. Once execution lands, the closure body must carry the directives so the runner can read them.
4. **Agent execution is stubbed** — `AgentDef.run()` throws `UnsupportedOperationException`. The langchain4j-driven `nf-agent` plugin (spec §3, §4) is entirely future work.
5. **Tool bridge is unimplemented** — `ModuleSpec` → langchain4j `ToolSpecification` adapter, channel-materialisation of tool args, output serialisation back to JSON.
6. **`agent` config scope** — `nextflow.config` `agent { defaultModel = ..., maxIterationsDefault = ... }` block is not parsed yet.

### Next steps

In rough order of dependency:

1. ~~**Plan B: Agent body resolution**~~ — **DONE (2026-06-01)**, see [`docs/superpowers/plans/2026-06-01-agent-body-resolution.md`](../docs/superpowers/plans/2026-06-01-agent-body-resolution.md). `AgentDsl` + `VariableScopeVisitor.visitAgent` landed; the directive-rich agent now resolves and loads. The next blocker is Plan C (closure population), since the lowered closure is still empty.
2. **Plan C: Closure population** — change `ScriptToGroovyVisitor.visitAgent` to embed the directives, inputs, outputs, and prompt into the agent's body closure (modelled on how `ProcessToGroovyVisitorV2` lowers a process). Without this, the runtime cannot read the agent's declarative content.
3. **Plan D: nf-agent plugin** — new `plugins/nf-agent/` module with `AgentRunner` SPI implementation, langchain4j `ChatModel` factory, `ToolDispatcher` that bridges JSON args ↔ channel values ↔ JSON results, and per-provider modules (`langchain4j-open-ai`, `langchain4j-anthropic`, etc.) selected via the `provider/model` directive prefix.
4. **Plan E: End-to-end demo** — a runnable example pipeline that invokes a real LLM and demonstrates tool dispatch through a real `nf-core/fastqc`-shaped module (mocked or real).
5. **Polish backlog** (any time): make `agentPrompt` grammar-optional for better errors; add `agent` config scope; add an `@author` tag style sweep across new files; revisit `simpleName` invariants when scoped sub-workflows invoke agents; add a one-line comment in `VariableScopeVisitor.visitAgent` documenting why agent outputs are intentionally not wrapped in an output DSL scope (divergence from `visitProcessV2`); harden the resolution test to assert no unused-variable warning is emitted for an input referenced only in the `prompt:` (the `findVariableDeclaration` suppression line is currently uncovered).

## Links

- Design spec: [`docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md`](../docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md)
- Implementation plan: [`docs/superpowers/plans/2026-05-04-agent-syntax-and-model.md`](../docs/superpowers/plans/2026-05-04-agent-syntax-and-model.md)
- Branch: `agent-syntax-and-model`
- Inspiration / inverse pattern: [`colbyford/nf-foundry-workflow`](https://github.com/colbyford/nf-foundry-workflow) — Microsoft Foundry agents calling Seqera Cloud MCP to launch Nextflow pipelines (the external-orchestrator alternative this design inverts)
- Runtime engine: [`langchain4j`](https://docs.langchain4j.dev/tutorials/agents)
- Module spec source of truth: [`adr/module-spec-schema.json`](module-spec-schema.json) (introduced in commit `571274552`)
