# Agent Demo & Docs Implementation Plan (Plan E)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make the agent feature usable and demonstrable: a user-facing docs page with a runnable example, and a gated end-to-end test that exercises the real OpenAI integration (skipped when no API key is present). This is the final phase — after it, the no-tools agent path is wired, documented, and demonstrably runnable.

**Architecture:** No new runtime code. A `docs/agent.md` guide (added to the Sphinx toctree) documents the DSL surface, how to enable the `nf-agent` plugin, the required `OPENAI_API_KEY`, a complete example, and the v1 limitations. A gated Spock test in the plugin (`@Requires({System.getenv('OPENAI_API_KEY')})`) drives `LangChainAgentRunner` against the real OpenAI API — the one integration point that unit tests mock — and is automatically skipped in keyless environments (CI, this dev box).

**Tech Stack:** Markdown/Sphinx docs, Spock (gated integration test), langchain4j/OpenAI.

---

## Honest scope note

A real LLM call requires `OPENAI_API_KEY` + network and **cannot be executed in this keyless environment** — the gated test is therefore *written and verified-as-skipped*, not run green here. The full stack is already proven in layers by Plan D: core dataflow execution (`AgentRunIntegrationTest`, mock runner), the plugin runner's message assembly + response handling (`LangChainAgentRunnerTest`, mock `ChatModel`), and plugin discovery via the standard PF4J `@Extension`/`getPriorityExtensions` mechanism. Plan E adds the user-facing artifacts and the one test that closes the real-API gap when a key is available.

Out of scope (unchanged): tools, non-OpenAI providers, multi-input, `agent {}` config scope.

---

## Reference

- Docs toctree: `docs/index.md` — "Developing pipelines" group (lines ~73-86, currently lists `script`, `process`, `process-typed`, `workflow`, `workflow-typed`, ...). Add `agent` there.
- Existing guide pages for tone/format: `docs/process.md`, `docs/process-typed.md`.
- Gating idiom: `@Requires({System.getenv('OPENAI_API_KEY')})` (see `modules/nextflow/src/test/groovy/nextflow/datasource/SraExplorerTest.groovy:221-222` for the `@IgnoreIf`/`@Requires` pattern).
- The runner under test: `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy` (`run(AgentRunnerRequest)→String`); request ctor order `(model, instruction, prompt, maxIterations, tools)` from `nextflow.agent.AgentRunnerRequest`.
- Plugin enable syntax in config: `plugins { id 'nf-agent' }` (see other plugins' docs).

---

## Task 1: User documentation — `docs/agent.md`

**Files:**
- Create: `docs/agent.md`
- Modify: `docs/index.md` (toctree)

- [ ] **Step 1.1: Write `docs/agent.md`**

Create `docs/agent.md`. It must include, in this order: a 1-2 sentence intro; an "Enabling the plugin" section; a complete runnable example; a "Directives" table; the I/O + prompt contract; and a "Limitations (v1)" section. Use this content (adjust prose to match the house style of `docs/process.md`, but keep all sections and the example intact):

````markdown
(agent-page)=

# Agents

:::{warning}
The `agent` construct is an experimental feature. Its syntax and behavior may change in future releases.
:::

An *agent* is a process-shaped primitive that wraps an LLM-driven step. Each invocation receives one input record, renders a prompt, runs a language model, and emits one output record. Agents compose with processes and other agents through the standard channel/workflow model.

Agents require the `nf-agent` plugin and the typed DSL.

## Enabling the plugin

Add the plugin and enable typed syntax in your configuration:

```groovy
// nextflow.config
plugins {
    id 'nf-agent'
}
```

```groovy
// at the top of your script
nextflow.enable.types = true
```

The model provider reads its API key from the environment. For OpenAI:

```bash
export OPENAI_API_KEY="sk-..."
```

## Example

```groovy
nextflow.enable.types = true

agent summarizer {
    model 'openai/gpt-5-mini'
    instruction 'You are a concise scientific assistant.'
    maxIterations 5

    input:
        question: String

    output:
        answer: String

    prompt:
    """
    Answer briefly: ${question}
    """
}

workflow {
    channel.of('What is FASTQ format?')
        | summarizer
        | view
}
```

Run it:

```bash
nextflow run main.nf
```

## Directives

| Directive | Required | Meaning |
|---|---|---|
| `model` | yes | Model identifier as `provider/model`, e.g. `openai/gpt-5-mini`. The provider prefix selects the chat-model backend. |
| `instruction` | yes | System prompt describing the agent's role. |
| `tools` | yes | Modules the agent may call as tools. *Not yet dispatched in v1 — declare `tools()`.* |
| `maxIterations` | no | Cap on the LLM tool-calling loop (default 20). |

## Inputs, outputs and prompt

- `input:` declares the agent's typed inputs. **v1 supports exactly one `val` input**, interpolated into the `prompt:` template via `${name}`.
- `output:` declares a single typed `val` that receives the model's final message as a string.
- `prompt:` is the templated user message sent to the model.

## Limitations (v1)

- Exactly one input and one output.
- Only the `openai` provider is supported.
- Tool dispatch is not yet implemented: `tools()` is accepted but the declared tools are not called.
- No `agent { }` configuration scope, streaming, or structured outputs yet.
````

- [ ] **Step 1.2: Add `agent` to the docs toctree**

In `docs/index.md`, in the "Developing pipelines" `{toctree}`, add `agent` after `workflow-typed` (or after `workflow` if `workflow-typed` isn't present — match the actual file). The edit changes only that toctree block.

- [ ] **Step 1.3: Sanity-check the docs reference**

Run a grep to confirm the page is wired and there are no obvious broken local links:
```bash
grep -n "agent" docs/index.md
```
Expected: shows the new `agent` toctree entry. (A full Sphinx build is optional and may require the docs toolchain; if `docs/make-html.sh` is trivially available, run it, otherwise skip — do not block on the docs build.)

- [ ] **Step 1.4: Commit**

```bash
git add docs/agent.md docs/index.md
git commit -s -m "docs: add agent feature guide"
```

---

## Task 2: Gated end-to-end test against the real OpenAI API

**Files:**
- Create: `plugins/nf-agent/src/test/nextflow/agent/AgentEndToEndTest.groovy`

- [ ] **Step 2.1: Write the gated test**

Create `plugins/nf-agent/src/test/nextflow/agent/AgentEndToEndTest.groovy`. It is skipped unless `OPENAI_API_KEY` is set; when set, it drives the real runner against the real API and asserts a non-empty answer:

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * Licensed under the Apache License, Version 2.0 ... (full header — copy from a sibling test)
 */
package nextflow.agent

import spock.lang.Requires
import spock.lang.Specification

/**
 * End-to-end test exercising the real OpenAI integration through the
 * langchain4j runner. Skipped automatically when OPENAI_API_KEY is not set
 * (CI and keyless dev environments).
 */
@Requires({ System.getenv('OPENAI_API_KEY') })
class AgentEndToEndTest extends Specification {

    def 'should get a real answer from the model'() {
        given:
        def runner = new LangChainAgentRunner()
        def req = new AgentRunnerRequest(
            'openai/gpt-5-mini',
            'You are a terse assistant. Answer in one word.',
            'What is the capital of France?',
            5,
            [])

        when:
        def answer = runner.run(req)

        then:
        answer != null
        !answer.trim().isEmpty()
        answer.toLowerCase().contains('paris')
    }
}
```

> Adjust the model id to a currently-available OpenAI model if `gpt-5-mini` is not valid at run time — but keep it `openai/<model>`. The `contains('paris')` assertion is a reasonable smoke check; if it proves flaky against the live model, relax it to just the non-empty assertions and leave a comment. The point is to prove the real call path works when a key is present.

- [ ] **Step 2.2: Verify it COMPILES and is SKIPPED without a key**

Run:
```bash
./gradlew :plugins:nf-agent:test --tests "nextflow.agent.AgentEndToEndTest" --rerun-tasks
```
Expected: BUILD SUCCESSFUL with the test reported as **skipped** (the `@Requires` predicate is false with no `OPENAI_API_KEY`). Confirm the build did not fail and the test did not actually call the network. (Inspect the test report or `--info` to confirm it was skipped, not executed.)

- [ ] **Step 2.3: Commit**

```bash
git add plugins/nf-agent/src/test/nextflow/agent/AgentEndToEndTest.groovy
git commit -s -m "test: gated end-to-end agent test against real OpenAI API"
```

---

## Task 3: Verification

- [ ] **Step 3.1: Full plugin suite**

```bash
./gradlew :plugins:nf-agent:test --rerun-tasks
```
Expected: BUILD SUCCESSFUL — `ChatModelFactoryTest`, `LangChainAgentRunnerTest` pass; `AgentEndToEndTest` skipped (no key). Report the skipped count.

- [ ] **Step 3.2: Full agent regression (core + lang + plugin)**

```bash
./gradlew :nf-lang:test --tests "nextflow.script.parser.AgentParserTest" --tests "nextflow.script.control.*" \
          :nextflow:test --tests "nextflow.script.Agent*" --tests "nextflow.agent.*" --tests "nextflow.script.parser.v2.AgentScriptLoadingTest" \
          :plugins:nf-agent:test --rerun-tasks
```
Expected: BUILD SUCCESSFUL across all agent suites — no regressions from the docs/test additions.

---

## What this plan does NOT deliver
- A keyless proof of the *full* DSL→plugin-discovery→real-LLM path (requires `OPENAI_API_KEY`); the gated test closes this when a key is present, and the layered Plan D tests cover the stack otherwise.
- Tools, additional providers, multi-input, config scope (all still deferred).

## Self-Review
- **Spec coverage:** Delivers ADR "Plan E: End-to-end demo" — a runnable example (in docs) and a real-LLM E2E test (gated). Documentation now exists, closing the ADR's "documenting a stub is misleading" deferral.
- **Placeholder scan:** Complete content; the two ">" notes are runtime-adjustment reminders (model id validity; assertion flakiness), not deferred work.
- **Type consistency:** `new LangChainAgentRunner().run(new AgentRunnerRequest(model, instruction, prompt, maxIterations, tools))` matches the implemented runner + `@Canonical` request ctor order from Plan D. Docs example uses the exact directive/IO/prompt syntax proven by `AgentScriptLoadingTest`.
