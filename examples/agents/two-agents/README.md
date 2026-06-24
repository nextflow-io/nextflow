# two-agents — composing agents over a channel

Two no-tool agents chained end to end: the first proposes a hypothesis, the
second peer-reviews it. Shows how agents compose with the rest of a pipeline
through the normal channel/workflow model.

## Purpose / What it demonstrates

Agents are process-shaped: one typed input per channel item in, one typed output
out. That means they **compose like processes** — the output channel of one
agent can feed straight into another. This example chains two agents and
highlights the one thing that makes the chain type-safe:

> The first agent's **output record type** is the second agent's **input record
> type**. That shared record is the entire contract between them — no tools, no
> glue code, no manual parsing.

Both agents are **no-tool** agents (single-shot LLM calls, no tool-call loop) and
both use **structured output** (their output record type is reflected into a JSON
schema the model is constrained to). So structured data flows agent → agent over
a plain Nextflow channel, exactly like any other dataflow value. It builds
directly on [`structured-output`](../structured-output) (a single record-typed
agent) by wiring two of them together.

## How it works

1. **`Hypothesis`** is the shared contract — the record type that stage 1 emits
   and stage 2 consumes:
   - `statement: String` — the proposed hypothesis
   - `rationale: String` — a brief justification
   - `confidence: Double` — the author's self-reported confidence (0–1)

2. **`Review`** is the pipeline's final output:
   - `verdict: String` — e.g. "plausible", "doubtful"
   - `critique: String` — a concise peer-review note

3. **Stage 1 — `hypothesizer`**: `input: question: String` →
   `output: hypothesis: Hypothesis`. Given a question, it proposes one testable
   hypothesis as a structured `Hypothesis` record.

4. **Stage 2 — `critic`**: `input: hypothesis: Hypothesis` →
   `output: review: Review`. It reads the `Hypothesis` fields in its prompt and
   returns a structured `Review`.

5. **The workflow** passes the `ChannelOut` of `hypothesizer` straight into
   `critic`; the `Hypothesis` records flow over the channel and a `Review` is
   emitted per input question:

   ```groovy
   def hypotheses = hypothesizer(channel.of( ...questions... ))
   critic(hypotheses).view { r -> "${r.verdict}: ${r.critique}" }
   ```

   Under `nextflow.enable.types = true` the call form is used (the bare `|` pipe
   is not used in the typed DSL).

## Key concepts

| Concept | In this example |
|---|---|
| Agent composition | Two agents chained over a channel |
| Shared record contract | Stage-1 output type **==** stage-2 input type (`Hypothesis`) |
| Structured output | Each agent's output record becomes a JSON-schema contract |
| No tools | Both are single-shot calls — no `tools`, `goal`, or `maxIterations` |
| Typed DSL | `nextflow.enable.types = true`; call form, not `|` |

## Running it

**Requirements:** an OpenAI API key and the `nf-agent` plugin. No container
runtime or data file.

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output — one `Review` per input question (values vary):

```
verdict  : plausible
critique : Testable via a randomized crossover trial; control for habitual
           caffeine intake and time-of-day effects …
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
