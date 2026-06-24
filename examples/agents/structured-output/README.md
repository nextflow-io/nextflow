# structured-output — structured output from a record-typed agent

The simplest possible agent: one LLM call with no tools, demonstrating how a
named `record` output type turns the model's response into a typed, validated
Nextflow value.

## Purpose / What it demonstrates

This example is the entry point to the `agent` primitive. It shows the one
feature that distinguishes agents from plain processes: **structured output via
record types**.

When an agent's output is a named `record`, Nextflow reflects that record's
fields and types into a JSON schema and passes it to the model as an OpenAI
structured-output contract. The model *must* return JSON that matches the
schema — no parsing heuristics, no post-processing. The returned JSON is
validated and bound to a record instance, which is emitted on the output
channel exactly like any other typed value.

Concretely, this example sends a scientific question to `gpt-5-mini` and gets
back an `Analysis` record with four fields: a plain-text summary, a numeric
confidence score, a boolean actionability flag, and a list of key points. Those
fields are immediately readable in the workflow without any string parsing.

No tools, no container, no input data file — the only requirement is an OpenAI
key. It is safe to run locally as a first test of the `nf-agent` plugin.

## How it works

1. **The `Query` input record** carries the question and an optional context
   string. Declaring `context` as `String?` makes it optional — the agent
   gracefully handles a missing value on the input side.

2. **The `Analysis` output record** has four fields:
   - `summary: String` — a concise plain-text answer
   - `confidence: Double` — the model's self-reported confidence (0–1)
   - `actionable: Boolean` — whether the answer implies a concrete next step
   - `key_points: List<String>` — a bullet-list of the main takeaways

   These four types (`String`, `Double`, `Boolean`, `List<String>`) cover the
   full set of scalar and collection types the v1 JSON-schema deriver supports.

3. **The `analyst` agent** is declared with:
   - `model 'openai/gpt-5-mini'` — the LLM to call.
   - `instruction` — a one-line system prompt that fixes the model's persona.
   - `input: query: Query` / `output: result: Analysis` — typed I/O. Because
     the output is a `record`, structured output is enabled automatically.
   - `prompt:` — the per-input message, interpolating `query.question`. The
     full `query` record is also serialized as JSON and appended to the model
     message, so the model also sees the `context` field even though the
     prompt template does not reference it explicitly.

4. **The workflow** constructs a `Query` value inline with `record(...)`,
   passes it to `analyst(channel.of(...))`, and chains `.view { }` to print
   each field of the returned `Analysis`.

   Under `nextflow.enable.types = true` the pipe operator (`|`) is replaced by
   the call form (`analyst(...)`) and `.view { }` for chaining operators.

## Key concepts

| Concept | In this example |
|---|---|
| `record` output | `result: Analysis` — enables structured output |
| JSON-schema contract | Nextflow derives the schema from `Analysis`'s fields; OpenAI enforces it |
| Optional input field | `context: String?` on `Query` |
| Typed DSL | `nextflow.enable.types = true` required for `record` syntax |
| `model` directive | `'openai/gpt-5-mini'` — provider/model string |
| `instruction` directive | System-prompt persona, set once per agent |
| `prompt:` block | Per-input message; input fields are interpolated with `${}` |
| No tools | Plain single-turn call; no `goal`, no `tools`, no `maxIterations` |

## Running it

**Requirements:**

- An OpenAI API key.
- The `nf-agent` plugin (declared in `nextflow.config`).
- No container runtime, no input data file needed.

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output (values vary by model run):

```
summary    : FASTQ is a plain-text format storing sequences and per-base quality scores.
confidence : 0.98
actionable : false
key_points : [Each record has four lines., Quality scores are Phred-encoded ASCII., ...]
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
