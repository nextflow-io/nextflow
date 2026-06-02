(agent-page)=

# Agents

:::{warning}
The `agent` construct is an experimental feature. Its syntax and behavior may change in future releases.
:::

An *agent* is a process-shaped primitive that wraps an LLM-driven step. Each invocation receives one input record, renders a prompt, runs a language model under a structured-output contract, and emits one output record. Agents compose with processes and other agents through the standard channel/workflow model.

Agent inputs and outputs are **named record types**. The input record is interpolated into the prompt template and also serialized as JSON and appended to the model message; the output record type is reflected into a JSON schema that constrains the model's response, and the returned JSON is parsed and bound to an output record emitted on the channel.

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

record Question {
    text: String
    context: String?
}

record Answer {
    answer: String
    confidence: Double
}

agent qa {
    model 'openai/gpt-5-mini'
    instruction 'You are a concise scientific assistant.'
    maxIterations 5

    input:
        q: Question

    output:
        a: Answer

    prompt:
    """
    Answer briefly: ${q.text}
    """
}

workflow {
    qa(channel.of(record(text: 'What is FASTQ format?')))
        .view { answer -> "${answer.answer} (confidence=${answer.confidence})" }
}
```

The `Question` record is bound under `q` in the prompt template (`${q.text}`) and also serialized as JSON and appended to the model message. The `Answer` record type is reflected into a JSON schema that the model must satisfy, so the emitted value is an `Answer` record whose fields (`answer`, `confidence`) you can read directly.

:::{note}
Under the typed DSL (`nextflow.enable.types = true`), invoke the agent with the call form and chain operators with the method/closure form (`.view { ... }`) as shown, rather than the bare pipe form (`| view`).
:::

Run it:

```bash
nextflow run main.nf
```

## Directives

| Directive | Required | Meaning |
|---|---|---|
| `model` | yes | Model identifier as `provider/model`, e.g. `openai/gpt-5-mini`. The provider prefix selects the chat-model backend. Required at run time. |
| `instruction` | no | System prompt describing the agent's role. Omitted if not set. |
| `tools` | no | Modules the agent may call as tools. *Not yet dispatched in v1 — declaring `tools()` has no effect.* |
| `maxIterations` | no | Cap on the LLM tool-calling loop (default 20). *No effect in the v1 single-shot runner.* |

Only the `prompt:` block is structurally required; `model` is additionally required when the agent runs.

## Inputs, outputs and prompt

Agent inputs and outputs must be **named record types** (declared with `record`). Scalars (e.g. `String`) and the destructured `record(...)` form are rejected at resolution.

- `input:` declares the agent's single typed input as `name: RecordType`. The input record is bound under `name` in the `prompt:` template (e.g. `${name.field}`) and is also serialized as JSON and appended to the model message so the model sees the full structured input.
- `output:` declares the agent's single typed output as `name: RecordType`. The output record type is reflected into a JSON schema that the model must satisfy (OpenAI strict structured output). The model's JSON response is parsed and bound to a record of that type, which is emitted on the channel.
- `prompt:` is the templated user message sent to the model, with the input record fields available for interpolation.

Supported output record field types are `String`, integer (`Integer`/`Long`), number (`Double`/`Float`), `Boolean`, nested records, and lists of those. `Path` fields are not allowed in output records.

:::{warning}
Under OpenAI strict structured output, langchain4j promotes **all** output properties to required. As a result, an optional (`?`) field in an agent **output** record is effectively required — the model must always emit it. Optional fields are honored on **input** records, but not on outputs.
:::

## Limitations (v1)

- Inputs and outputs must be **named record types** — scalars (e.g. `String`) are no longer allowed, and the destructured `record(...)` form is not yet supported.
- Exactly one input and one output.
- `Path` fields are not allowed in output records.
- Optional (`?`) fields on output records are effectively required under OpenAI strict structured output (see the note above). Optional input-record fields are honored.
- Only the `openai` provider is supported.
- Tool dispatch is not yet implemented: `tools()` is accepted but the declared tools are not called.
- No `agent { }` configuration scope or streaming yet.
