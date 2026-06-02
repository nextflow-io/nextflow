(agent-page)=

# Agents

:::{warning}
The `agent` construct is an experimental feature. Its syntax and behavior may change in future releases.
:::

An *agent* is a process-shaped primitive that wraps an LLM-driven step. Each invocation receives one input value, renders a prompt, runs a language model, and emits one output value. Agents compose with processes and other agents through the standard channel/workflow model.

Agent inputs and outputs use the same typed I/O as processes: a scalar (`val`-style, e.g. `String`), a `path`, or a named `record` type. The input is interpolated into the prompt template and also serialized as JSON and appended to the model message. When the **output** is a named record type, structured output is enabled: the record type is reflected into a JSON schema that constrains the model's response, and the returned JSON is parsed and bound to an output record emitted on the channel. When the output is a scalar (e.g. `String`), the model's text response is emitted verbatim.

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

agent qa {
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
    qa(channel.of('What is FASTQ format?'))
        .view { answer -> answer }
}
```

The `question` input is bound under `question` in the prompt template (`${question}`) and also serialized as JSON and appended to the model message. The output is a scalar (`String`), so the model's text response is emitted verbatim.

### Structured output (record types)

Declare a named `record` **output** type to opt into structured output. The record type is reflected into a JSON schema that the model must satisfy, and the returned JSON is parsed and bound to a record whose fields you can read directly:

```groovy
nextflow.enable.types = true

record Answer {
    answer: String
    confidence: Double
}

agent qa {
    model 'openai/gpt-5-mini'
    instruction 'You are a concise scientific assistant.'

    input:
        question: String

    output:
        a: Answer

    prompt:
    """
    Answer briefly: ${question}
    """
}

workflow {
    qa(channel.of('What is FASTQ format?'))
        .view { answer -> "${answer.answer} (confidence=${answer.confidence})" }
}
```

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

Agent inputs and outputs use process-style typed I/O: a scalar (e.g. `String`), a `path`, or a named `record` type (declared with `record`). The destructured `record(...)` form is not yet supported and is rejected at resolution.

- `input:` declares the agent's single typed input as `name: Type`. The input is bound under `name` in the `prompt:` template (e.g. `${name}` for a scalar, `${name.field}` for a record) and is also serialized as JSON and appended to the model message so the model sees the full input.
- `output:` declares the agent's single typed output as `name: Type`. When `Type` is a named record, structured output is enabled: the record type is reflected into a JSON schema that the model must satisfy (OpenAI strict structured output), and the model's JSON response is parsed and bound to a record emitted on the channel. When `Type` is a scalar (e.g. `String`), the model's text response is emitted verbatim.
- `prompt:` is the templated user message sent to the model, with the input available for interpolation.

Supported output record field types are `String`, integer (`Integer`/`Long`), number (`Double`/`Float`), `Boolean`, nested records, and lists of those. `Path` fields are not allowed in output records.

:::{warning}
Under OpenAI strict structured output, langchain4j promotes **all** output properties to required. As a result, an optional (`?`) field in an agent **output** record is effectively required — the model must always emit it. Optional fields are honored on **input** records, but not on outputs.
:::

## Limitations (v1)

- Inputs and outputs use process-style typed I/O (scalar, `path`, or named `record`); the destructured `record(...)` form is not yet supported.
- Structured output (JSON schema binding) requires a named `record` output type; scalar outputs emit the model's text verbatim.
- Exactly one input and one output.
- `Path` fields are not allowed in output records.
- Optional (`?`) fields on output records are effectively required under OpenAI strict structured output (see the note above). Optional input-record fields are honored.
- Only the `openai` provider is supported.
- Tool dispatch is not yet implemented: `tools()` is accepted but the declared tools are not called.
- No `agent { }` configuration scope or streaming yet.
