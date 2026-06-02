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
    summarizer(channel.of('What is FASTQ format?'))
        .view { answer -> answer }
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

- `input:` declares the agent's typed inputs. **v1 supports exactly one `val` input**, interpolated into the `prompt:` template via `${name}`.
- `output:` declares a single typed `val` that receives the model's final message as a string.
- `prompt:` is the templated user message sent to the model.

## Limitations (v1)

- Exactly one input and one output.
- Only the `openai` provider is supported.
- Tool dispatch is not yet implemented: `tools()` is accepted but the declared tools are not called.
- No `agent { }` configuration scope, streaming, or structured outputs yet.
