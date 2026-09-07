# LangChain4j agent runner plugin for Nextflow

## Summary

This plugin provides the `langchain4j` agent runner, one of the two runners that execute the Nextflow `agent` primitive.

An agent is a process-shaped primitive: it declares typed inputs and outputs, renders a prompt, calls a language model — optionally letting the model call Nextflow modules as tools — and emits the result on a channel. Each invocation runs as an ordinary Nextflow task, so work directories, retries, parallelism, resume, and lineage all apply.

The `langchain4j` runner calls the model from the driver JVM, which makes it the lighter of the two runners: no container is pulled and no runtime is shipped. It speaks the OpenAI wire protocol, so any OpenAI-compatible endpoint can be reached by keeping the `openai/` model prefix and pointing `agent.baseUrl` at it. Agents that need a shell tool or container isolation should use the `pi` runner (`nf-agent-pi`) instead.

> Agents are a preview feature. Their syntax and behavior may change in future releases.

## Get Started

Select the runner in your `nextflow.config`:

```groovy
agent.runner = 'langchain4j'
```

The plugin is loaded automatically when the configuration declares an `agent` scope. Declare it explicitly to pin a version:

```groovy
plugins {
    id 'nf-agent@0.1.0'
}
```

Provide the model credential, either through the environment:

```bash
export OPENAI_API_KEY="sk-..."
```

or through configuration, preferably via a pipeline secret:

```groovy
agent {
    runner = 'langchain4j'
    apiKey = secrets.LLM_KEY
}
```

## Examples

### A single agent

```nextflow
agent qa {
    model 'openai/gpt-5-mini'
    instruction 'You are a concise scientific assistant.'

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
    qa('What is FASTQ format?').view()
}
```

### Structured output

Declare a record type as the output. It is passed to the model as a JSON schema, and the response is validated against it:

```nextflow
record Classification {
    topic: String
    confidence: Float
}

agent classify {
    model 'openai/gpt-5-mini'

    input:
    text: String

    output:
    result: Classification

    prompt:
    """
    Classify the topic of this abstract: ${text}
    """
}
```

### An OpenAI-compatible endpoint

The `openai/` prefix names the wire protocol, not the vendor, so a compatible gateway is reached by setting `baseUrl`:

```groovy
agent {
    runner = 'langchain4j'
    baseUrl = 'https://openrouter.ai/api/v1'
    apiKey = secrets.OPENROUTER_KEY
    model = 'openai/gpt-5-mini'
}
```

### Filesystem tools

The model can be given tools, including the built-in filesystem family. On this runner the tools execute in the driver JVM, restricted to the agent's work directory, its staged `Path` inputs, and the module-output paths returned by module tools; only the work directory is writable:

```nextflow
agent reviewer {
    model 'openai/gpt-5-mini'
    tools 'fs:read', 'fs:grep'

    input:
    report: Path

    output:
    summary: String

    prompt:
    """
    Read ${report} and summarize its findings.
    """
}
```

## Resources

- [Nextflow agents documentation](https://docs.seqera.io/nextflow/agent)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
