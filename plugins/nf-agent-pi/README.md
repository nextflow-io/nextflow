# Pi agent runner plugin for Nextflow

## Summary

This plugin provides the `pi` agent runner, one of the two runners that execute the Nextflow `agent` primitive.

An agent is a process-shaped primitive: it declares typed inputs and outputs, renders a prompt, calls a language model — optionally letting the model call Nextflow modules as tools — and emits the result on a channel. Each invocation runs as an ordinary Nextflow task, so work directories, retries, parallelism, resume, and lineage all apply.

The `pi` runner executes the agent inside a container rather than in the driver JVM. That container boundary is what makes the `shell:bash` tool available, which the `langchain4j` runner (`nf-agent`) rejects. The plugin itself ships no runtime: the agent harness and the `agent-rpc` proxy live in a runner image published alongside each Nextflow release, and the plugin jar embeds the coordinate of the image built from the same commit. Provider credentials reach the container over RPC and never enter the task environment, the task script, or the runner's credential store.

> Agents are a preview feature. Their syntax and behavior may change in future releases.

## Get Started

Select the runner in your `nextflow.config`:

```groovy
agent.runner = 'pi'
```

The plugin is loaded automatically when the configuration declares an `agent` scope. Declare it explicitly to pin a version:

```groovy
plugins {
    id 'nf-agent-pi@0.5.1'
}
```

Provide the model credential, preferably as a pipeline secret:

```groovy
agent {
    runner = 'pi'
    apiKey = secrets.LLM_KEY
    model = 'openai/gpt-5-mini'
}
```

A container engine is required, since every agent invocation runs in the runner image.

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

### Shell and filesystem tools

`shell:bash` is supported only on this runner. The runner's file tools are rooted at the work directory with the container as the outer bound; `shell:bash` has no boundary inside the container:

```nextflow
agent inspector {
    model 'openai/gpt-5-mini'
    tools 'fs:*', 'shell:bash'

    input:
    sample: Path

    output:
    report: String

    prompt:
    """
    Inspect ${sample} and report what kind of data it holds.
    """
}
```

### Pinning the runner image

By default the runner uses the image published with the Nextflow release it ships in. Override it to run a mirrored or custom image:

```groovy
agent {
    runner = 'pi'
    container = 'public.cr.seqera.io/nextflow/nf-agent-pi:0.5.1'
}
```

### TypeSafe Jev decisions

Use the `typesafe` provider with a Jev model to call TypeSafe's System One endpoint instead of a
chat model. Configure the API key through `agent.apiKey` (or `NXF_AGENT_API_KEY`):

```groovy
agent {
    runner = 'pi'
    model = 'typesafe/jev-latest'
    apiKey = secrets.TYPESAFE_API_KEY
}
```

Jev derives one typed question per agent output and evaluates all questions against the same input
state in one request. Enum outputs become Choice questions, `Double` outputs become normalized
Score questions from 0 (low) to 1 (high), and `Boolean` outputs become Noul questions (probability
at least 0.5 is `true`). For example:

```nextflow
enum Route { BILLING, TECHNICAL }
record Decision { route: Route; severity: Double; urgent: Boolean }

agent triage {
    input:
    message: String

    output:
    decision: Decision

    prompt:
    "Classify the support request, score its severity, and decide whether it is urgent: ${message}"
}
```

Set `agent.baseUrl` to override the default `https://api.typesafe.ai/v1` endpoint. A
`TYPESAFE_API_KEY` variable already present inside the runner container is also recognized.

### Reaching the driver from the container

Tool calls travel back to the driver over RPC. The driver host is inferred where possible; override it when the inferred address is not reachable from the container:

```groovy
agent {
    runner = 'pi'
    rpc {
        remoteHost = '10.0.1.7'
        port = 8099
    }
}
```

## Resources

- [Nextflow agents documentation](https://docs.seqera.io/nextflow/agent)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
