# Record-typed agent example

Demonstrates the `agent` primitive with **record-typed structured I/O**:

- `input: query: Query` — the input record is interpolated into the prompt (`${query.question}`) and also serialized as JSON and appended to the model message.
- `output: result: Analysis` — the output record type is reflected into a JSON schema used as the LLM's structured-output contract; the returned JSON is parsed and bound to an `Analysis` record, whose fields (`summary`, `confidence`, `actionable`, `key_points`) are read directly.

## Requirements

- An OpenAI API key:
  ```bash
  export OPENAI_API_KEY="sk-..."
  ```
- The `nf-agent` plugin (declared in `nextflow.config`).

## Run (released Nextflow)

```bash
nextflow run main.nf
```

## Run from this repo (development build)

The plugin must be built and discovered in dev mode:

```bash
# from the repo root
make compile
./gradlew :plugins:nf-agent:jar :plugins:nf-agent:copyPluginManifest :plugins:nf-agent:copyPluginLibs

# then, from this directory
NXF_PLUGINS_MODE=dev OPENAI_API_KEY="$OPENAI_API_KEY" \
  ../../launch.sh run main.nf
```

## Expected output

A structured analysis bound from the model's JSON, e.g.:

```
summary    : FASTQ is a plain-text format ...
confidence : 0.95
actionable : false
key_points : [text-based, four lines per record, ...]
```

## Notes (v1 limits)

- Inputs and outputs must be **named record types** (scalars and destructured `record(...)` are rejected).
- Exactly one input and one output; OpenAI provider only.
- `Path` fields are not allowed in output records.
- Under OpenAI strict structured output, optional (`?`) fields in the **output** record are effectively required (the model must emit them). Optional **input** fields are honored.

See `docs/agent.md` for the full reference.
