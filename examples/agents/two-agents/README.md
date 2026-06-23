# Two agents, no tools, connected by a channel

Demonstrates **agent composition**: two no-tool agents chained over a channel,
the output of the first feeding the input of the second.

- `hypothesizer` — `input: question: String` → `output: hypothesis: Hypothesis`. Proposes a testable hypothesis.
- `critic` — `input: hypothesis: Hypothesis` → `output: review: Review`. Peer-reviews it.

Neither agent declares `tools`, so each runs a **single-shot** LLM call (no tool-call loop). Both use **structured output** (the output record type is reflected into a JSON schema the model is constrained to).

The composition works because the first agent's **output type** (`Hypothesis`) is the second agent's **input type** — that shared record is the only contract between them. In the workflow the `ChannelOut` from `hypothesizer` is passed straight into `critic`; the `Hypothesis` records flow over the channel like any other dataflow value.

```groovy
def hypotheses = hypothesizer(channel.of( ...questions... ))
critic(hypotheses).view { r -> "${r.verdict}: ${r.critique}" }
```

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

```bash
# from the repo root
make compile
./gradlew :plugins:nf-agent:jar :plugins:nf-agent:copyPluginManifest :plugins:nf-agent:copyPluginLibs

# then, from this directory
NXF_PLUGINS_MODE=dev OPENAI_API_KEY="$OPENAI_API_KEY" \
  ../../launch.sh run main.nf
```

## Expected output

One `Review` per input question, bound from the second agent's JSON, e.g.:

```
verdict  : plausible
critique : Testable via a randomized crossover trial; control for habitual
           caffeine intake and time-of-day effects. Confidence seems high ...
```

## Notes

- This is the **no-tools** path for both agents — contrast with `../isolate-triage` (tools that run real modules) and `../main.nf` (a single record-typed agent).
- Exactly one input and one output per agent; OpenAI provider only.
- Chaining requires stage-1's output type to equal stage-2's input type.

See `docs/agent.md` for the full reference.
