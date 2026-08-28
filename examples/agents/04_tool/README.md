# tool — the simplest tool-calling agent

An agent that calls a single canonical Nextflow process as a tool. The smallest
possible demonstration of the `nf:module_run` tool family.

## Purpose / What it demonstrates

The previous examples (`structured-output`, `two-agents`) are no-tool agents:
one LLM call in, one structured value out. This example introduces **tools** —
the ability for the model to *run Nextflow processes* mid-conversation and use
their results.

It shows the simplest form: a plain in-scope process (`uppercase`) is exposed to
the model as a tool via `tools 'nf:module_run'`. No `include` statement is needed —
`nf:module_run` automatically discovers every process defined in (or included into)
the script and advertises **each one as its own tool**, named after the process.
The tool's parameter schema is derived from the process's declared inputs, so
the model is constrained to call it with the right arguments.

This is the foundation every other tool example builds on. The process uses a
portable shell `script:` body, so it can run through the local executor or be
offloaded through any canonical Nextflow executor.

## How it works

1. **The `uppercase` process** takes a `text: String` input and returns
   `result: String` captured from a portable shell `script:` block. It is
   a normal Nextflow process; nothing about it is agent-specific.

2. **The `shouty` agent** declares `tools 'nf:module_run'`. At run time the harness
   discovers `uppercase` and advertises a tool named `uppercase` whose parameter
   schema is `{text: string}` (required) — taken straight from the process's
   input declaration. The model cannot misname or omit the field.

3. **The tool-call loop:** given the prompt *"uppercase the word hello"*, the
   model decides to call `uppercase({"text": "hello"})`. The harness runs the
   process as a real dataflow node (executor, work dir, caching), serializes its
   output back to the model as JSON, and the model produces its final answer.

4. **Exact scalar output:** the agent declares `answer: String`. The Pi harness
   returns it through a structured final-answer contract so explanatory prose
   cannot contaminate the task value.

## Key concepts

| Concept | In this example |
|---|---|
| `tools 'nf:module_run'` | Exposes each in-scope process as its own tool |
| Auto-discovery | No `include` needed for a locally-defined process |
| Tool schema | Derived from the process inputs (`{text}`, required) |
| Tool-call loop | The model decides when to call the tool, then replies |
| Exact scalar output | Tool-backed `String` result uses the final-answer contract |
| Executor-portable tool | `script:` process can run locally or on a remote backend |
| No input data | The shell tool only needs a basic POSIX container remotely |

## Running it

**Requirements:** an OpenAI API key and the `nf-agent-pi` plugin (in
`nextflow.config`). No data file, and no image for the tool process — only the
container engine and `pi` runner image the agent task itself needs (see
[the examples README](../README.md#requirements)).

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output:

```
ANSWER=HELLO
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
