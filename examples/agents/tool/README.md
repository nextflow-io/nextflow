# tool — the simplest tool-calling agent

An agent that calls a single local Nextflow process as a tool. The smallest
possible demonstration of the `module_run` capability.

## Purpose / What it demonstrates

The previous examples (`structured-output`, `two-agents`) are no-tool agents:
one LLM call in, one structured value out. This example introduces **tools** —
the ability for the model to *run Nextflow processes* mid-conversation and use
their results.

It shows the simplest form: a plain in-scope process (`uppercase`) is exposed to
the model as a tool via `tools 'module_run'`. No `include` statement is needed —
`module_run` automatically discovers every process defined in (or included into)
the script and advertises **each one as its own tool**, named after the process.
The tool's parameter schema is derived from the process's declared inputs, so
the model is constrained to call it with the right arguments.

This is the foundation every other tool example builds on. It runs fully offline
(the process is a trivial `exec:` block — no container, no data file).

## How it works

1. **The `uppercase` process** takes a `text: String` input and returns
   `result: String` (computed in an `exec:` block — `text.toUpperCase()`). It is
   a normal Nextflow process; nothing about it is agent-specific.

2. **The `shouty` agent** declares `tools 'module_run'`. At run time the harness
   discovers `uppercase` and advertises a tool named `uppercase` whose parameter
   schema is `{text: string}` (required) — taken straight from the process's
   input declaration. The model cannot misname or omit the field.

3. **The tool-call loop:** given the prompt *"uppercase the word hello"*, the
   model decides to call `uppercase({"text": "hello"})`. The harness runs the
   process as a real dataflow node (executor, work dir, caching), serializes its
   output back to the model as JSON, and the model produces its final answer.

4. **Plain output:** because the agent declares `tools`, its output must be a
   plain type (`answer: String`) — the model's free-text reply. Structured
   `record` output and tools are mutually exclusive.

## Key concepts

| Concept | In this example |
|---|---|
| `tools 'module_run'` | Exposes each in-scope process as its own tool |
| Auto-discovery | No `include` needed for a locally-defined process |
| Tool schema | Derived from the process inputs (`{text}`, required) |
| Tool-call loop | The model decides when to call the tool, then replies |
| Plain output required | `tools` ⇒ non-record output (`String`) |
| Fully offline | `exec:` process — no container, no input data |

## Running it

**Requirements:** an OpenAI API key and the `nf-agent` plugin (in
`nextflow.config`). No container runtime or data file.

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
