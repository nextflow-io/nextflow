# tool-structured — tools + structured output together

A tool-calling agent that **also** returns a record-typed structured output. It
builds directly on [`tool/`](../04_tool/main.nf) by changing the output type from a
plain `String` to a `record`.

## Purpose / What it demonstrates

Earlier, declaring `tools` (or `skills`) forced a plain (non-record) output — the
tool loop emitted the model's free text. This example shows the M5 capability:
**tools and a structured record output combined**.

## How it works

1. **The `uppercase` process** is exposed as a tool via `tools 'nf:module_run'`,
   exactly as in the `tool/` example.

2. **The tool loop:** the model calls `uppercase({"text": "hello"})`, the harness
   runs the real process and feeds the JSON result back. No response format is
   forced on the loop itself.

3. **A terminating `final_answer` tool:** because the agent declares a `record`
   output (`Shout`), the harness offers one extra tool whose parameter schema *is*
   the record schema, and tells the model it is the only valid way to finish. Its
   arguments bind to the `Shout` record and are emitted on the output channel.

## Trade-offs (read before using)

- **No second pass.** The record comes out of the loop's own final tool call, so
  there is no extra request and nothing is re-encoded from free text.
- **The model has to call it.** If it answers in prose instead, the runner allows
  one corrective turn; if it still does not call the tool, the task fails.
- **Keep the schema close to the tools' results.** The model fills the record from
  what the loop produced, so a schema demanding facts no tool returned invites
  invention.

## Running it

**Requirements:** an OpenAI API key and the `nf-agent-pi` plugin (in
`nextflow.config`). No data file, and no tool image — only the container engine
and `pi` runner image the agent task needs (see
[the examples README](../README.md#requirements)).

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output (a bound `Shout` record):

```
ANSWER=[result:HELLO]
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
