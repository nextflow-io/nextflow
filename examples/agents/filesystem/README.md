# filesystem — the `filesystem` capability alongside `module_run`

An agent that runs a tool, then **writes and reads files** in its own sandboxed
work directory. Fully offline.

## Purpose / What it demonstrates

This example adds the second built-in capability, **`filesystem`**, and shows it
working together with `module_run`:

```groovy
tools 'module_run', 'filesystem'
```

`'filesystem'` gives the agent four sandboxed file operations —
read / write / list / exists — scoped to its **per-invocation work directory**
(plus the output paths of any modules it ran). The agent can persist
intermediate artifacts, re-read them, and reference their paths in its answer,
without ever escaping the sandbox.

The whole example runs offline: the only "module" is a trivial `exec:` process
(`word_stats`), so no container runtime or input data is required — it is a
self-contained way to see the filesystem tools in action.

## How it works

1. **The `word_stats` process** takes `text: String` and returns a small JSON
   string `{"words":N,"chars":M}` (an `exec:` block). Returning JSON means the
   agent can read the numbers straight out of the tool result.

2. **The `analyst` agent** declares both capabilities and is given a
   step-by-step `instruction`:
   1. call the `word_stats` tool to get the counts;
   2. **write** them as `report.json` in the work directory (filesystem write);
   3. **read** `report.json` back to confirm it exists (filesystem read);
   4. reply with a one-line summary that includes the absolute report path.

3. **The sandbox:** the `filesystem` writes land in the agent's work dir, so the
   final answer can quote a real, existing absolute path. Writes outside the
   sandbox are rejected.

4. **Plain output:** as with any tool-using agent, the output is a plain
   `summary: String` (the model's reply), not a record.

## Key concepts

| Concept | In this example |
|---|---|
| `'filesystem'` capability | Sandboxed read / write / list / exists |
| Work-dir sandbox | Files live in the agent's per-invocation work dir |
| Combining capabilities | `tools 'module_run', 'filesystem'` together |
| JSON tool output | `word_stats` returns JSON the model reads directly |
| Fully offline | `exec:` process — no container, no input data |

## Running it

**Requirements:** an OpenAI API key and the `nf-agent` plugin. No container
runtime or data file.

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output (path varies):

```
RESULT: Stats: 9 words, 43 chars. Report written to /…/work/<hash>/report.json
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
