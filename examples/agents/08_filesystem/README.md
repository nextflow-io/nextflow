# filesystem — the `fs:` tool family alongside `nf:module_run`

An agent that runs a tool, then **writes and reads files** in its own sandboxed
work directory. No tool image, no input data.

## Purpose / What it demonstrates

This example adds the second built-in tool family, **`fs:`**, and shows it
working together with `nf:module_run`:

```groovy
tools 'nf:module_run', 'fs:*'
```

`'fs:*'` selects the whole family — the six sandboxed file tools `read`,
`write`, `edit`, `ls`, `grep` and `find` — scoped to its **per-invocation work
directory** (plus the output paths of any modules it ran). The agent can persist
intermediate artifacts, re-read them, and reference their paths in its answer,
without ever escaping the sandbox. Naming the leaves individually
(`'fs:read', 'fs:write'`) hands over fewer.

The tool side needs nothing: the only "module" is a trivial `exec:` process
(`word_stats`), so no tool image and no input data are required — it is a
self-contained way to see the filesystem tools in action. The agent task itself
runs in the `pi` runner image, as every agent does.

## How it works

1. **The `word_stats` process** takes `text: String` and returns a small JSON
   string `{"words":N,"chars":M}` (an `exec:` block). Returning JSON means the
   agent can read the numbers straight out of the tool result.

2. **The `analyst` agent** declares both families and is given a
   step-by-step `instruction`:
   1. call the `word_stats` tool to get the counts;
   2. **write** them as `report.json` in the work directory (the `write` tool);
   3. **read** `report.json` back to confirm it is there (the `read` tool);
   4. reply with a one-line summary that includes the absolute report path.

3. **The sandbox:** the `fs:` writes land in the agent's work dir, so the
   final answer can quote a real, existing absolute path. Writes outside the
   sandbox are rejected.

4. **Plain output:** as with any tool-using agent, the output is a plain
   `summary: String` (the model's reply), not a record.

## Key concepts

| Concept | In this example |
|---|---|
| `'fs:*'` family | Sandboxed `read`, `write`, `edit`, `ls`, `grep`, `find` |
| Work-dir sandbox | Files live in the agent's per-invocation work dir |
| Combining families | `tools 'nf:module_run', 'fs:*'` together |
| JSON tool output | `word_stats` returns JSON the model reads directly |
| No tool image | `exec:` process — no tool container, no input data |

## Running it

**Requirements:** an OpenAI API key, the `nf-agent-pi` plugin, and the container
engine plus `pi` runner image every agent task needs (see
[the examples README](../README.md#requirements)). No data file.

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

## Resume

The agent caches like any other. Its sandbox reaches only the task work directory
and the files `nf:module_run` produced during the invocation — both already covered by
the resume key — so a second run replays the stored generation with no model call:

```bash
nextflow run main.nf -resume
[SUCCESS] completed=0 failed=0 cached=2
```

Add `agent.cache = false` to force a fresh generation each run.
