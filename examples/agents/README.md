# Nextflow agent examples

An **agent** is a process-shaped primitive that wraps an LLM-driven step: it
takes one typed input per channel item, renders a prompt, runs a language model
(optionally calling tools in a loop), and emits one typed output. Agents compose
with processes and other agents through the normal channel/workflow model.

These examples build up from a single no-tool agent to multi-module,
tool-calling, sandboxed agents.

## The model in one minute

- **`model`** — the LLM, as `provider/model` (e.g. `openai/gpt-5-mini`).
- **`instruction`** — the agent's role/persona (system prompt).
- **`goal`** *(optional)* — a high-level objective that steers the multi-turn
  loop; advisory (it never raises `maxIterations`). See `isolate-triage`.
- **`input:` / `output:`** — process-style typed I/O: a scalar, a `path`, or a
  named `record`. A **record output** opts into *structured output* (the record
  type becomes the model's JSON-schema contract). A plain output (e.g. `String`)
  emits the model's text.
- **`tools`** — a list of **capabilities** (not module names):
  - **`'module_run'`** — exposes **each** process in scope as its **own** tool,
    named after the module: `include`d modules **and** locally-defined
    processes. Each tool's `parameters` schema IS that module's flattened input
    schema (required fields, `additionalProperties:false`, the nf-core `meta.id`
    convention), so OpenAI function-calling enforces the field names and the
    model cannot omit or rename a field. The LLM picks which to call; it executes
    as a real dataflow node (container / work dir / cache) and its outputs return
    as JSON (files as absolute path handles; small text/JSON outputs are inlined
    so the model can reason over them).
  - **`'filesystem'`** — a sandboxed read/write/list/exists tool scoped to the
    agent's per-invocation work directory (plus the output paths of modules it
    ran). Writes stay inside the sandbox.

  Declaring `tools` requires a **plain** (non-record) output — tools and
  structured record output are mutually exclusive in v1.

## The examples (simple → complex)

| Example | Tools | Demonstrates |
|---|---|---|
| [`main.nf`](main.nf) | none | A single agent with **record-typed structured output** (`Query` → `Analysis`). |
| [`two-agents/`](two-agents/main.nf) | none | Two no-tool agents **chained over a channel** — stage-1's output record type is stage-2's input type. |
| [`tool/`](tool/main.nf) | `module_run` | The simplest tool agent: a **local in-scope process** auto-discovered and run via `module_run`. |
| [`skesa/`](skesa/main.nf) | `module_run` | An **`include`d nf-core module** (`nf-core/skesa`) run as a tool; the LLM bridges the agent's input record to the module's tuple input. |
| [`filesystem/`](filesystem/main.nf) | `module_run`, `filesystem` | Both capabilities together, **fully offline** (a trivial `exec:` process): run a module, then write/read a report file in the sandbox. |
| [`isolate-triage/`](isolate-triage/main.nf) | `module_run`, `filesystem` | A real-world **adaptive** agent: a `goal`, three nf-core modules, a data-driven QC gate (inline JSON stats), and a conditional annotation branch. |

## Requirements

- An OpenAI API key:
  ```bash
  export OPENAI_API_KEY="sk-..."
  ```
- The `nf-agent` plugin (declared in each example's `nextflow.config`).
- The tool examples that use `nf-core/*` modules (`skesa`, `isolate-triage`)
  also need a container runtime (Docker/Wave) — the modules run as real
  containerized tasks — **and an input FASTQ at `<example>/data/sample.fastq`**
  (the `data/` dir is gitignored). Fetch the sarscov2 test dataset used for demos:
  ```bash
  mkdir -p data
  curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
    | gunzip > data/sample.fastq
  ```
  **Note for `isolate-triage`:** the sarscov2 test FASTQ assembles to ~N50=310,
  ~6 contigs, which **correctly FAILS** the lenient QC gate (N50 < 500 bp), so the
  agent returns a FAIL verdict and skips PROKKA annotation. This is the expected
  and intended gate behaviour for this small viral read set. To exercise the PASS +
  PROKKA annotation branch, use a real bacterial isolate FASTQ or adjust the
  thresholds in `isolate-triage/main.nf`.
  The `tool/` and `filesystem/` examples use `exec:` processes and run **fully
  offline** (no container, no input data).

## Run (released Nextflow)

```bash
cd <example-dir>
nextflow run main.nf
```

## Run from this repo (development build)

The plugin must be built and discovered in dev mode:

```bash
# from the repo root
make compile
./gradlew :plugins:nf-agent:jar :plugins:nf-agent:copyPluginManifest :plugins:nf-agent:copyPluginLibs

# then, from an example directory
NXF_PLUGINS_MODE=dev OPENAI_API_KEY="$OPENAI_API_KEY" \
  ../../launch.sh run main.nf
```

## v1 limits / notes

- Exactly one input and one output per agent; OpenAI provider only.
- Structured (record) output requires a named `record` type; `Path` fields are
  not allowed in output records; under OpenAI strict mode, optional (`?`) fields
  in an **output** record are effectively required (optional **input** fields are
  honored).
- `module_run` tool calls are **serialized** (one runs to completion before the
  next); a failing tool *task* aborts the run (only dispatch-level errors — bad
  module name, malformed args — are returned to the model for recovery).
- `bash` / `nextflow run` tools are not available in this release. The legacy
  `tools '<module-ref>'` string form still works but is deprecated in favor of
  `include` + `'module_run'`.

See [`docs/agent.mdx`](../../docs/agent.mdx) for the full reference.
