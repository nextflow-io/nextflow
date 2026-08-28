# Nextflow agent examples

An **agent** is a process-shaped primitive that wraps an LLM-driven step: it
takes one typed input per channel item, renders a prompt, runs a language model
(optionally calling tools in a loop), and emits one typed output. Agents compose
with processes and other agents through the normal channel/workflow model.

These examples build up from a single no-tool agent to multi-module,
tool-calling, sandboxed, and skill-equipped agents.

## A minimal agent

```nextflow
nextflow.enable.types = true

record Analysis {
    summary: String
    confidence: Float
}

agent analyst {
    model 'openai/gpt-5-mini'
    instruction 'You are a precise scientific analyst. Be concise and honest about uncertainty.'

    input:
        question: String

    output:
        result: Analysis

    prompt:
    """
    Analyze the following question and return a structured analysis.

    Question: ${question}
    """
}

workflow {
    analyst(channel.of('Is FASTQ a binary or a text format?'))
        .view { r -> "summary: ${r.summary} (confidence: ${r.confidence})" }
}
```

An `agent` looks like a process: typed `input:`/`output:`, plus a `model`, an
`instruction` (the system prompt), and a `prompt` rendered per channel item.
Here the record-typed `output:` makes the model return structured JSON. Run it
over a channel and it composes with the rest of your workflow. Everything below
builds on this shape.

## The model in one minute

- **`model`** — the LLM, as `provider/model` (e.g. `openai/gpt-5-mini`).
- **`instruction`** — the agent's role/persona (system prompt).
- **`goal`** *(optional)* — a high-level objective that steers the multi-turn
  loop; advisory (it never raises `maxIterations`). See `09_goal-directed/`,
  `10_convergence-loop/`, and `11_contig-filter/`.
- **`input:` / `output:`** — process-style typed I/O: a scalar, a `path`, or a
  named `record`. A **record output** opts into *structured output* (the record
  type becomes the model's JSON-schema contract). A plain output (e.g. `String`)
  emits the model's text.
- **`tools`** — a list of **namespaced tool refs** (not module names). Every entry
  is `family[:group]:name`; a ref that names a *non-leaf* selects its whole
  subtree, so `'nf:module_run'` means `'nf:module_run:*'`. A ref that selects
  nothing is an error, never a silent no-op.
  - **`'nf:module_run'`** — exposes **each** process in scope as its **own** tool,
    named after the module: `include`d modules **and** locally-defined
    processes. Each tool's `parameters` schema IS that module's flattened input
    schema (required fields, `additionalProperties:false`, the nf-core `meta.id`
    convention), so OpenAI function-calling enforces the field names and the
    model cannot omit or rename a field. The LLM picks which to call; it executes
    as a real dataflow node (container / work dir / cache) and its outputs return
    as JSON (files as absolute path handles; small text/JSON outputs are inlined
    so the model can reason over them). Narrow it by naming one process
    (`'nf:module_run:SKESA'`) or a group of them (`'nf:module_run:SAMTOOLS_*'`);
    the glob is case-sensitive and only ever trails.
  - **`'fs:*'`** — the sandboxed file tools — `read`, `write`, `edit`, `ls`,
    `grep` and `find` — scoped to the agent's per-invocation work directory (plus
    the output paths of modules it ran). Writes stay inside the sandbox. Take a
    subset by naming the leaves you want (`'fs:read', 'fs:write'`).
- **`skills`** *(optional)* — a list of **skills** (Anthropic-style `SKILL.md`
  folders) the agent may use. A skill packages expert instructions (plus optional
  reference files) the model loads *on demand* through runner tools (
  `activate_skill` / `read_skill_resource`) — where `tools` let the agent *run
  code*, `skills` give it *instructions*. Each entry is a **local** name (resolved
  under `skills/<name>/` beside the script) or a **remote** GitHub ref
  (`github.com/<org>/<repo>[@rev]`) cloned and cached into that same `skills/`
  directory. No code runs at inference. See `03_skills/`.

  Declaring `tools` **or `skills`** can be combined with a **record-typed
  structured output**. The Pi runner exposes the schema as a terminating
  `final_answer` tool; the JSON result binds directly to the record. If a model
  initially omits that tool, the runner gives it one corrective turn. Multiple
  named outputs use the same schema wrapper and split into N output channels.
  See `06_tool-structured/`.

## The examples (simple → complex)

| Example | Tools | Demonstrates |
|---|---|---|
| [`01_structured-output/`](01_structured-output/main.nf) | none | A single agent with **record-typed structured output** (`Query` → `Analysis`). |
| [`02_two-agents/`](02_two-agents/main.nf) | none | Two no-tool agents **chained over a channel** — stage-1's output record type is stage-2's input type. |
| [`03_skills/`](03_skills/main.nf) | `skills` | A **local skill** (`sequence-report`): the model activates the `SKILL.md` on demand and follows its instructions — emitting a fixed report format it would not otherwise produce. No module and no tool image — the agent's own runner image is the only one pulled. |
| [`04_tool/`](04_tool/main.nf) | `nf:module_run` | The simplest tool agent: an **executor-portable in-scope script process** auto-discovered and run via `nf:module_run`. |
| [`05_tool-parallel/`](05_tool-parallel/main.nf) | `nf:module_run` | A tool agent **fanning out over a queue**: one `TaskRun` per record, running in parallel, each driving its own tool calls. |
| [`06_tool-structured/`](06_tool-structured/main.nf) | `nf:module_run` | **Tools + structured output together**: the tool loop runs, then the model ends it by calling the terminating `final_answer` tool, whose arguments are the record (`Shout`). |
| [`07_module-as-tool/`](07_module-as-tool/main.nf) | `nf:module_run` | An **`include`d nf-core module** (`nf-core/skesa`) run as a tool; the LLM bridges the agent's input record to the module's tuple input. |
| [`08_filesystem/`](08_filesystem/main.nf) | `nf:module_run`, `fs:*` | Both families together, with **no tool image** (a trivial `exec:` process): run a module, then write/read a report file in the sandbox. |
| [`09_goal-directed/`](09_goal-directed/main.nf) | `nf:module_run`, `fs:*` | The **`goal` directive**: from a high-level `goal` (no step list, tools not named in the prompt) the agent plans and chains the tools itself — assemble → assembly-stats → judge against an N50 bar. Tool *composition* of distinct tools (contrast the iterative `10_convergence-loop/`). |
| [`10_convergence-loop/`](10_convergence-loop/main.nf) | `nf:module_run` | A **convergence loop** with no tool image: the agent re-runs the **same** tool many times, varying one numeric parameter, reading the metric each call, and converging (coarse scan → refine) on the threshold that maximises F1. The tunable knob is a **declared input**, so iterative tuning needs no core change. |
| [`11_contig-filter/`](11_contig-filter/main.nf) | `nf:module_run` | A **convergence loop driving real nf-core modules**: the agent titrates depth — re-running `SEQTK_SAMPLE(sample_size)` + `SKESA` with a varying subsample size to find the smallest subsample whose assembly still clears an N50 bar. `sample_size` is a declared module input, so `nf:module_run` exposes it as the tunable knob. |
| [`12_isolate-triage/`](12_isolate-triage/main.nf) | `nf:module_run`, `fs:*` | A real-world **adaptive** agent: a `goal`, three nf-core modules, a data-driven QC gate (inline JSON stats), and a conditional annotation branch. |
| [`13_samplesheet-builder/`](13_samplesheet-builder/main.nf) | none | **Agentic map-reduce**: a deterministic process *shards* (fetches raw ENA metadata per accession), a **map** agent normalizes each messy run into a samplesheet row (structured output — picks R1/R2, single/paired, drops ENA orphan files), and a **reduce** agent `collect()`s them all to reconcile replicates into one `samplesheet.csv`. Agents supply only the reasoning; the queue/shard/gather are plain Nextflow. |
| [`14_data-labelling/`](14_data-labelling/main.nf) | `skills` | **Agentic data labelling** (map-reduce): shard fetches raw ENA metadata, a **map** agent labels each sample with controlled-vocabulary ontology terms (organism/tissue/disease/assay) guided by a **skill** carrying the vocabulary (structured output — reads free text, resolves cell lines, returns `unknown` rather than hallucinate), and a **reduce** agent unifies terms across the cohort and splits confident labels from a human `review_queue`. Sibling of `13_samplesheet-builder/`, showing `skills` + structured output as the map node. |
| [`15_map-reduce/`](15_map-reduce/main.nf) | none | **Fully-agentic map-reduce**: every phase is an agent — `planner` shards a brief, `mapper` answers each shard as one **parallel** TaskRun, `reducer` fans in the whole `Bag` via `collect()`. Exercises parallel map + multiple structured I/O + fan-in + resume end to end. |
| [`16_fan-in-parity/`](16_fan-in-parity/main.nf) | none | A canonical `process` and an `agent` reduce the **same** `collect()`ed `Bag<Finding>` — both fire exactly **once** over the bag and cache identically on `-resume`, showing the agent inherits canonical Nextflow fan-in cardinality. |
| [`17_agent-module/`](17_agent-module/main.nf) | `skills`, `nf:module_run` | An **agent authored as a module** and `include`d under an alias (`reporter as qc`). The module directory bundles its own two `skills/` and its own `tools/` module: skills, relative tool paths and `moduleDir` all resolve from the **defining** module dir, and `nf:module_run` sees only the module's own scope. No network; the tool declares a container, run by the same engine the agent's runner image needs. |
| [`19_shell-tools/`](19_shell-tools/main.nf) | `fs:*`, `shell:bash` | The **runner-native** families, which are executed by the runner rather than brokered back to the driver — on `pi`, the SDK builtins inside the agent container, so no tool call crosses the RPC link. The agent computes exact assembly statistics (GC%, N50) over 15 kB of sequence with `awk`: arithmetic no model can do by reading tokens, so an agent without a shell can only estimate. Expected numbers are fixed and verifiable. `shell:bash` is **`pi`-only** — on `langchain4j` it is rejected at agent-build time, because that loop runs in the driver JVM. |

## Requirements

- An OpenAI API key:
  ```bash
  export OPENAI_API_KEY="sk-..."
  ```
  The model is called from *inside* the agent container, which does not inherit your shell —
  exporting it is nevertheless enough: Nextflow resolves the credential once in the driver and
  hands it to the agent task in band, on the TLS-protected RPC link, so nothing is forwarded
  through the container environment and no tool process an agent invokes ever sees it. Set
  `agent.rpc.tls = false` while debugging and the credential is withheld instead, and the agent
  then needs an out-of-band channel — `agent.containerOptions = '-e OPENAI_API_KEY'` on Docker or
  Podman, the `env` scope or a `secret` elsewhere. Repinning `agent.container` to an older image is
  not a way around it: the harness enforces the RPC protocol version, so anything below
  `nf-agent-pi 0.4.1` is refused at the start frame with `Unsupported protocol version` rather than
  running with a feature missing.
- The `nf-agent-pi` plugin and `agent.runner = 'pi'` (declared in each example's `nextflow.config`).
- **A container engine, and the `pi` runner image.** `nf-agent-pi` ships no runtime of its own:
  the agent proxy and the Node harness live in the runner image, so *every* example runs its agent
  as a containerized task — and nothing needs Node (or any other interpreter) on the machine
  running Nextflow. Each `nextflow.config` therefore enables Docker — and nothing more:
  **no example names an image**. The plugin declares the one it needs, derived from its own
  `plugins/nf-agent-pi/VERSION` and published by the Nextflow release;
  `plugins/nf-agent-pi/build-image.sh ref` prints it. Set `agent.container` only to override
  that — with a locally built image, for instance, which is what
  `plugins/nf-agent-pi/build-image.sh build -l` prints. Using another engine is a one-line swap
  of `docker.enabled` for `podman.enabled`, `singularity.enabled`, and so on.
- The tool examples that use `nf-core/*` modules (`07_module-as-tool/`, `09_goal-directed/`,
  `11_contig-filter/`, `12_isolate-triage/`) additionally need Wave (enabled in their configs) to
  provision the *module* images — those modules run as real containerized tasks —
  **and an input FASTQ**. Each of those four carries a `data` **symlink** to
  [`examples/data/`](../data/README.md), so the fixture is fetched **once** and shared by all
  of them rather than copied per example. From the repository root:
  ```bash
  mkdir -p examples/data
  curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
    -o examples/data/sample.fastq.gz
  gunzip -c examples/data/sample.fastq.gz > examples/data/sample.fastq
  ```
  That single command produces both forms, and both are needed: `11_contig-filter/` reads
  `data/sample.fastq.gz` because its `SEQTK_SAMPLE` step names its output after the input and
  emits `*.fastq.gz`, while the other three read the uncompressed `data/sample.fastq`.
  The fixture files are gitignored; the `examples/data/` directory and the symlinks are not,
  so a fresh clone has the wiring in place and only the download missing.
  **Note for `12_isolate-triage/`:** the sarscov2 test FASTQ assembles to ~N50=310,
  ~6 contigs, which **correctly FAILS** the lenient QC gate (N50 < 500 bp), so the
  agent returns a FAIL verdict and skips PROKKA annotation. This is the expected
  and intended gate behaviour for this small viral read set. To exercise the PASS +
  PROKKA annotation branch, use a real bacterial isolate FASTQ or adjust the
  thresholds in `isolate-triage/main.nf`.
  The `08_filesystem/` and `10_convergence-loop/` examples need no input data and no
  *tool* image: their tools are local `exec:` processes, so the runner image is the only
  one pulled. `04_tool/` uses a portable `script:` process, which needs no image of its
  own locally and can be assigned one when offloaded to a remote executor. `03_skills/`
  runs no tool at all (a pure-LLM agent with a local skill — the OpenAI key and the
  runner image are all it takes). A **remote** skill ref additionally needs network
  access to clone the GitHub repo on first use (it is then cached locally).

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
./gradlew :plugins:nf-agent-pi:jar :plugins:nf-agent-pi:copyPluginManifest :plugins:nf-agent-pi:copyPluginLibs

# then, from an example directory
NXF_PLUGINS_MODE=dev OPENAI_API_KEY="$OPENAI_API_KEY" \
  ../../../launch.sh run main.nf
```

## v1 limits / notes

- Multiple inputs and multiple named structured outputs are supported (a queue
  input maps per-item, a value/singleton input fans in); at least one output is
  required (zero outputs are not yet supported). Provider support follows the Pi SDK.
- Structured output requires a named `record` type; `Path` fields are not
  allowed in output records. Optional (`?`) fields are preserved in the portable
  schema for both input and output records.
- A tool agent runs on the **task path** like any other agent, so it fans out over
  a queue in parallel (one `TaskRun` per record — see [`05_tool-parallel/`](05_tool-parallel/main.nf)).
  Within a single agent the tool calls are **sequential** — the model waits for each
  result before asking for the next — but calls from different agent tasks are
  dispatched concurrently, each into its own cloned process graph. A tool agent
  **caches on `-resume`** like any other agent: each
  tool's schema and backing script are folded into the resume key, so a replay is
  only ever served for the exact tools it was produced with.
  A failing tool *task* aborts the run (only dispatch-level errors — bad module name,
  malformed args — are returned to the model for recovery).
- No `nextflow run` tool in this release. `shell:bash` exists on the `pi` runner
  only — the `langchain4j` loop runs in the driver JVM, so declaring it there is
  rejected before the run starts. The `tools '<module-ref>'` string form is
  **gone** — a module path or a registry ref is no longer a tool ref: `include`
  the module, then name it with `'nf:module_run:<NAME>'` (or take every process
  in scope with `'nf:module_run'`).
- **Skills** run in *Tool Mode* only — instructions plus bundled resources, no code
  execution at inference. A **remote** skill's `SKILL.md` becomes model
  instructions, so treat it as untrusted and **pin a commit SHA** rather than a
  moving branch; cached clones land in the local `skills/` dir (add them to
  `.gitignore`). Skills resolve from `skills/` beside the script.

See [`docs/agent.mdx`](../../docs/agent.mdx) for the full reference.
