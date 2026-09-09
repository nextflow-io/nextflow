# isolate-triage — a real-world adaptive triage agent (capstone)

The most complete example: a goal-directed agent that uses three nf-core modules
plus the filesystem, with a **data-driven QC gate** and a **conditional**
annotation branch. The route through the tools is decided at run time, not
hard-wired in the DAG.

## Purpose / What it demonstrates

A surveillance / clinical-microbiology lab receives short reads from a bacterial
isolate. The agent is asked, in plain language, to assemble the genome, decide
whether the assembly is good enough to annotate, and only then run the
(expensive) annotation step. This pulls together everything from the earlier
examples into one adaptive pipeline:

- **`goal` + `instruction`** drive a multi-turn loop with branching logic.
- **Three nf-core modules** (`SKESA`, `ASSEMBLYSCAN`, `PROKKA`) are `include`d
  and surfaced by `nf:module_run`, each as its own registry-schema'd tool.
- **The `fs:` tools** are used to write a short triage summary JSON.
- **A QC gate decides the path:** the model reads the assembly stats and chooses
  whether to annotate — the conditional is a *reasoning step*, not a hard-coded
  edge in a DAG.

### Two kinds of tool output

This example is the clearest illustration of how module outputs reach the model:

- **Bulk / binary artifacts** the model only forwards between tools — the SKESA
  contigs, the PROKKA annotations — come back as **opaque absolute-path
  handles**. The model never reads their bytes; it just passes the path on.
- **Small text/JSON outputs** — ASSEMBLYSCAN's stats — are **inlined into the
  tool result**, so the model can *reason over the numbers* (N50, contig count)
  and gate on them. This inlining is inferred from the output file's format and
  size; you do not annotate anything for it.

## How it works

1. **`Isolate`** carries `sample_id`, `organism`, and `reads` (a FASTQ handle).
2. **The `triage` agent** declares `tools 'nf:module_run', 'fs:*'`, a `goal`
   (assemble → QC-gate → annotate-only-if-passing → report PASS/FAIL), and an
   `instruction` spelling out the QC-gate logic:
   - assemble with `SKESA`;
   - compute stats with `ASSEMBLYSCAN`;
   - **QC gate** — if the assembly is too fragmented (**N50 < 500 bp** OR
     **> 1000 contigs**) it **FAILs**: write a summary JSON via the `write` tool and
     report `FAIL …`, skipping annotation;
   - otherwise it **PASSes**: annotate with `PROKKA`, write the summary, and
     report `PASS …` with the annotation path.
   The output is a plain `verdict: String`.

3. **With the sarscov2 test data** the assembly reaches ~N50 = 310, ~6 contigs,
   which **fails** the gate — so the agent returns a `FAIL` verdict and **skips
   PROKKA**. That is the correct, deterministic outcome for this small viral read
   set.

## Exercising the PASS + PROKKA branch

To see the annotation branch run, use a real bacterial isolate FASTQ (reads that
assemble above the gate), or tighten the thresholds in the `instruction` to match
the test data. Note that **PROKKA's container is large (~6 GB)** — its first pull
takes a while, and because tool calls are serialized the run will look idle while
Docker fetches the image.

## Key concepts

| Concept | In this example |
|---|---|
| Adaptive branching | The QC gate is a reasoning step, not a DAG edge |
| `goal` + `instruction` | Objective plus the gate logic |
| Three module tools | `SKESA`, `ASSEMBLYSCAN`, `PROKKA` via `nf:module_run` |
| Path handles vs inlined output | Bulk artifacts are paths; small stats are inlined |
| `fs:*` for artifacts | Writes the triage summary JSON to the sandbox |
| Plain output | `verdict: String` (this example emits the model's text; tools *can* also return a record — see `tool-structured/`) |

## Running it

**Requirements:** an OpenAI API key, the `nf-agent-pi` plugin, a container runtime
(Docker/Wave), and an input FASTQ at `data/sample.fastq` (gitignored):

```bash
mkdir -p data
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
  | gunzip > data/sample.fastq

export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output with the sarscov2 test data:

```
TRIAGE: FAIL isolate_001: fragmented (N50=310, contigs=6), needs manual review
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
