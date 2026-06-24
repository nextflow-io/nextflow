# skesa — running an included nf-core module as a tool

An agent that calls a real **nf-core module** (`nf-core/skesa`) to assemble a
genome, with the LLM bridging the agent's input record to the module's input
schema.

## Purpose / What it demonstrates

The `tool` example exposed a trivial local process. This one exposes a real,
containerized **nf-core module** pulled in with an `include` statement — and
shows the feature that makes that practical: the tool's schema is **fetched from
the module's registry / `meta.yml` metadata**, not written by hand.

Two things to take away:

1. **`include` + `module_run` = a tool.** `include { SKESA } from 'nf-core/skesa'`
   brings the module into scope; `tools 'module_run'` then advertises it to the
   model as a tool named `SKESA`, whose parameter schema (`{meta:{id}, fastq}`,
   required) is derived automatically from the module's metadata. OpenAI
   function-calling enforces those exact field names.

2. **The LLM bridges mismatched shapes.** The agent's input record
   (`AssemblyRequest{sample_id, reads}`) looks *nothing* like the module's input
   (`{meta, fastq}`). There is no field-to-field mapping in the code — the model
   reads the prompt and the registry-derived tool schema and **synthesizes** the
   correct call, mapping `sample_id → meta.id` and `reads → fastq`. The
   `instruction` never spells out the field structure; that is the whole value of
   the feature (the schema comes from the registry, not the prompt).

## How it works

1. **`AssemblyRequest`** carries `sample_id` and `reads` (an absolute FASTQ
   path, treated as an opaque handle).

2. **The `assembler` agent** declares `tools 'module_run'` and a high-level
   `instruction` ("use the SKESA tool to assemble the provided reads, then report
   the path to the assembled contigs") — no field names.

3. **The model's tool call**, observed in the run, is e.g.
   `SKESA({"meta":{"id":"sample1"},"fastq":"/…/data/sample.fastq"})` — produced
   purely from the registry-derived schema plus the prompt's `reads` path.

4. **Execution:** SKESA runs as a real dataflow node in its container
   (provisioned by Wave), and its `fasta` output comes back to the model as JSON
   with the contigs path as an absolute handle.

5. **Plain output:** the agent output is `assembly_path: Path` (a plain type) —
   required because the agent declares `tools`.

## Key concepts

| Concept | In this example |
|---|---|
| `include` + `module_run` | An nf-core module becomes a tool named `SKESA` |
| Registry-derived schema | `{meta:{id}, fastq}` comes from `meta.yml`/registry, enforced by OpenAI |
| Record→tuple bridging | The LLM maps `sample_id`/`reads` to `meta.id`/`fastq` — no code mapping |
| No structure in the prompt | The `instruction` stays high-level; the schema carries the shape |
| Containerized module | SKESA runs in its container via Wave |
| Plain output | `assembly_path: Path` (tools ⇒ non-record output) |

## Running it

**Requirements:** an OpenAI API key, the `nf-agent` plugin, a container runtime
(Docker/Wave), and an input FASTQ at `data/sample.fastq` (the `data/` dir is
gitignored). Fetch the sarscov2 test dataset:

```bash
mkdir -p data
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
  | gunzip > data/sample.fastq

export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output (path varies):

```
ASSEMBLY=/…/work/<hash>/sample1.fa
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
