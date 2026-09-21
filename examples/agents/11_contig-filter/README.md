# contig-filter — a convergence loop driving real nf-core modules

This example shows an agent running a **convergence loop** over **real nf-core
modules**: it re-runs the *same* tools many times, varying one parameter and
reading the metric each round, to converge on an optimum.

It is the bioinformatics counterpart to [`convergence-loop`](../10_convergence-loop)
(the same loop pattern on a static local tool that needs no image) — here the tuned
parameter is a **declared input of a stock nf-core module**, so the loop turns a
knob on a real containerized tool.

## The real-world task

Sequencing usually produces **more reads than an assembly needs**. Deeper data
costs money and compute, and past a point adds little assembly quality. So a
practical question is **depth right-sizing**: how few reads can you subsample and
still get a good assembly?

Quality is measured by **N50**: the contig length at which half the total
assembly sits in contigs that long or longer — higher N50 means a less
fragmented, "better" assembly. Subsample too aggressively and the assembly
fragments (N50 falls); subsample too little and you waste data. So the example
poses an optimization: **what is the smallest read subsample whose assembly still
clears an N50 bar?**

## What the agent does

Driven only by the `goal` (no scripted steps), for each candidate subsample size
the model:

1. **Subsamples** — calls the `SEQTK_SAMPLE` tool (`{meta, reads, sample_size}`)
   with a `sample_size` fraction → a smaller FASTQ. A stock nf-core module in a
   container; `sample_size` is a **declared input**, so `nf:module_run` exposes it
   as a tunable tool parameter.
2. **Assembles** — calls the `SKESA` tool (`{meta, fastq}`) on the subsample → a
   contigs FASTA. Also a real nf-core module in a container.
3. **Measures** — calls `assembly_stats(contigs)` (a cheap local `exec:` process
   in `main.nf`) → the **N50, total length, and contig count**.
4. **Steers the search** — reads the N50 from each round and picks the next
   `sample_size` (a coarse scan, then refinement), **converging** on the
   *smallest* subsample whose assembly still has N50 ≥ 300 bp.
5. **Reports** the optimal `sample_size`, the resulting N50, contig count, and the
   path to the contigs.

**Steps 1–4 are the convergence loop**: the *same* tools (`SEQTK_SAMPLE` +
`SKESA`) executed many times, the agent driving `sample_size` from the metric the
previous round returned — unlike [`goal-directed`](../09_goal-directed), which
chains *different* tools once (tool composition, not iteration).

## Example run

On the sarscov2 test FASTQ the agent binary-searched the subsample fraction,
running `SEQTK_SAMPLE → SKESA → assembly_stats` each round:

```
fraction → N50:  0.50 → 287   0.75 → 354   0.625 → 287   0.6875 → 364   0.65625 → 364
```

> 0.65625 is the smallest tested fraction that still yields N50 ≥ 300 bp.

```
Optimal subsample fraction: 0.65625 (65.6%)
N50: 364 bp   Contigs: 7   Total length: 2424 bp
```

Exact numbers depend on the read set and the seqtk seed; the point of the demo is
the loop — real nf-core modules (`SEQTK_SAMPLE`, `SKESA`) invoked many times with
a varying declared input — not a specific verdict. On a real bacterial isolate
the depth-vs-N50 curve is smooth and the right-sizing is genuinely meaningful —
same code, richer data.

## Two design notes

- **Why it needs no core change.** The tunable knob `sample_size` is a **declared
  input** of the stock `SEQTK_SAMPLE` module, and `nf:module_run` feeds declared
  inputs per call — no wrapper needed. What the agent still **cannot** tune are a
  module's *flag-only* parameters (SKESA's k-mer size, seqtk's seed, etc.): those
  live in Nextflow's `task.ext.args`, which is not part of the tool schema.
  Surfacing `ext.args` in the schema is a planned follow-up.
- **Why `assembly_stats` is local.** N50 could come from an nf-core module
  (`assemblyscan`, as in [`goal-directed`](../09_goal-directed)), but that would
  spin a third container every iteration. Measuring is cheap and deterministic, so
  it is a small local `exec:` process — the expensive containerized work
  (`SEQTK_SAMPLE`, `SKESA`) is what the loop actually re-runs.

## Running it

Needs an OpenAI key, a container runtime (Docker/Wave for SEQTK_SAMPLE + SKESA),
and a **gzipped** input FASTQ at `data/sample.fastq.gz` (the `data/` dir is
gitignored). `SEQTK_SAMPLE` names its output after the input and emits
`*.fastq.gz`, so the input must be gzipped. Fetch the sarscov2 test dataset:

```bash
mkdir -p data
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
  -o data/sample.fastq.gz

export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

See the top-level [examples/agents/README.md](../README.md) for the dev-build
(run-from-repo) instructions.
