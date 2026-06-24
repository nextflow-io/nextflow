# contig-filter — a convergence loop on real module output

This example shows an agent running a **convergence loop**: it executes the
*same* tool many times, varying one parameter and reading the metric each run,
to converge on an optimum — operating on the real output of an nf-core module.

It is the bioinformatics counterpart to [`convergence-loop`](../convergence-loop)
(which is the same loop pattern, fully offline, on a static dataset).

## The real-world task

Assembling sequencing reads produces a set of **contigs** (contiguous
reconstructed sequences) of varying lengths. A common quality metric is **N50**:
the contig length at which half the total assembly sits in contigs that long or
longer — a higher N50 means a less fragmented, "better" assembly. A standard
cleanup step is to **drop very short contigs**, but it is a trade-off:

- Raising the minimum-length filter **increases N50** (it discards the small
  contigs dragging it down)…
- …but **decreases the retained sequence** (it throws away real bases).

So the example poses an optimization: **what minimum contig length maximizes N50
while still retaining at least 70% of the assembly?**

## What the agent does

Driven only by the `goal` (no scripted steps), the model:

1. **Assembles once** — calls the `SKESA` tool (`{meta, fastq}`) → a contigs
   FASTA. This is a real stock nf-core module running in a container.
2. **Gets back the contigs path.**
3. **Loops the filter** — repeatedly calls `filter_contigs(contigs, min_len)`
   with different `min_len` values. `filter_contigs` is a small local process
   (defined in `main.nf`): it reads the FASTA, drops contigs shorter than
   `min_len`, and returns the **N50, retained length, % retained, and contig
   counts** of what remains.
4. **Steers the search** — reads the metrics from each call and picks the next
   `min_len` (a coarse scan, then refinement), **converging** on the largest
   `min_len` whose `% retained ≥ 70%`.
5. **Reports** the optimal `min_len`, the resulting N50, % retained, and the
   path to the contigs.

**Step 3–4 is the convergence loop**: the *same* tool (`filter_contigs`)
executed many times, the agent driving the parameter from the metric the
previous run returned — unlike [`goal-directed`](../goal-directed), which chains
*different* tools once (tool composition, not iteration).

## Example run

SKESA produced 6 contigs (total 1814 bp, N50 310). The agent swept
`min_len = 0, 100, 200, 500, 1000, 1500, 2000, 3000, 5000` (several deliberately
too aggressive → everything dropped), then refined `250, 300, 350, 400, 450`,
and concluded:

> N50 stays 310 for min_len 0–250 (88.9% retained); at min_len=300 N50 rises to
> 330 but retained drops below 70%, so **250** is the most stringent filter that
> maximizes N50 while retaining ≥70% of the assembly.

```
Optimal min_len: 250
Resulting N50: 310
Percent retained: 88.9% (1612 / 1814)
Contigs kept: 5 of 6
```

## Two design notes

- **Why it needs no core change.** The tunable knob `min_len` is a **declared
  input** of `filter_contigs`, and `module_run` feeds declared inputs per call.
  The custom `filter_contigs` process is exactly the "thin wrapper that exposes
  the knob as a declared input" pattern. The agent still **cannot** tune SKESA's
  *own* assembly parameters (k-mer size, etc.) — those live in Nextflow's
  `task.ext.args`, which is not part of the tool schema. Surfacing `ext.args` in
  the schema (so stock nf-core modules become tunable without a wrapper) is a
  planned follow-up.
- **Caveat: small search space.** Six contigs is a coarse landscape, so the
  optimum here is borderline — N50 barely moves within the 70% constraint. On a
  real bacterial isolate (hundreds of contigs) the N50-vs-`min_len` curve is
  smooth and the optimization is genuinely meaningful — same code, richer data.

## Running it

Needs an OpenAI key, a container runtime (Docker/Wave for SKESA), and an input
FASTQ at `data/sample.fastq` (the `data/` dir is gitignored). Fetch the sarscov2
test dataset:

```bash
mkdir -p data
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
  | gunzip > data/sample.fastq

export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

See the top-level [examples/agents/README.md](../README.md) for the dev-build
(run-from-repo) instructions.
