# samplesheet-builder — agentic map-reduce over public metadata

Build a valid FASTQ samplesheet from messy public sequencing metadata, using the
classic **map-reduce** shape where only the reasoning steps are agents and
everything else is plain Nextflow.

## Purpose / What it demonstrates

This is the `agent` primitive used as the **map** and **reduce** nodes of a
map-reduce pipeline (the pattern from Devin's *agentic map-reduce*, which is a
natural fit for Nextflow's dataflow model):

```
accessions ──► ENA_FETCH ──► normalize ──► collect() ──► consolidate ──► samplesheet.csv
   queue       SHARD          MAP           gather        REDUCE
             (process)       (agent)                      (agent)
             deterministic   agentic                      agentic
```

- **SHARD is deterministic.** `ENA_FETCH` is an ordinary Nextflow process that
  pulls one run's raw metadata from the ENA `filereport` API. The channel is the
  work queue; the process fans out one task per accession. No LLM.
- **MAP is agentic.** `normalize` turns one run's *messy* metadata into one clean
  samplesheet candidate. This is real judgment, not regex: it decides
  single- vs paired-end, picks the correct `_1`/`_2` mate files out of ENA's
  semicolon-joined URL list (ENA sometimes lists an extra **orphan** FASTQ that
  must be dropped), and derives a filesystem-safe sample name from a free-text
  title.
- **REDUCE is agentic.** `consolidate` receives *all* candidates at once (via
  `collect()`) and reconciles them: it groups runs that belong to the same
  biological sample under one consistent `sample` name (so nf-core merges them as
  technical replicates), guarantees unique names, and flags conflicts.

The point is the division of labour: Nextflow already provides the deterministic
scaffolding that agentic map-reduce frameworks have to build by hand — the work
queue (channel), the shard (process fan-out), and the gather (`collect()`) — so
the agents only supply reasoning.

## Why an agent (and not a process)?

The MAP and REDUCE steps are exactly the parts a bioinformatician does by hand
because they resist clean rules:

| Messy case (all present in the example accessions) | Handled by |
|---|---|
| ENA lists a stray orphan `.fastq.gz` next to `_1`/`_2` (e.g. `SRR1039513`) | MAP drops it, keeps the mates |
| Mixed single-end and paired-end runs | MAP sets `fastq_2` only when paired |
| Free-text titles → safe, consistent sample names | MAP derives the name |
| Several runs are technical replicates of one biosample | REDUCE merges them under one name |
| A sample accidentally mixing single/paired runs | REDUCE flags it in `notes` |

## How it works

1. **`ENA_FETCH`** (process, `exec:`) fetches the ENA filereport JSON for one
   accession and emits it as a `String`. Runs locally — no container of its own.
2. **`normalize`** (agent) has a `Candidate` **record** output, so it uses
   structured output: the model must return JSON matching the record's fields.
   No tools, so it stays within the v1 *tools-XOR-structured* rule.
3. **`collect()`** gathers every `Candidate` into a single list — the standard
   Nextflow gather — which becomes the single input of the reduce agent.
4. **`consolidate`** (agent) takes `List<Candidate>` and returns a `Samplesheet`
   record (`csv` + `notes`). The workflow writes `csv` to `samplesheet.csv`.

## Running it

**Requirements:**

- An OpenAI API key (`export OPENAI_API_KEY="sk-..."`).
- The `nf-agent-pi` plugin (declared in `nextflow.config`).
- Network access to the ENA API and OpenAI. **No local data**, and no image for the
  fetch process — only the container engine and `pi` runner image the agent tasks
  need (see [the examples README](../README.md#requirements)).

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

It writes `samplesheet.csv` (gitignored) and prints it plus the reduce agent's
`notes`. Expected shape (exact names/URLs vary by model run):

```
sample,fastq_1,fastq_2
N61311_untreated,ftp://ftp.sra.ebi.ac.uk/.../SRR1039508_1.fastq.gz,ftp://.../SRR1039508_2.fastq.gz
...
S2_DRSC_Untreated_1,ftp://.../SRR031708.fastq.gz,
S2_DRSC_Untreated_1,ftp://.../SRR031712.fastq.gz,
```

(The two `S2_DRSC_Untreated_1` rows are technical replicates the REDUCE agent
merged under one sample name; the `SRR1039513` orphan file was dropped by MAP.)

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.

## Notes / caveats

- **The map runs in parallel.** A tool-free agent now runs as a real task, so the
  `normalize` calls fan out concurrently (one `TaskRun` per run), bounded by
  `maxForks` / `executor.cpus` — like the `ENA_FETCH` *shard*. Each call is also
  cached, so `-resume` replays without re-calling the model.
- **Not perfectly reproducible.** ENA metadata can change upstream; the pinned
  accession list keeps it stable in practice. An invalid accession makes
  `ENA_FETCH` fail the run.
- **To change the cohort,** edit the `channel.of(...)` accession list in
  `main.nf`.
