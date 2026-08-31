# data-labelling — agentic data labelling over public metadata

Label messy public sequencing-sample metadata (from ENA/SRA/GEO) with
**controlled-vocabulary ontology terms** — organism, tissue/cell type, disease,
assay — using the classic **map-reduce** shape where only the reasoning steps are
agents and everything else is plain Nextflow.

This is the sibling of [`samplesheet-builder/`](../13_samplesheet-builder/main.nf):
same shard → map → reduce skeleton, but the payload is *semantic labelling to an
ontology* rather than file plumbing, and the ontology is injected as a **skill**.

## Purpose / What it demonstrates

Data labelling is the textbook fit for the `agent` primitive: a label is a
*structured record*, you need *one label per item in parallel*, the labels must
obey *your* controlled vocabulary, and low-confidence items should be *routed to
a human*. Each of those maps onto a feature:

```
accessions ─► META_FETCH ─► classify ──► toSortedList() ─► audit ──► harmonized.tsv
   queue       SHARD         MAP            gather           REDUCE     + review queue
             (process)      (agent)                          (agent)
             deterministic  agentic + skill                  agentic
```

| Labelling need | Feature | Where |
|---|---|---|
| A label is structured (term + ontology ID + confidence + rationale) | **record output** → JSON-schema contract | `SampleLabels` |
| Label N samples, one per item | **map agent** (one call per run) | `classify` |
| Labels must use *your* ontology, not the model's prior | **`skills`** — a `SKILL.md` carrying the vocabulary + rules | `metadata-ontology` |
| Reconcile the cohort; split confident from review-needed | **reduce agent** over the gathered list | `audit` |
| Re-labelling a cohort is expensive | **resume** replays the map from cache | see below |

## Why an agent (and not a script that calls an LLM)?

The labelling itself is exactly the part that resists clean rules:

| Messy case (all present in the example cohort) | Handled by |
|---|---|
| `library_strategy=OTHER` but the title says `ATACseq` (`SRR891268`) | MAP reads the free text → labels assay `ATAC-seq` |
| `GM12878` cell-line name → its cell type is not stated | MAP knows `GM12878` → `B lymphocyte [CL:0000236]` |
| The same cell line appears under two runs (`SRR891268`, `SRR307898`) | REDUCE unifies them to one term |
| A bare replicate title (`N61311_untreated`, `Set7KD_rep1`) evidences no tissue | MAP returns `unknown` (does not hallucinate), lowers confidence |
| Low-confidence samples must not silently enter the dataset | REDUCE puts them in a `review_queue` |

The differentiator over a plain OpenAI call: the *structured contract*, *resumable*
execution over the cohort as first-class dataflow, and the *skill-injected
vocabulary* — all composing with the rest of a Nextflow pipeline.

## How it works

1. **`META_FETCH`** (process, `exec:`) fetches one accession's ENA filereport
   JSON and emits it as a `String`. Runs locally — no container of its own.
2. **`classify`** (agent) has a `SampleLabels` **record** output → structured
   output, and declares `skills 'metadata-ontology'`. The skill (under
   `skills/metadata-ontology/SKILL.md`) carries the allowed terms, the
   `term [ONTOLOGY:ID]` format, and the labelling rules; the model activates it on
   demand, then ends the turn by calling `final_answer`, whose arguments are the
   record.
3. **`toSortedList { it.accession }`** gathers every label into one list, ordered
   by accession (stable fan-in → deterministic output + hashable reduce input).
4. **`audit`** (agent) takes the `List<SampleLabels>`, unifies equivalent terms
   across the cohort, splits confident labels from a `review_queue`, and returns a
   `HarmonizedCohort` (`tsv` + `review_queue` + `summary`). The workflow writes
   `tsv` to `harmonized.tsv`.

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

It writes `harmonized.tsv` (gitignored) and prints it plus the review queue and a
summary. Expected shape (exact terms/confidences vary by model run):

```
accession   organism                                  tissue                    disease   assay      confidence
SRR031708   Drosophila melanogaster [NCBITaxon:7227]  unknown                   none [—]  RNA-seq    0.75
SRR1039508  Homo sapiens [NCBITaxon:9606]             unknown                   unknown   RNA-seq    0.50
SRR307898   Homo sapiens [NCBITaxon:9606]             B lymphocyte [CL:0000236]  none [—] RNA-seq    0.75
SRR3192396  Homo sapiens [NCBITaxon:9606]             B lymphocyte [CL:0000236]  none [—] RNA-seq    0.80
SRR5150592  Homo sapiens [NCBITaxon:9606]             unknown                   unknown   RNA-seq    0.50
SRR891268   Homo sapiens [NCBITaxon:9606]             B lymphocyte [CL:0000236]  none [—] ATAC-seq   0.50
```

The tells that the **skill** fired: labels use the curated ontology IDs
(`CL:0000236`, `NCBITaxon:…`), `GM12878` runs resolve to `B lymphocyte`, and
`SRR891268` is labelled `ATAC-seq` from its title despite `library_strategy=OTHER`.
The **reduce** then flags `SRR5150592` (confidence 0.50) into the review queue.
Add `-with-agent-trace` to see the `activate_skill` call.

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.

## Notes / caveats

- **The map is parallel, even though it uses a skill.** A `skills` agent lowers to
  a real task (skills carry no dataflow coupling), so `classify` fans out up to
  `executor.cpus` with progress rows and resume caching — exactly like a tool-free
  map agent (cf. [`samplesheet-builder/`](../13_samplesheet-builder/main.nf),
  [`map-reduce/`](../15_map-reduce/main.nf)). Only `tools` agents still take the legacy
  serial path.
- **`-resume` caches end to end.** The map agents (`classify`, one LLM labelling
  call per sample) are served from cache with **zero LLM calls** on resume. With
  `classify` on the task path its output is a task, so the downstream `audit` reduce
  has a stable resume hash and caches too — like [`map-reduce/`](../15_map-reduce/main.nf).
  A changed `SKILL.md` invalidates the cache (skill identity is folded into the
  resume key), so an updated vocabulary correctly re-runs the labelling.
- **Skills + structured output costs +1 turn and can be lossy** (the tool/skill
  loop runs schema-free, then a structuring turn re-encodes to JSON). Keep the
  record close to the answer and the vocabulary small — as here.
- **Confidence is the gate.** The `confidence` field lets you hard-route review
  items in the workflow with a `.branch { it.confidence < 0.6 }` if you want them
  handled downstream rather than only listed by the reduce. (LLMs occasionally
  emit an out-of-range value; the `audit` agent flags anomalies it spots.)
- **Swap in your real ontology.** `skills/metadata-ontology/SKILL.md` is a small
  illustrative subset — replace it with your institution's controlled vocabulary,
  or point `skills` at a pinned remote `SKILL.md` (`github.com/<org>/<repo>@<sha>`).
- **Not perfectly reproducible.** ENA metadata can change upstream; the pinned
  accession list keeps it stable in practice. An invalid accession fails the run.
- **To change the cohort,** edit the `channel.of(...)` accession list in `main.nf`.
