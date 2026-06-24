# goal-directed — let the agent plan the steps from a `goal`

You state an objective; the agent works out *which* tools to run and *in what
order* to achieve it. A declarative `goal` instead of an imperative recipe.

## Purpose / What it demonstrates

The earlier tool examples tell the agent exactly what to do. This one does the
opposite: it hands the agent a **`goal`** — a high-level statement of the
outcome — and lets the model figure out the plan.

`goal` is a directive, distinct from `instruction`. `instruction` fixes the
agent's role and constraints ("you are a QC assistant; never guess a number");
`goal` states the objective to reach ("assess whether this assembly meets the
quality bar"). The objective is folded into the system message and the model
pursues it across as many turns as it needs, deciding each tool call from what
it has learned so far. The prompt deliberately contains **no step list and names
no tools** — the model discovers the available tools (and their input schemas)
from `module_run` and assembles the plan itself.

Given the goal "assess the assembly quality," the model arrives at this plan on
its own:

1. assemble the reads with the **SKESA** tool → contigs;
2. compute statistics on those contigs with the **ASSEMBLYSCAN** tool;
3. read the N50 and contig count and compare them to the quality bar (N50 ≥ 500);
4. report a verdict.

The interesting part is step 1 → 2: the model takes SKESA's output and feeds it
as ASSEMBLYSCAN's input. That data dependency is reasoned out at run time, not
wired in code.

## Composition vs. iteration

This example **composes distinct tools** — each module runs once and the model
reasons over the results. It does **not** re-run a tool with adjusted parameters
to drive a metric. That stricter "same tool, many runs, converge on an optimum"
pattern is a different thing, shown in
[`convergence-loop`](../convergence-loop) and [`contig-filter`](../contig-filter).
Keeping the two ideas separate is exactly why this example is named
*goal-directed* rather than a "loop."

## A note on the `filesystem` capability

The agent declares `tools 'module_run', 'filesystem'`, but in this run it does
not actually need a file read: ASSEMBLYSCAN produces a *small* stats file, and
Nextflow **inlines small text/JSON tool outputs into the tool result**, so the
model reads N50 and the contig count straight from ASSEMBLYSCAN's return value.
The `filesystem` capability is there for persisting artifacts; this particular
task simply doesn't exercise it.

## How it works

- **`Sample`** carries `sample_id` and `reads` (a FASTQ path handle).
- **The `qc` agent** declares `model`, a role-only `instruction`, a `goal`
  describing the objective, `tools 'module_run', 'filesystem'`, and
  `maxIterations 15` as a safety cap. Its output is a plain `report: String`
  (the model's free-text verdict).
- `SKESA` and `ASSEMBLYSCAN` are `include`d, so `module_run` advertises each as
  its own tool with a registry-derived input schema.
- In a verified run the model called `SKESA`, then `ASSEMBLYSCAN` on the
  resulting contigs, and reported that the assembly (~N50 = 310, ~6 contigs)
  does not meet the N50 ≥ 500 bar — the correct outcome for the small sarscov2
  test set.

## Key concepts

| Concept | In this example |
|---|---|
| `goal` directive | A declared objective folded into the system message (advisory) |
| `goal` vs `instruction` | Outcome to reach vs role/constraints |
| Model-planned steps | The tool sequence isn't scripted; the model derives it |
| Data dependency at run time | SKESA's output becomes ASSEMBLYSCAN's input |
| Small outputs inlined | Stats are read from the tool result, not via a file read |
| `maxIterations` | Hard cap on how many turns the agent may take |

## Running it

**Requirements:** an OpenAI API key, the `nf-agent` plugin, a container runtime
(Docker/Wave), and an input FASTQ at `data/sample.fastq` (the `data/` dir is
gitignored):

```bash
mkdir -p data
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
  | gunzip > data/sample.fastq

export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output (the model's free-text verdict):

```
QC=… N50 ≈ 310, contigs ≈ 6 — does NOT meet the quality bar (N50 ≥ 500 bp).
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
