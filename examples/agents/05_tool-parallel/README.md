# tool-parallel — agent tools on the task path (parallel + structured)

A tool-using agent now lowers to a real `TaskProcessor`/`TaskRun`, exactly like a
tool-free or skills-only agent. So it **fans out over a queue in parallel**, gets
per-task work dirs and progress rows, and can pair `tools` with a **structured
(record) output**. Before this milestone a tool agent ran on a serial legacy
operator path (one item at a time); that path has been removed and `runLegacy` is
gone — every agent runs as a task.

```
seqs (queue of 4) ─► annotate (agent + `revcomp` tool, PARALLEL map) ─► SeqAnnotation record
```

## What it demonstrates

| Capability | Where |
|---|---|
| A tool agent runs as a real task (parallel fan-out, progress rows, work dirs) | `annotate` over a 4-item queue → 4 concurrent `annotate` TaskRuns |
| `tools` + a **structured record** output | `output: annotation: SeqAnnotation` (the model ends the tool loop by calling the terminating `final_answer` tool, whose arguments are the record) |
| The tool does the exact work the LLM shouldn't | `revcomp` (reverse-complement) is deterministic string surgery an LLM gets subtly wrong; the agent is instructed to **call** it |
| Multiple could-be-wrong inputs handled | the `ACGTACGTNACG` sequence (an `N`) is flagged in `note` |

## Why an agent (and not a plain `revcomp` process)?

The mechanical part (reverse-complement) is a process — deterministic, cached,
exact. The *judgment* is the agent's: validate the sequence, decide whether to
flag a non-ACGT base, carry context into a structured record. The agent delegates
the exact computation to the tool and keeps only the reasoning.

## Run

```bash
nextflow run main.nf     # needs OPENAI_API_KEY
```

Expected (order varies — the map is parallel):

```
[PROCESS 76/a83f94] annotate (1)
[PROCESS 2d/7d0065] annotate (2)
[PROCESS b8/d97c45] annotate (3)
[PROCESS 32/d92056] annotate (4)
[PROCESS 70/a98504] revcomp (1)
...
len=12  note=ok
   seq     = ATGCATTAGCCG
   revcomp = CGGCTAATGCAT
len=12  note=contains non-ACGT base(s): N
   seq     = ACGTACGTNACG
   revcomp = CGTNACGTACGT
...
[SUCCESS] completed=8 failed=0 cached=0
```

## Resume behaviour

Everything caches. Each tool's schema and backing script source are folded into the
agent's resume key, so a replay is only ever served for the exact tools it was
produced with — editing `revcomp` re-runs `annotate`. On `-resume` both the agent
tasks and the deterministic `revcomp` tool tasks are served from cache, with no
model calls:

```
$ nextflow run main.nf -resume
[SUCCESS] completed=0 failed=0 cached=8   # 4 annotate + 4 revcomp, all cached
```

Resume for an agent means **replay**, not reproduce: the stored generation is
returned rather than the model being asked again.

## Notes

- Within one agent the tool calls are sequential — the model waits for each result
  before asking for the next — but the agent *tasks* run in parallel and their
  calls are dispatched concurrently, each into its own cloned process graph.
- A failing tool task aborts the run; dispatch-level errors (bad args, unknown
  tool) are returned to the model so it can recover.
