# map-reduce — fully-agentic map-reduce (planner → mapper → reducer)

The agentic map-reduce pattern where **every phase is an `agent`**, composed with
plain Nextflow channels. This is the workflow the *agent process-parity* work
unlocks: agents now run as real tasks, so the map fans out in parallel, every call
is cached/resumable, and each step appears in the progress table and lineage.

## The shape

```
brief ──► planner ──► plan.flatMap{ it.shards } ──► mapper ──► collect() ──► reducer ──► report
          PLAN         SHARD                          MAP                     REDUCE
        (agent)      (deterministic channel op)     (agent, parallel)       (agent, fan-in)
```

- **PLAN** — `planner` turns one brief into a `Plan` with a list of independent
  research shards (one LLM call, single value input → runs once).
- **SHARD** — `plan.flatMap { it.shards }` is ordinary Groovy on a channel: the
  shards become a queue.
- **MAP** — `mapper(shards)` runs **one task per shard**, up to `maxForks` /
  `executor.cpus` **in parallel**. Each is an independent, cached, resumable task.
- **REDUCE** — `reducer(findings.collect())` receives the whole `Bag` of findings
  as one value input, so it fires **exactly once** (fan-in) and synthesises a
  single `Report`.

Contrast [`samplesheet-builder`](../13_samplesheet-builder/main.nf), where the shard
step is a deterministic *process* (fetching metadata); here the plan itself is
agentic.

## Why this needs agent-as-task

Each capability the workflow relies on is a process-parity feature the agent now has:

| Workflow line | Capability |
|---|---|
| `mapper(shards)` over a queue | **parallel map** — one `TaskRun` per shard, concurrent |
| `reducer(findings.collect())` | **fan-in** — a `Bag`-typed value input fires the reduce once |
| `planner(...)` → `Plan{shards: List<Shard>}` | **structured record output** |
| `nextflow run … -resume` | **resume** — every cached generation replays with no model call |

An agent that only reads its declared inputs and returns records (no `tools`/`skills`)
runs on the **task path**; tool/skill agents keep the legacy behaviour.

## Running it

**Requirements:** an OpenAI API key; the `nf-agent-pi` plugin (in `nextflow.config`);
the container engine and `pi` runner image every agent task needs (see
[the examples README](../README.md#requirements)). No local data.

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf                 # planner(1) + mapper(N) + reducer(1) tasks
nextflow run main.nf -resume         # every task cached — no LLM calls
nextflow run main.nf --brief 'Your own research brief here'
```

The single-output agents auto-unwrap to their channel, so `planner(...)` **is** the
`Plan` channel (no `.plan` accessor) — mirror this when adapting the example.

## Notes

- **Parallelism is configured, not free.** On the local executor, concurrency is
  gated by CPU count and drops to serial on a 1-vCPU machine. `nextflow.config`
  raises `executor.cpus` so the per-shard mappers run concurrently; the number of
  shards is decided by the planner at run time.
- **Resume replays a stored generation** (input-keyed memoization), so it is
  reproducible but can be stale if the model changes server-side. The example uses a
  floating alias (no dated `gpt-5*` snapshot is available to the Pi runner yet), so the
  cache key follows the alias rather than the concrete model behind it.
- **Output order over the shard queue is not preserved** (parallel tasks); each
  `Finding` carries its `shardId` so results correlate by key, not position.
