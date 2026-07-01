# `benchmark` process directive for per-task performance metrics

- Authors: Edmund Miller
- Status: draft
- Deciders:
- Date: 2026-05-14
- Tags: directive, trace, metrics, benchmark

## Summary

Introduce a `benchmark` process directive that emits per-task runtime metrics
(wall time, peak memory, CPU usage, I/O) to a user-specified file, with an
optional `repeats` option for noise-reduction via repeated execution.

## Problem Statement

Nextflow already collects rich per-task runtime metrics through `TraceRecord` —
`realtime`, `peak_rss`, `peak_vmem`, `%cpu`, `rchar`, `wchar`, context switches,
and roughly 40 other fields — and surfaces them through the workflow-wide
`trace.txt`, the HTML `report.html`, and `timeline.html`. These artifacts are
useful for whole-pipeline post-mortem analysis but are not designed for
per-process benchmarking workflows:

- The trace file mixes all tasks across the whole run; isolating measurements
  for a single process requires post-processing.

- There is no first-class way to capture a per-task metrics artifact at a
  stable, predictable path that downstream processes can consume as input.

- There is no built-in support for repeated execution of a task to obtain
  multiple samples and reduce measurement noise — a common requirement when
  benchmarking algorithms or hardware.

- Users who want to benchmark a third-party pipeline (e.g. nf-core workflows)
  currently have no way to enable metrics capture for a specific process
  without forking the pipeline source.

Snakemake addresses the same need with a `benchmark` rule keyword that records
wall time, memory, I/O and CPU usage to a per-job TSV (or JSONL) file, with an
optional `repeat(file, N)` wrapper for multiple samples. Nextflow has no
equivalent.

## Goals

- Provide a declarative, per-process way to capture performance metrics as a
  named pipeline artifact.

- Make the directive settable from `nextflow.config` (via `withName` and
  `withLabel` selectors) so users can benchmark an existing pipeline without
  modifying its source.

- Reuse the existing `TraceRecord` collection — no new probes, no changes to
  the bash wrapper.

- Support repeated execution by spawning N real task executions and aggregating
  their metrics, not by looping inside a single task.

- Make the benchmark file a regular workflow artifact with a deterministic
  schema, usable as input to downstream processes.

- Multi-format output driven by file extension (`.tsv` by default,
  `.jsonl` opt-in).

## Non-goals

- Aggregating benchmarks across tasks — the existing `trace.txt` and
  `report.html` already cover whole-workflow summaries.

- Collecting metrics beyond what `TraceRecord` already records. Snakemake's
  `max_uss` and `max_pss` are not exposed by Nextflow's `nxf_stat` probe and
  are documented as gaps rather than added in this ADR.

- Reliable benchmarking inside `exec:` blocks. As with Snakemake's `run:`
  body, the Groovy-side execution time is intertwined with Nextflow's process
  supervision and cannot be measured cleanly.

- Changes to the trace, report, or timeline subsystems.

## Decision

Introduce a `benchmark` process directive that accepts either a path string or
a map of options. Each task writes its own metrics file at the resolved path.
When `repeats: N` is set, Nextflow submits N independent task executions and
aggregates their `TraceRecord`s into a single output file.

## Core Capabilities

### Syntax

The `benchmark` directive accepts a string shorthand or a map of options,
mirroring the `publishDir` convention:

```groovy
// process definition — shorthand
process align {
    cpus 8
    memory '16 GB'
    benchmark "benchmarks/align/${sample}.tsv"

    input:
    tuple val(sample), path(reads)

    script:
    """
    bwa mem -t ${task.cpus} ref.fa ${reads} > out.bam
    """
}
```

```groovy
// process definition — map form
process align {
    benchmark file: "benchmarks/align/${sample}.jsonl", repeats: 3
    // ...
}
```

The path is templated with task variables (`task.process`, `task.tag`,
`task.hash`, plus any input variables in scope) and resolved relative to
`workflow.launchDir`, the same as `publishDir`.

### Config-level usage

Because `benchmark` is registered in the standard process directive list, it
is settable from `nextflow.config` via the usual `withName` and `withLabel`
selectors. This is a primary use case: benchmarking an existing pipeline (for
example, an nf-core workflow) without editing its source.

```groovy
// nextflow.config — benchmark one named process, 5 repeats
process {
    withName: 'ALIGN' {
        benchmark = [file: "benchmarks/align/${task.tag}.tsv", repeats: 5]
    }
}
```

```groovy
// nextflow.config — benchmark every process matching a label
process {
    withLabel: 'heavy' {
        benchmark = "benchmarks/${task.process}/${task.tag}.tsv"
    }
}
```

```groovy
// nextflow.config — benchmark all processes (high cost; documented caveat)
process {
    benchmark = [file: "benchmarks/${task.process}/${task.tag}.jsonl", repeats: 3]
}
```

### Options

| Option | Type | Default | Description |
|--|--|--|--|
| `file` | String (templated) | _required_ | Destination path. Templated with task variables; resolved against `workflow.launchDir`. |
| `repeats` | Integer | `1` | Number of independent task executions. See _Repeat semantics_ below. |
| `fields` | List\<String\> | sensible default subset of `TraceRecord.FIELDS` | Which columns to emit. Unknown names are rejected at parse time. |
| `format` | `'tsv'` \| `'jsonl'` | inferred from `file` extension | Explicit override when the extension is ambiguous. |

### Schema

Column names match `TraceRecord.FIELDS` directly. This keeps benchmark files
consistent with `trace.txt` and `report.html` so users can correlate them
without translation.

Default columns:

`task_id`, `hash`, `process`, `tag`, `attempt`, `realtime`, `%cpu`,
`peak_rss`, `peak_vmem`, `rchar`, `wchar`, `vol_ctxt`, `inv_ctxt`

Users can override via `fields:`. Any name in `TraceRecord.FIELDS` is
accepted.

TSV files have a header row followed by one data row per repeat. JSONL files
contain one JSON object per repeat, with keys matching the column names.

### Repeat semantics

`repeats: N` does not loop inside the wrapper script. Nextflow submits N
independent task executions (each with its own workdir, its own
`.command.trace`, its own `TraceRecord`), and aggregates them into a single
benchmark file. This stays within Nextflow's existing execution model and
avoids introducing new wrapper-script behavior.

Implications:

- Each repeat is a full, isolated task run — independent scheduling, retry,
  and failure handling under the configured `errorStrategy`.

- The benchmark file is written exactly once per logical task invocation,
  with N rows (TSV) or N JSON records (JSONL), one per repeat. Each row's
  `task_id` and `hash` identifies the underlying execution so users can
  cross-reference with `trace.txt`.

- For output channels, the first successful repeat's outputs are emitted
  downstream; other repeats' workdirs are retained on disk for inspection but
  their outputs are not emitted. This avoids inventing new semantics around
  duplicate outputs and matches Snakemake's behavior where repeated runs do
  not change the rule's output graph.

- All N repeats must reach a terminal state (success or terminal failure)
  before the benchmark file is written.

### Field coverage compared to Snakemake

Snakemake's columns map cleanly onto Nextflow's native fields. The ADR
chooses to expose Nextflow names directly rather than aliasing, but the
mapping is recorded here for users porting workflows:

| Snakemake | Nextflow `TraceRecord` |
|--|--|
| `s` | `realtime` |
| `h:m:s` | derived from `realtime` |
| `max_rss` | `peak_rss` |
| `max_vms` | `peak_vmem` |
| `max_uss`, `max_pss` | not collected — gap, documented |
| `io_in` | `rchar` |
| `io_out` | `wchar` |
| `mean_load`, `cpu_time` | derived from `%cpu` and `realtime` |

### Cloud and remote executors

`TraceRecord` collection already works for AWS Batch, Google Batch, and
Kubernetes through the existing bash wrapper, so no executor-specific work
is required. The benchmark file is materialized to the resolved path using
the same path-handling code that `publishDir` uses, which supports remote
launch directories (S3, GCS, Azure Blob) without special-casing.

## Implementation surface

The change is small and self-contained:

- `modules/nextflow/src/main/groovy/nextflow/script/dsl/ProcessBuilder.groovy` —
  add `'benchmark'` to `DIRECTIVES`.

- `modules/nextflow/src/main/groovy/nextflow/script/dsl/ProcessConfigBuilder.groovy` —
  parse the string-shorthand and map forms.

- New `modules/nextflow/src/main/groovy/nextflow/processor/Benchmark.groovy` —
  mirrors `PublishDir`; fields: `file`, `repeats`, `fields`, `format`.

- `modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy` —
  after all repeats complete, serialize the collected `TraceRecord`s to
  TSV/JSONL at the resolved path.

- `modules/nextflow/src/main/groovy/nextflow/trace/TraceRecord.groovy` — no
  schema changes; `FIELDS` is reused.

- `docs/reference/process.md` — new `benchmark` section in the directive
  reference.

No edits to `command-trace.txt` or any wrapper template.

## Open questions

- Exact mechanism for fanning out N repeats: reuse the existing
  attempt/retry path (treating repeats as forced retries with success) or a
  dedicated fan-out scheduler hook. The ADR commits to the user-visible
  contract; the mechanism is left to follow-up implementation design.

- Whether to warn when a workflow-level `benchmark` directive is combined
  with caching (`cache true`) — repeated runs against a cached task will
  return cached metrics rather than fresh measurements. Tentative default:
  emit a warning the first time a benchmarked task hits a cache.

## Links

- Snakemake benchmark rules: <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#benchmark-rules>
- Related directive ADR: [hints process directive](20260323-hints-process-directive.md)
