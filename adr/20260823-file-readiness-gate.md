# `FileReadinessGate` plugin extension point for head-node file reads

- Authors: Rob Syme
- Status: proposed
- Deciders: Paolo Di Tommaso, Ben Sherman, Rob Syme
- Date: 2026-08-23
- Tags: plugin, extension-point, file-staging, cold-storage

Technical Story: follow-up to [20260516-task-readiness-gate](20260516-task-readiness-gate.md); see the "doesn't work for staging files into work/stage" finding on [nextflow-io/nextflow#7151](https://github.com/nextflow-io/nextflow/pull/7151).

## Summary

Introduce a `FileReadinessGate` plugin extension point that lets a plugin make remote files readable before Nextflow reads them on the head node. `TaskReadinessGate` gates task submission, but foreign-file staging (`FilePorter`) and content-based cache hashing both read remote objects *before* a task is scheduled, so a cold Glacier input kills the run before any gate is consulted. Plugins implement a blocking `void prepare(Collection<Path>)`; core invokes it from `FilePorter.transfer()` and, for content-hashing cache modes, before `TaskHasher.compute()`.

## Problem Statement

`TaskReadinessGate` covers exactly one window: after a `TaskHandler` is created and scheduled, before the executor submits it. Testing against real Glacier-stored inputs showed that Nextflow reads remote files in two places upstream of that window, both on the head node, in `TaskProcessor.invokeTask()`:

1. **Foreign-file staging.** `resolveTaskInputs()` classifies any input whose scheme differs from the work dir's scheme as foreign (`Executor.isForeignFile`), and `session.filePorter.transfer(foreignFiles)` downloads them synchronously on the processor thread. With an `s3://` input and a local, NFS, `az://`, or `gs://` work dir — i.e. the local executor, every grid/HPC executor, and any cross-cloud run — a cold Glacier object throws `InvalidObjectState` during `GetObject`. The error is non-transient, so `FilePorter`'s retries are futile, and the run dies with a `ProcessStageException` before the `TaskHandler` exists. No `TaskReadinessGate` implementation can help, because staging strictly precedes scheduling.

2. **Content-based cache hashing.** `cache 'deep'` and `cache 'sha256'` read full input content in `TaskHasher.compute()`, which also runs before scheduling. This fails on cold objects even in the otherwise-covered S3-work-dir case. The default cache mode is unaffected (metadata only; `HeadObject` succeeds on Glacier objects).

The covered cases work precisely because the first content read happens on the worker *after* the gate: `nxf_s3_download` at task start with an S3 work dir, or lazy FUSE reads under Fusion. Without a fix for the head-node reads, the Glacier feature works on cloud batch executors but fails on local and HPC runs — a fragmentation of the feature across run types.

A third class of head-node reads exists — channel operators and script-level file APIs (`splitCsv`, `file(...).text`, etc.) that never involve a task at all. It shares the root cause but not the remedy surface; see Non-goals.

## Goals

- Close the foreign-file staging gap: a cold remote input must be restorable before `FilePorter` attempts `GetObject`, on every executor.

- Close the content-hashing gap for `cache 'deep'` / `cache 'sha256'`.

- Reuse the design language of `TaskReadinessGate`: blocking call, interrupt-based cancellation, throw-for-permanent-failure, plugin-owned timeout policy, zero cost when no plugin registers the extension.

- Keep the primitive path-keyed and context-free so future call sites (operator-level reads) can adopt it without interface changes.

## Non-goals

- **Operator-level reads (no task exists).** `splitCsv`, `countLines`, `file('s3://...').text`, a Glacier-stored samplesheet. Gating these means touching many splitter/file-API entry points or intercepting inside the S3 `FileSystemProvider` — both invasive, and severable from this change. The extension point is deliberately shaped so these call sites can be added later. Note this does not reopen the prior ADR's rejection of filesystem-level interception: that rejection was about gating *Batch task submission*, which a head-node filesystem hook cannot see; head-node reads are exactly what such a hook does see.

- **Re-validation of restore expiry.** A restore window can lapse between the gate passing and the worker reading. Accepted lifecycle behavior, unchanged by this ADR.

- **Per-process or per-path opt-out.** The gate receives paths only, no task context. Where staging previously failed fatally, the only useful policy is "restore it"; per-process policy for the scheduler-side gate remains available via `hints` on `TaskReadinessGate`.

- **Replacing or extending `TaskReadinessGate`.** The two extension points coexist and cover disjoint windows (see Rationale).

- **Core-enforced timeouts.** Same position as the prior ADR: plugins own timeout policy.

## Considered Options

- **Option A — proactive path gate invoked by `FilePorter.transfer()` (this proposal).**
- **Option B — reactive recovery hook on transfer failure.**
- **Option C — task-level pre-stage hook in `invokeTask()`.**
- **Option D — reorder `invokeTask()` so staging happens after scheduling.**

Within Option A, a sub-decision on the signature: paths-only vs. paths plus a context object (call site, process name, hints). Paths-only was chosen — see Rationale.

## Pros and Cons of the Options

### Option A — proactive path gate in `FilePorter`

A new extension point, `void prepare(Collection<Path>)`, called synchronously at the top of `FilePorter.transfer(Batch)` with the batch's source paths, before any download is submitted.

- Good, because `FilePorter.transfer()` is the single choke point every foreign staging operation passes through — one call site closes the whole staging gap on every executor.
- Good, because the plugin sees the whole batch at once: one scan, bulk `RestoreObject`, one poll loop.
- Good, because no `FilePorter` transfer-pool permits are held while a restore waits — warm downloads for other tasks are unaffected.
- Good, because the same primitive is callable from the hasher (this ADR) and from operator entry points (future work) without interface changes.
- Bad, because the processor thread blocks for the restore duration. Mitigated: `transfer()` already blocks that thread for the full download duration today, the task cannot proceed regardless, and the call is interruptible.

### Option B — reactive recovery hook on transfer failure

An extension point consulted when a `FileTransfer` fails: `boolean recover(Path, Exception)`; returning true retries the download.

- Good, because there is no proactive scan — zero cost when nothing is cold.
- Bad, because a restore-in-progress transfer holds a transfer-pool permit for hours, starving unrelated warm downloads, unless the permit handling is restructured.
- Bad, because recovery is per-file: the plugin loses batch visibility and must coalesce restores itself.
- Bad, because it entangles the recovery contract with `FilePorter`'s existing transient-error retry logic.

### Option C — task-level pre-stage hook in `invokeTask()`

A second task-keyed hook before `filePorter.transfer()`, keyed on `TaskRun` (no handler exists yet).

- Good, because gaps 1 and 2 are covered at one call site with task context (hints) available.
- Bad, because it duplicates the gate concept at a second lifecycle point with a different key type.
- Bad, because the existing handler-keyed gate must still exist for the non-foreign case — two overlapping task-level extension points.
- Bad, because it does nothing for operator-level reads, ever.

### Option D — reorder `invokeTask()` so staging happens after scheduling

Move `filePorter.transfer()` past `checkCachedOrLaunchTask()` so the existing `TaskReadinessGate` window covers it.

- Bad, because staging is load-bearing for cache hashing: `resolveTaskInputs` rewrites foreign inputs to their stage targets and the hash is computed over them, so resume semantics depend on the current order.
- Bad, because it is a large behavioral change to the task lifecycle in exchange for what one added call in `FilePorter` achieves.

Rejected without prototyping. Documented to avoid revisiting.

## Solution or decision outcome

Adopt Option A with a paths-only signature. Introduce a `FileReadinessGate` interface in `nextflow.file`, invoked from `FilePorter.transfer()` (staging) and from `TaskProcessor.invokeTask()` before `TaskHasher.compute()` when the cache mode reads content.

## Core capabilities

### Interface

```groovy
package nextflow.file

import org.pf4j.ExtensionPoint

/**
 * Plugin extension point that makes remote files readable before Nextflow
 * reads them on the head node (foreign-file staging, content-based cache
 * hashing).
 *
 * <p>Implementations are invoked synchronously on the caller's thread and may
 * block freely; returning normally means every path in the collection is
 * readable now. Throw to signal permanent failure; the error is routed through
 * the standard staging failure path ({@code ProcessStageException}).
 *
 * <p>Implementations must ignore paths they are not responsible for
 * (non-matching schemes, already-readable objects) cheaply — every foreign
 * staging batch in the run flows through this call, most of them needing no work.
 *
 * <p>Implementations must honor {@code Thread.interrupt()} so that workflow
 * abort can unblock {@code prepare} promptly. Core does not enforce a
 * wall-clock deadline — plugins own timeout policy.
 */
interface FileReadinessGate extends ExtensionPoint {
    void prepare(Collection<Path> paths) throws InterruptedException
}
```

### Contract

| Aspect | Contract |
|---|---|
| Return | Method returns normally → all paths readable now |
| Permanent failure | Throw; wrapped in `ProcessStageException` unless already a `ProcessException` |
| Cancellation | `InterruptedException` propagates; gate must honor `Thread.interrupt()` |
| Threading | Synchronous on the caller's thread (processor/operator thread); blocking is expected |
| Irrelevant paths | Must be ignored cheaply (scheme check, no network call) |
| Aggregation across gates | All gates invoked sequentially; one failure aborts the read |
| Cross-call dedup | Plugin-internal (e.g. restore state keyed by path) |

### Call site 1 — foreign-file staging

At the top of `FilePorter.transfer(Batch)`, before `submitStagingActions`, the batch's *source* paths (the foreign originals, not the stage targets) are passed to each registered gate. Only then are downloads submitted to the transfer pool.

```groovy
void transfer(Batch batch) {
    if( batch.size() ) {
        prepareForeignFiles(batch)   // invokes each FileReadinessGate with source paths
        log.trace "Stage foreign files: $batch"
        submitStagingActions(batch.actions)
    }
}
```

Multiple processor threads call `transfer()` concurrently and the same reference file can appear in many batches. Core does not dedup across calls; a plugin keys restore state by path so the second caller joins the first caller's wait rather than issuing a second restore.

### Call site 2 — content-based cache hashing

In `TaskProcessor.invokeTask()`, between `filePorter.transfer(foreignFiles)` and `new TaskHasher(task).compute()`, guarded on the task's hash mode:

```groovy
session.filePorter.transfer(foreignFiles)

if( hashModeReadsContent(task) ) {   // DEEP or SHA256 only
    final remote = task.getInputFilesMap().values().findAll { !it.fileSystem.equals(FileSystems.default) }
    fileGates.each { it.prepare(remote) }
}
final hash = new TaskHasher(task).compute()
```

The filter matters: foreign inputs were already rewritten to local stage targets by `resolveTaskInputs`, so the remote paths remaining here are exactly the same-scheme case (e.g. S3 work dir) that `transfer()` never sees. The default hash mode skips the call entirely — metadata hashing works on cold objects today and stays zero-cost.

### Error handling

A gate throwing inside `transfer()` propagates like any staging failure: non-`ProcessException` causes are wrapped in `ProcessStageException`, so `invokeTask` routes the error through `errorStrategy` exactly as a failed download does today. The hashing call site wraps the same way. `InterruptedException` propagates untouched — the thread is being torn down on session abort — and core re-asserts the interrupt flag if a gate returns after being interrupted.

### Discovery and lifecycle

- Gates are resolved via `Plugins.getExtensions(FileReadinessGate)` lazily on first use and cached for the run.
- No manager class, no executor, no future tracking. Unlike `TaskReadinessGate`, the caller's thread is *supposed* to wait here — `transfer()` already blocks it for the full download duration — so the synchronous call is the simple and honest shape.
- With no gates registered, the cost is an empty-list check per batch. Behavior is bit-identical to upstream.

### Coexistence with `TaskReadinessGate`

The two extension points cover disjoint windows and both remain necessary:

| Window | Read location | Covered by |
|---|---|---|
| Foreign staging (work dir scheme differs from input scheme) | Head node, `FilePorter` | `FileReadinessGate` |
| Content-based cache hashing | Head node, `TaskHasher` | `FileReadinessGate` |
| Same-scheme inputs, first read at task start | Worker (`nxf_s3_download`, Fusion) | `TaskReadinessGate` |

A cold-storage plugin implements both against a shared restore manager.

## Rationale & discussion

### Why paths-only, no context object

A context parameter (call site, process name, hints) was considered and dropped. Where staging previously failed fatally, per-process opt-out has no useful semantics — the alternative to restoring is a dead run — so the hints machinery has nothing to decide here. Paths-only also makes the primitive trivially reusable from call sites that have no task (the hasher already, operator reads later) without fabricating context. As a side effect, the staging path needs no per-process annotation at all, removing the adoption gap where a missed `hints` line degraded to a fatal staging error.

### Why proactive, not reactive

Restore latency dominates (hours), so proactive scanning and reactive recovery have the same wall-clock cost — the difference is structural. The reactive shape pins a `FilePorter` transfer-pool permit per cold file for the restore duration, starving warm downloads for unrelated tasks, and forces the plugin to coalesce per-file recoveries into bulk restores itself. The proactive shape does its waiting before any permit is acquired and hands the plugin the whole batch.

### Why synchronous, when `TaskReadinessGate` got a managed executor

The prior ADR's executor exists because the *scheduler* thread must never block — it multiplexes every pending task. Here the caller is a processor thread that is already committed to blocking until this specific task's inputs are staged; moving the wait to another thread would add machinery only to have the caller wait on that thread instead. The `FilePorter` precedent cited in the prior ADR's review cuts the same way: blocking work on the thread that needs the result.

### Why the hashing call site is guarded, not unconditional

Only `cache 'deep'` and `cache 'sha256'` read content before scheduling. Guarding on hash mode keeps the default path at zero added calls, and avoids a per-task `HeadObject` scan for the overwhelmingly common metadata-hash case.

### Why operator-level reads are deferred

Their natural interception points are per-operator (many call sites across `nextflow.splitter` and the file APIs) or inside the S3 filesystem provider (an nf-amazon change with its own blast radius). Both are severable: the extension point defined here is the primitive either would call. Shipping staging plus hashing now fixes every run type for task inputs; the samplesheet-in-Glacier case remains a documented limitation until a follow-up.

## Testing

- `FilePorterTest` — a recording stub gate proving gate-before-download ordering, the empty-extension fast path, exception wrapping into `ProcessStageException`, and interrupt propagation.
- `TaskProcessorTest` — the hash-mode guard: gate invoked for `deep`/`sha256` only; staged-local paths filtered out of the collection.
- Consumer plugin — unit tests with a mocked S3 client (cold object restored then `prepare` returns; restore failure throws), plus an end-to-end scenario on the local executor with an `s3://` input: the exact configuration that dies with `ProcessStageException` today.

## Links

- Extends [20260516-task-readiness-gate](20260516-task-readiness-gate.md)
- ADR discussion: [nextflow-io/nextflow#7151](https://github.com/nextflow-io/nextflow/pull/7151)
- `TaskReadinessGate` implementation: [nextflow-io/nextflow#7158](https://github.com/nextflow-io/nextflow/pull/7158)
