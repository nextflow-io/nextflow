# `TaskReadinessGate` plugin extension point for deferred task submission

- Authors: Rob Syme
- Status: draft
- Deciders: Paolo Di Tommaso, Ben Sherman, Rob Syme
- Date: 2026-05-16
- Tags: plugin, extension-point, scheduler, task-submission

## Summary

Introduce a `TaskReadinessGate` plugin extension point that lets a plugin defer task submission until an external precondition is met (e.g. an S3 object has been restored from Glacier), removing the need for plugins to subclass an executor and its task handler.

## Problem Statement

Some plugins need to hold task submission until an external precondition is satisfied. A concrete example is restoring S3 objects from Glacier: a plugin needs to issue a `RestoreObject` request and wait for the object to be available before an AWS Batch job tries to stage it. The same shape applies to other cold-storage backends and to any "warm-up before run" workflow.

The only existing hook is `TaskHandler.isReady()`, which is consulted by `TaskPollingMonitor.canSubmit()` on every polling tick. To override it, a plugin must subclass both an executor and its task handler. This shape has three problems:

1. **Opt-in cost for users.** Pipelines must set `process.executor = '<plugin-specific-name>'`, instead of the plugin being a drop-in via `plugins { id '...' }`.

2. **Coupling to executor internals.** Subclassing an executor's task handler tangles plugin code with executor implementation details. In practice, plugins that subclass `AwsBatchTaskHandler` cannot use `@CompileStatic` because the parent's proxy dispatch breaks under static compilation — i.e. the subclassing approach is already hitting its limits for a single executor.

3. **Single-cloud lock-in.** The same shape recurs for any executor that stages remote inputs (AWS Batch, Google Batch, Azure Batch, Kubernetes, local) and any cold-storage backend (Glacier, GCS Coldline, on-demand warm-up, lazy datasets). Each combination would need a new executor subclass under the current design.

A generic extension point — "consult these plugins before submitting any task" — addresses all three problems with a small additive change to core.

## Goals

- Provide a way for plugins to defer task submission based on external state, without subclassing executors or task handlers.

- Work uniformly across every executor that submits tasks via `TaskPollingMonitor` (i.e. all existing executors).

- Be safe-by-default: a plugin that ships a gate cannot accidentally stall the scheduler thread for the rest of the run.

- Allow gates to kick off async preparation work the moment a task is queued, not when an executor slot frees up — cold-storage restores can take hours and must not wait for capacity.

- Allow plugins to signal permanent task failure cleanly, integrating with the standard `errorStrategy` machinery.

## Non-goals

- **Per-process gate scoping.** A gate is consulted for every task. Selectors that restrict a gate to specific processes are a config-layer concern that can be added later.

- **Gate ordering or priority.** Evaluation order across gates is unspecified. An `int order()` default method can be added non-breakingly if a real use case appears.

- **An `AbstractAsyncReadinessGate` helper class.** A base class that runs blocking work on a dedicated executor and exposes results via per-task `CompletableFuture`s would lower the cliff for plugin authors. We are deliberately not shipping it in v1: the contract is documented in javadoc, and consumers with their own cross-task state (dedup, rate limiting, batched API calls) typically would not use a per-task-future helper. It can be extracted later when a clear use case appears.

- **Replacing `TaskHandler.isReady()`.** The existing method stays and is still consulted; gates are AND-combined with it. No flag day for existing subclasses.

## Considered Options

- **Option A — `TaskReadinessGate` SPI in core (this proposal).**
- **Option B — Channel operator `glacierRestore()`.**
- **Option C — `FileSystemProvider` interception in nf-amazon's S3 NIO.**

## Pros and Cons of the Options

### Option A — `TaskReadinessGate` SPI in core

A new extension point consulted by `TaskPollingMonitor` before admitting a task for submission. Plugins implementing the interface are discovered via PF4J.

- Good, because it matches the actual problem ("defer task submission until external state is ready"), which is generic across executors and backends.
- Good, because users get drop-in behavior — `plugins { id '...' }` is enough; no `process.executor` override.
- Good, because the plugin stops being tied to a specific executor implementation and works with every executor in a single class.
- Good, because the change to core is small (one interface plus a few lines in `TaskPollingMonitor`) and behavior is bit-identical when no plugin registers a gate.
- Bad, because it requires a coordinated upstream change in Nextflow before the consumer plugin can be refactored.

### Option B — Channel operator `glacierRestore()`

A plugin-provided operator that holds items in a channel until each is verified restored, then emits them: `Channel.fromPath('s3://...') | glacierRestore | MY_PROCESS`.

- Good, because it ships independently of Nextflow with no core changes.
- Good, because it makes the restore step explicit and visible in the pipeline code.
- Bad, because it requires pipeline authors to rewrite channel plumbing.
- Bad, because it does not help users who declare `path('s3://...')` directly as a process input without an explicit upstream channel.
- Best treated as a complement to Option A, not a replacement.

### Option C — `FileSystemProvider` interception in nf-amazon's S3 NIO

Intercept S3 reads at the NIO layer and trigger restore on demand.

- Bad, because Nextflow never reads input bytes on the head node when using AWS Batch — the Batch worker does, via `aws s3 cp` or Fusion. A filesystem hook only catches head-node access (local executor, head-node `publishDir`, etc.) and cannot gate Batch submissions, which is the core requirement.

Rejected. Documented here to avoid revisiting.

## Decision

Adopt Option A. Introduce a `TaskReadinessGate` interface in `nextflow.processor` and consult registered implementations from `TaskPollingMonitor`. Option B may be pursued separately by individual plugins as a complementary surface.

## Core capabilities

### Interface

```groovy
package nextflow.processor

import nextflow.plugin.extension.ExtensionPoint

/**
 * Plugin extension point that defers task submission until external preconditions are met.
 *
 * <p>Invoked by the task scheduler:
 * <ul>
 *   <li>Once when a task is added to the pending queue (return value ignored). This lets
 *       the gate kick off any async preparation, e.g. issuing an S3 RestoreObject request.</li>
 *   <li>On every subsequent polling tick (~100ms) until it returns {@code true}.</li>
 * </ul>
 *
 * <p><b>Implementations must return promptly.</b> This method runs on the task scheduler
 * thread; blocking it stalls submission for every task in the run, not only the one being
 * gated. Kick off async work on the first call and return {@code false}; report status on
 * subsequent calls.
 *
 * <p>Throwing from this method marks the task as permanently failed, with the thrown
 * exception attached as the cause. The standard {@code errorStrategy} applies, so a
 * transient failure can be retried. Returning {@code false} indefinitely keeps the task
 * waiting forever — throw to signal definitive failure.
 *
 * <p>When multiple gates are registered, a task is admitted only when every gate returns
 * {@code true}. Evaluation order across gates is unspecified.
 */
interface TaskReadinessGate extends ExtensionPoint {
    boolean isReady(TaskRun task)
}
```

### Contract

| Aspect | Contract |
|---|---|
| Return value | `true` = ready to submit; `false` = not yet, poll again |
| Permanent failure | Throw an exception; routed through the standard task failure path |
| Latency | Must return promptly. No blocking I/O on the calling thread. |
| First-call semantics | Idempotent. Kick off async preparation, return `false` immediately. |
| Aggregation across gates | Logical AND. Order unspecified. |
| Per-task identity | `TaskRun.id` is stable for dedup and caching inside the gate |

### Call site integration

Three small edits to `nextflow.processor.TaskPollingMonitor`.

**Cache the gate list once per session** in `start()`:

```groovy
private List<TaskReadinessGate> readinessGates = Collections.emptyList()

@Override
TaskMonitor start() {
    readinessGates = Plugins.getExtensions(TaskReadinessGate)
    // ... existing start logic
    return this
}
```

**Fire-once on enqueue** so the gate can start async work the moment the task is queued — not when an executor slot frees up:

```groovy
void schedule(TaskHandler handler) {
    // ... existing enqueue logic
    if( readinessGates ) {
        for( gate in readinessGates ) {
            try {
                gate.isReady(handler.task)
            }
            catch( Throwable t ) {
                handleGateFailure(handler, t)
                return
            }
        }
    }
}
```

The return value is intentionally discarded — this call exists for the side effect.

**AND-in gates with the existing admission checks** in `canSubmit`:

```groovy
protected boolean canSubmit(TaskHandler handler) {
    // The single-threaded scheduler invariant is load-bearing here: TaskReadinessGate
    // implementations are contractually required to return promptly.
    (capacity > 0 ? checkQueueCapacity(handler) : true) \
        && handler.canForkProcess() \
        && allGatesReady(handler) \
        && handler.isReady()
}

private boolean allGatesReady(TaskHandler handler) {
    if( !readinessGates ) return true
    for( gate in readinessGates ) {
        try {
            if( !gate.isReady(handler.task) ) return false
        }
        catch( Throwable t ) {
            handleGateFailure(handler, t)
            return false
        }
    }
    return true
}
```

`handleGateFailure(TaskHandler, Throwable)` routes the exception through the same path used today for submit-time `ProcessException`s: the task transitions to `FAILED` with the gate exception as cause, and the configured `errorStrategy` (`terminate` / `ignore` / `retry`) applies.

### Discovery and lifecycle

- Gates are resolved via `Plugins.getExtensions(TaskReadinessGate)` — the standard PF4J path used by `TraceObserverFactory`, `ContainerResolver`, and other plugin extension points.
- Instances are constructed once when the monitor starts and reused for the run.
- No teardown hook on the SPI. Plugins that need cleanup own that via their existing `BasePlugin.start/stop` lifecycle.
- No config flag. A plugin that ships a gate is automatically consulted for every task; opting out means not installing the plugin.

### Backward compatibility

- With no gates registered, `canSubmit` evaluates one empty-list check per pending task per tick. `schedule` adds one untaken branch. Behavior is bit-identical to today.
- `TaskHandler.isReady()` is unchanged and still consulted. Existing executor and task-handler subclasses keep working.
- The SPI is purely additive: no existing signatures change.

## Rationale & discussion

**Why polling rather than blocking.** `TaskPollingMonitor` runs a single scheduler thread that walks the pending queue every ~100ms. A blocking gate (one that waits inside `isReady` for the external precondition to be satisfied) would stall submission for every task in the run, not just the one being gated. For the driving use case — Glacier restores that can take up to 12 hours — this is unacceptable. The contract is therefore explicitly polling-based: the first call kicks off async work and returns `false`; subsequent calls cheaply check status. This also matches the AWS RestoreObject API, which is async-by-design (issue `RestoreObject`, then poll `HeadObject` for status).

**Why fire on enqueue, not only on `canSubmit`.** If the gate were only consulted from `canSubmit`, a task waiting behind a full executor queue or a `maxForks` limit would not trigger the gate, so the restore would not begin until capacity frees up. Adding the fire-once call in `schedule()` decouples "start the restore" from "is there a slot." The cost is one extra call per task per gate per run; the call is contractually idempotent so this is safe.

**Why exceptions for permanent failure.** Returning `false` indefinitely would keep a task waiting forever in cases where success is impossible (object permanently inaccessible, restore window timed out, AccessDenied). Adding a tri-state return value (`READY`/`WAITING`/`FAILED`) would enlarge the public API for the same effect Java already provides via exceptions. Routing through the existing failure path also means `errorStrategy 'retry'` naturally handles transient failures — a `ThrottlingException` bubbling out of the gate triggers the same retry behavior as a transient executor failure.

**Why no helper class.** An `AbstractAsyncReadinessGate` base class that runs blocking work on a dedicated executor and exposes results via per-task `CompletableFuture`s would make safe-by-default gates trivial to write. The reason we are skipping it in v1: realistic consumers (e.g. a Glacier restorer with cross-task dedup, rate limiting, and prefix expansion) tend to have their own async state and would not use a per-task-future helper. Shipping the helper now means designing it for hypothetical consumers whose shape we cannot yet see. Extracting it later is non-breaking.

**Why this lives in `nextflow.processor`, not `nf-commons`.** Other task-scoped extension points (`TaskTipProvider`, the task handler itself) live in `nextflow.processor`. `nf-commons` hosts cross-cutting infrastructure (the `ExtensionPoint` marker itself, plugin loading machinery). Co-locating `TaskReadinessGate` with `TaskRun` and `TaskHandler` matches the existing layering.

## Testing

New Spock specs in `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorTest.groovy`:

- No gates registered → `canSubmit` matches current behavior; `schedule` does not touch gate code.
- One gate returning `true` → task admitted as soon as capacity and forks allow.
- One gate returning `false` then `true` → task waits in pending across ticks; admitted on the tick the gate flips.
- Gate invoked exactly once on enqueue; return value ignored even when `false`.
- Multiple gates → AND semantics; one `false` blocks admission.
- Gate throws on enqueue → task fails with the thrown cause; not added to running queue.
- Gate throws during poll → task fails with the thrown cause; removed from pending.
- Gate throws under `errorStrategy 'retry'` → standard retry path.

Gates are injected through a test double rather than a real PF4J context.

