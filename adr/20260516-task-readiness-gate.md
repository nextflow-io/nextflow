# `TaskReadinessGate` plugin extension point for deferred task submission

- Authors: Rob Syme
- Status: accepted
- Deciders: Paolo Di Tommaso, Ben Sherman, Rob Syme
- Date: 2026-05-20
- Tags: plugin, extension-point, scheduler, task-submission

## Summary

Introduce a `TaskReadinessGate` plugin extension point that lets a plugin defer task submission until an external precondition is met (e.g. an S3 object has been restored from Glacier), removing the need for plugins to subclass an executor and its task handler. Plugins implement a blocking `void prepare(TaskHandler)` method; core runs each gate on a managed virtual-thread executor via a new `TaskGateManager`, and the standard `TaskPollingMonitor` delegates to the manager from `schedule`, `canSubmit`, and `evict`.

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

- Keep the plugin-author surface minimal — gate authors write straightforward blocking code, not an async state machine.

- Allow gates to kick off preparation work the moment a task is queued, not when an executor slot frees up — cold-storage restores can take hours and must not wait for capacity.

- Allow plugins to signal permanent task failure cleanly, integrating with the standard `errorStrategy` machinery (including retry routing through `ProcessRetryableException`).

## Non-goals

- **Per-process gate scoping as a new directive.** A gate is consulted for every task. Plugins implement per-process opt-out by reading the existing `hints` directive (e.g. `hints 'glacier/skip': true`).

- **An executor-level timeout.** Core does not enforce a wall-clock deadline on `prepare`. The right value depends entirely on which plugin is registered and on the workload (Glacier Standard is hours; Deep Archive Bulk is days). Plugins own timeout policy, typically reading a per-process hint (e.g. `hints 'glacier/maxWait': '5h'`) and throwing when the deadline elapses.

- **Gate ordering or priority.** Evaluation order across gates is unspecified; all gates run in parallel.

- **An `AbstractAsyncReadinessGate` helper class.** Since the SPI is blocking and core owns the executor, there is no async runtime for a helper to wrap. Plugin authors write blocking code directly.

- **Replacing `TaskHandler.isReady()`.** The existing method stays and is still consulted; gates are AND-combined with it. No flag day for existing subclasses.

## Considered Options

- **Option A — `TaskReadinessGate` SPI in core (this proposal).**
- **Option B — Channel operator `glacierRestore()`.**
- **Option C — `FileSystemProvider` interception in nf-amazon's S3 NIO.**

Within Option A, there was a sub-decision between two contract shapes:

- **A.1 polling** — `boolean isReady(TaskRun)`, called repeatedly on the scheduler thread; plugin manages its own async state.
- **A.2 blocking** — `void prepare(TaskHandler) throws InterruptedException`, called once per task on a core-managed virtual-thread executor; plugin writes straightforward blocking code.

The blocking variant (A.2) was chosen — see "Why blocking, not polling" below.

## Pros and Cons of the Options

### Option A — `TaskReadinessGate` SPI in core

A new extension point consulted by `TaskPollingMonitor` before admitting a task for submission. Plugins implementing the interface are discovered via PF4J.

- Good, because it matches the actual problem ("defer task submission until external state is ready"), which is generic across executors and backends.
- Good, because users get drop-in behavior — `plugins { id '...' }` is enough; no `process.executor` override.
- Good, because the plugin stops being tied to a specific executor implementation and works with every executor in a single class.
- Good, because the change to core is small (one interface, one manager class, ~14 lines in `TaskPollingMonitor`) and behavior is bit-identical when no plugin registers a gate.
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

Adopt Option A with the blocking-`prepare` contract (A.2). Introduce a `TaskReadinessGate` interface in `nextflow.processor` and a `TaskGateManager` that orchestrates registered implementations on behalf of `TaskPollingMonitor`. Option B may be pursued separately by individual plugins as a complementary surface.

## Core capabilities

### Interface

```groovy
package nextflow.processor

import org.pf4j.ExtensionPoint

/**
 * Plugin extension point that defers task submission until an external precondition is met.
 *
 * <p>Implementations are invoked once per task by the scheduler, on a managed background
 * thread (virtual thread when available). The task is admitted for submission when this
 * method returns; throw to mark the task as permanently failed and route the cause through
 * the task's {@code errorStrategy}.
 *
 * <p>Implementations may block freely — {@code Thread.sleep}, network calls, long polling.
 * The scheduler thread is never blocked by this call.
 *
 * <p>Implementations must honor {@code Thread.interrupt()} so that task eviction and
 * workflow abort can unblock {@code prepare} promptly. Core does not enforce a wall-clock
 * deadline — plugins own timeout policy (see the developer guide for the recommended
 * {@code hints}-based pattern).
 *
 * <p>When multiple gates are registered, a task is admitted only when every gate's
 * {@code prepare} method has returned successfully. Evaluation order across gates is
 * unspecified; all gates start in parallel on the managed executor.
 */
interface TaskReadinessGate extends ExtensionPoint {
    void prepare(TaskHandler handler) throws InterruptedException
}
```

### Contract

| Aspect | Contract |
|---|---|
| Return | Method returns normally → task is ready to submit |
| Permanent failure | Throw any `Throwable`; `ProcessException` is canonical |
| Cancellation | `InterruptedException` propagates as cancellation; gate must honor `Thread.interrupt()` |
| Threading | Runs on a core-managed virtual-thread executor; blocking I/O is fine |
| Aggregation across gates | All gates must complete successfully; one failure aborts the task |
| Order across gates | Unspecified; gates run in parallel |
| Per-task identity | `handler.task.id` stable for plugin-internal dedup/caching |
| Per-process opt-out / timeout | Plugin reads `handler.task.config.hints` and enforces its own policy |

### Architecture

The orchestration lives in a new `TaskGateManager` class that owns:

- The registered gates list (resolved once via `Plugins.getExtensions(TaskReadinessGate)`).
- A managed `ExecutorService` (virtual-thread per task when `Threads.useVirtual()`, cached pool otherwise) created lazily only when at least one gate is registered.
- A `ConcurrentMap<TaskHandler, List<Future<?>>>` tracking in-flight gate work per handler.
- The session shutdown hook that drains the executor.

The manager exposes three methods to `TaskPollingMonitor`:

```groovy
void submit(TaskHandler handler)   // submit prepare() futures for each gate
boolean isReady(TaskHandler handler)   // poll futures; throw on failure; true when all done
void evict(TaskHandler handler)        // cancel any in-flight futures
```

### Monitor integration

`TaskPollingMonitor` adds one field and four delegation lines. The total diff against upstream is ~14 lines.

```groovy
@PackageScope
TaskGateManager gateManager = new TaskGateManager(null, [])   // empty until start()

@Override
TaskMonitor start() {
    this.gateManager = new TaskGateManager(session)
    // ... existing start logic
}

void schedule(TaskHandler handler) {
    pendingLock.lock()
    try {
        pendingQueue << handler
        gateManager.submit(handler)
        // ... existing signal/notify
    } finally { pendingLock.unlock() }
}

protected boolean canSubmit(TaskHandler handler) {
    gateManager.isReady(handler)
        && handler.canForkProcess()
        && handler.isReady()
        && (capacity > 0 ? checkQueueCapacity(handler) : true)
}

boolean evict(TaskHandler handler) {
    if( !handler ) return false
    gateManager.evict(handler)
    // ... existing remove/signal
}
```

When `gateManager.isReady` throws, the existing catch in `submitPendingTasks` routes the cause through `handleException` → `resumeOrDie`, which dispatches via the task's `errorStrategy`. No new failure path.

### Exception handling

`TaskGateManager.isReady` unwraps `ExecutionException` from the future and preserves cause types so retry routing works:

- `ProcessException` (and subclasses) thrown by `prepare` propagate identity-preserved.
- Any other throwable is wrapped in `new ProcessException("Task readiness gate failed...", cause)`, so markers like `ProcessRetryableException` carried on a `RuntimeException` reach `TaskProcessor.resumeOrDie` via `error.cause`.
- When one gate throws, peer futures for the same handler are cancelled so background work does not outlive the failing task.

### Discovery and lifecycle

- Gates are resolved via `Plugins.getExtensions(TaskReadinessGate)` — the standard PF4J path used by `TraceObserverFactory`, `ContainerResolver`, and other plugin extension points.
- Instances are constructed once when the monitor starts and reused for the run.
- No teardown hook on the SPI. Plugins that need cleanup own that via their existing `BasePlugin.start/stop` lifecycle.
- The managed executor is shut down via `session.onShutdown` (`shutdownNow()` interrupts any in-flight gates).
- No config flag. A plugin that ships a gate is automatically consulted for every task; opting out means not installing the plugin (or using `hints` for per-process opt-out).

### Per-process opt-out and timeout via `hints`

Plugins read process-level overrides from the existing `hints` directive — no new core directive required. Namespaced keys (`<plugin>/<name>`) avoid collisions across plugins.

```nextflow
process FAST_PATH {
    hints 'glacier/skip': true
    // ...
}

process LONG_RESTORE {
    hints 'glacier/maxWait': '48h'
    // ...
}
```

Inside the gate:

```groovy
@Override
void prepare(TaskHandler handler) throws InterruptedException {
    final hints = handler.task.config.hints
    if( hints['glacier/skip'] == true ) return

    final maxWait = (hints['glacier/maxWait'] as Duration) ?: defaultMaxWait
    final deadline = System.currentTimeMillis() + maxWait.toMillis()

    // ... enforce deadline + interrupts
}
```

### Backward compatibility

- With no gates registered, `gateManager.submit` and `gateManager.evict` are no-ops and `gateManager.isReady` returns `true` immediately. No executor is created. Behavior is bit-identical to upstream.
- `TaskHandler.isReady()` is unchanged and still consulted. Existing executor and task-handler subclasses keep working.
- The SPI is purely additive: no existing signatures change.

## Rationale & discussion

### Why blocking, not polling

An earlier iteration of this ADR proposed a polling contract — `boolean isReady(TaskRun)` called repeatedly on the scheduler thread, with the plugin responsible for kicking off async work idempotently on the first call. Review pushed back on three grounds:

1. **Plugin-author cognitive load.** "Must return promptly, kick off async on first call, idempotent thereafter, manage your own state" is a lot of contract to get right. A blocking method with the standard `Thread.interrupt` cancellation contract is what plugin authors already know how to write.

2. **`FilePorter` precedent.** The in-tree `FilePorter` already uses a managed-pool + blocking-`Runnable` shape for similar "prepare external resource for a task" work. New extension points should follow.

3. **Virtual threads make the blocking cost negligible.** A gate parked on a virtual thread costs about 1 KB. The "blocking gate would stall the scheduler" objection that originally motivated polling no longer applies — the scheduler thread doesn't run the gate; the managed virtual-thread executor does.

The blocking design moves the async management cost from every plugin author into one place in core (`TaskGateManager`), and the plugin-facing API shrinks to a single method with no idempotency reasoning.

### Why no helper class

The earlier proposal noted an `AbstractAsyncReadinessGate` helper as a likely follow-up. With the blocking design, no helper is needed: plugin authors write blocking code directly, and core owns the executor. The `AbstractAsyncReadinessGate` non-goal is now permanent rather than deferred.

### Why no executor-level timeout

An earlier iteration carried an `executor.gateMaxWait` config option as a safety net for stuck gates. It was removed during review:

- The right deadline depends on which plugin is registered (Glacier Standard hours vs Deep Archive Bulk days). A one-size-fits-all executor-level default either fires on every legitimate restore or is functionally "no limit" for buggy gates.
- The safety net was always partial: a plugin that ignores `Thread.interrupt()` cannot be unblocked by a wall-clock timeout. Such plugins already break eviction, so the safety net never delivered on its premise.
- Discoverability is better when the knob lives with the plugin that defines what a sensible value is. A user reading `executor.gateMaxWait` in the config reference cannot tell what value is right; a user reading the plugin's docs for `hints 'glacier/maxWait'` can.

Plugins enforce their own deadlines. The hints-based pattern documented above gives natural per-process scoping for free.

### Why fire on `schedule`, not only on `canSubmit`

If gate work only started when `canSubmit` was first reached, a task waiting behind a full executor queue or a `maxForks` limit would not trigger preparation until capacity frees up. The fire-on-`schedule` design decouples "start the restore" from "is there a slot," so a 5-hour Glacier restore overlaps with whatever else is running rather than waiting in line for it.

### Why exceptions for permanent failure

Routing through the existing failure path means `errorStrategy 'retry'` naturally handles transient failures — a `ThrottlingException` bubbling out of the gate triggers the same retry behavior as a transient executor failure. Returning a tri-state value (`READY`/`WAITING`/`FAILED`) would enlarge the public API for the same effect Java already provides via exceptions.

Cause-type preservation is load-bearing: `TaskProcessor.resumeOrDie` checks `error.cause instanceof ProcessRetryableException` (and `CloudSpotTerminationException`) — it inspects the cause, not the thrown exception itself. The manager therefore wraps non-`ProcessException` causes in a `ProcessException` so the marker reaches `resumeOrDie` via `error.cause`.

### Why this lives in `nextflow.processor`

Other task-scoped extension points (`TaskTipProvider`, the task handler itself) live in `nextflow.processor`. `nf-commons` hosts cross-cutting infrastructure (plugin loading machinery). Co-locating `TaskReadinessGate` and `TaskGateManager` with `TaskRun` and `TaskHandler` matches the existing layering.

## Testing

- `TaskReadinessGateTest` — interface shape: extends `org.pf4j.ExtensionPoint`, declares `prepare(TaskHandler) throws InterruptedException`.
- `TaskGateManagerTest` — full behavioural coverage targeting the manager directly: no-gates no-op, single-gate happy path, multi-gate all-must-complete, fail-fast peer cancellation, eviction interrupt propagation, `ProcessException` identity preservation, `ProcessRetryableException` routing via cause, `RuntimeException` wrapping, checked-exception wrapping, external cancellation.
- The standard `TaskPollingMonitorTest` and `ParallelPollingMonitorTest` continue to pass with a one-line stub adjustment for the `canSubmit` AND-chain reorder.

## References

- Implementation PR: [nextflow-io/nextflow#7158](https://github.com/nextflow-io/nextflow/pull/7158).
- This ADR PR: [nextflow-io/nextflow#7151](https://github.com/nextflow-io/nextflow/pull/7151).
