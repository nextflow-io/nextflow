(task-readiness-gate)=

# `TaskReadinessGate`

`TaskReadinessGate` is a plugin extension point that defers task submission until an external precondition is met — for example, restoring an S3 object from Glacier before an AWS Batch worker tries to stage it. It works uniformly with every executor that uses Nextflow's `TaskPollingMonitor` and removes the need for plugins to subclass an executor and its task handler purely to override `TaskHandler.isReady()`.

## Interface

A gate implements one method:

```groovy
package nextflow.processor

import org.pf4j.ExtensionPoint

interface TaskReadinessGate extends ExtensionPoint {
    void prepare(TaskHandler handler) throws InterruptedException
}
```

The plugin registers an implementation via the standard PF4J `@Extension` mechanism, the same way `TraceObserverFactory`, `CacheFactory`, and other extension points are discovered.

## Contract

- **Blocking is allowed.** `prepare` runs on a managed virtual-thread executor inside `TaskPollingMonitor`. Calling `Thread.sleep`, blocking I/O, or long-polling APIs is fine. The scheduler thread is never blocked.
- **Throwing fails the task.** Any exception marks the task as permanently failed and routes the cause through the task's `errorStrategy` directive. `ProcessException` (and subclasses) propagate identity-preserved. Other throwables are wrapped in a `ProcessException` with the original attached as `cause`, so retry markers like `ProcessRetryableException` reach `TaskProcessor.resumeOrDie` intact.
- **Interrupts must be honored.** Task eviction, workflow abort, and the `executor.gateMaxWait` safety net cancel the in-flight `prepare` by interrupting its thread. Use interruptible primitives (`Thread.sleep`, blocking I/O on NIO channels, `Future.get`).
- **Multiple gates compose.** When several plugins register gates, all must complete successfully before the task is admitted. Gates run in parallel; order is unspecified.

## Configuration

`executor.gateMaxWait` bounds the time a gate may spend preparing a single task (default `24h`). Tasks whose gates exceed this limit are cancelled and fail with a `ProcessException` that the task's `errorStrategy` can handle (including `retry`).

## Per-process opt-out

Use the existing `hints` directive — no core change required:

```nextflow
process FAST_PATH {
    hints 'glacier/skip': true
    // ...
}
```

Inside the gate, inspect `handler.task.config.hints` and short-circuit:

```groovy
@Override
void prepare(TaskHandler handler) throws InterruptedException {
    if( handler.task.config.hints['glacier/skip'] == true ) return
    // ... real work
}
```

The `hints` directive uses dot-separated and prefix-separated keys; plugins should namespace their hint keys (e.g. `glacier/skip`, `mycorp.cold-storage/skip`) to avoid collisions.

## Example

A minimal gate that issues an S3 Glacier restore request and waits for completion:

```groovy
@CompileStatic
class GlacierReadinessGate implements TaskReadinessGate {

    private final GlacierRestoreManager manager

    GlacierReadinessGate() {
        this.manager = new GlacierRestoreManager(/* config from session */)
    }

    @Override
    void prepare(TaskHandler handler) throws InterruptedException {
        if( handler.task.config.hints['glacier/skip'] == true ) return
        for( S3Path path : extractS3Inputs(handler.task) ) {
            manager.issueRestoreIfNeeded(path)        // idempotent
            while( !manager.isRestored(path) ) {
                if( manager.isPermanentlyFailed(path) )
                    throw new ProcessException("Glacier restore failed for ${path}")
                Thread.sleep(60_000)                  // virtual thread parks
            }
        }
    }
}
```

## When `prepare` runs

`TaskPollingMonitor.schedule()` submits `prepare` to the managed executor the moment the task is enqueued — before any executor slot has freed up, before `canForkProcess()` is consulted. This means restore work for tasks queued behind a full executor begins immediately, not when slots free up.

`canSubmit()` then polls the resulting `Future` on every monitor tick. The task is admitted as soon as every gate's future has completed successfully and the standard `canForkProcess` / `isReady` / capacity checks pass.
