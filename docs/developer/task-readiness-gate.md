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
- **Retry markers** (`ProcessRetryableException`, `CloudSpotTerminationException`) are recognised by `resumeOrDie` on the *cause* of the thrown exception, not on the exception itself. If you want `errorStrategy 'retry'` to fire for a transient failure, throw the marker as-is (it will be wrapped) — do not pre-wrap it in a `ProcessException`, since the `ProcessException` would propagate identity-preserved and the marker would be lost.
- **Interrupts must be honored.** Task eviction and workflow abort cancel the in-flight `prepare` by interrupting its thread. Use interruptible primitives (`Thread.sleep`, blocking I/O on NIO channels, `Future.get`). Core does not enforce a wall-clock deadline — the plugin owns timeout policy.
- **Multiple gates compose.** When several plugins register gates, all must complete successfully before the task is admitted. Gates run in parallel; order is unspecified.

## Timeouts and per-process overrides via `hints`

Core deliberately does not ship an executor-level timeout. The right deadline depends on which plugin is registered and on the workload (e.g. Glacier Standard is hours; Deep Archive Bulk is days). Plugins enforce their own deadlines via the existing `hints` directive — namespaced under the plugin's name to avoid collisions.

```nextflow
process FAST_PATH {
    hints 'glacier/skip': true        // bypass the gate entirely
    // ...
}

process LONG_RESTORE {
    hints 'glacier/maxWait': '48h'    // per-process deadline override
    // ...
}
```

Inside the gate, read `handler.task.config.hints` and enforce the deadline yourself:

```groovy
@Override
void prepare(TaskHandler handler) throws InterruptedException {
    final hints = handler.task.config.hints
    if( hints['glacier/skip'] == true ) return

    final maxWait = (hints['glacier/maxWait'] as Duration) ?: defaultMaxWait
    final deadline = System.currentTimeMillis() + maxWait.toMillis()

    for( S3Path path : extractS3Inputs(handler.task) ) {
        manager.issueRestoreIfNeeded(path)
        while( !manager.isRestored(path) ) {
            if( System.currentTimeMillis() > deadline )
                throw new ProcessException("Glacier restore exceeded ${maxWait} for ${path}")
            if( manager.isPermanentlyFailed(path) )
                throw new ProcessException("Glacier restore failed for ${path}")
            Thread.sleep(60_000)
        }
    }
}
```

The `hints` directive uses prefix-separated keys (`<plugin>/<name>`); plugins should namespace their hint keys (e.g. `glacier/skip`, `mycorp.cold-storage/skip`) to avoid collisions with other plugins.

## When `prepare` runs

`TaskPollingMonitor.schedule()` submits `prepare` to the managed executor the moment the task is enqueued — before any executor slot has freed up, before `canForkProcess()` is consulted. This means restore work for tasks queued behind a full executor begins immediately, not when slots free up.

`canSubmit()` then polls the resulting `Future` on every monitor tick. The task is admitted as soon as every gate's future has completed successfully and the standard `canForkProcess` / `isReady` / capacity checks pass.
