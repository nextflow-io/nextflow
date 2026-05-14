# `service` process directive for long-running auxiliary tasks

- Authors: Edmund Miller
- Status: draft
- Deciders: TBD
- Date: 2026-05-14
- Tags: directive, process, lifecycle, executor

## Summary

Introduce a `service` process directive that marks a process as a long-running auxiliary task whose declared outputs are emitted to downstream channels as soon as they appear in the task work directory, and which is automatically terminated (SIGTERM) once all consumer processes have finished.

## Problem Statement

Some pipeline steps depend on a co-running auxiliary process rather than on a file that is produced once and then consumed. Examples include:

- A local database server (`postgres`, `redis`, `chromadb`) that several downstream tools query.
- A long-lived inference server (e.g. `vllm`, `triton`, `ollama`) that batches multiple consumer tasks against a single GPU model load.
- A named pipe, ramdisk-backed shared workspace, or socket-based message bus.
- A local HTTP service exposing test fixtures (`localstack`, `minio`, mock APIs) to integration-style consumer tasks.

In Nextflow today there is no first-class way to express this dependency. Process outputs are bound to downstream channels *only* when the task script exits, so a script that creates a socket then blocks (the usual shape of a server) never unblocks its consumers. Users currently work around this with `nohup` + manual `process.start.sh` files, with `beforeScript` hooks that leak processes across resume, or by serializing the whole pipeline into a single process — all brittle and non-portable.

Snakemake offers exactly this pattern via "service rules" (Snakemake 7+). Bringing the same semantic to Nextflow gives users a clean, portable way to model the co-running-resource case.

## Goals or Decision Drivers

- Allow a process to expose a resource (socket, port, file, ramdisk) while still running, so consumer processes can use it before the producer exits.
- Tie the producer's lifetime to its consumers — terminate it automatically when no consumer needs it anymore.
- Reuse existing Nextflow plumbing (directives, `TaskHandler.kill()`, `TraceObserverV2`) rather than introducing parallel lifecycles.
- Be opt-in and have zero impact on processes that don't use it.
- Keep the MVP small and reviewable; expand executor coverage in follow-ups.

## Non-goals

- Cluster-wide or cross-node co-location of producer and consumers. The MVP runs on the **local executor only** and assumes producer and consumers share a POSIX work directory.
- Group-local service instances (one service per job group, as in Snakemake's `foo.{groupid}.socket` pattern). Single-instance services only for the MVP; multi-instance variants can be layered on later if there is demand.
- Networked service discovery. The user is responsible for choosing a stable port / path and communicating it to consumers via the declared output.
- Per-output `service(...)` markers in the Snakemake style. Nextflow's directive model is process-level, and the simpler `service true` directive captures the same intent.

## Considered Options

- **Option 1**: New process-level directive `service true`, output emission triggered when declared output paths appear, termination via observer + `TaskHandler.kill()`.
- **Option 2**: Per-output marker, e.g. `output: path 'foo.sock', service: true`, mirroring Snakemake's `service("foo.socket")`.
- **Option 3**: No new primitive — document a recipe using `beforeScript`, `nohup`, and a sentinel file. Status quo.

## Pros and Cons of the Options

### Option 1 — process-level `service true` (chosen)

- Good, because it fits Nextflow's existing directive vocabulary (`cpus`, `memory`, `debug`, etc.).
- Good, because it reuses `TaskHandler.kill()`, `Session` event notifications, and the v2 trace-observer infrastructure — no new lifecycle states.
- Good, because declarative termination via the DAG matches user intuition ("kill when no one needs this anymore").
- Bad, because every output of a service process is implicitly treated as a service resource; if a service legitimately produces both a long-lived socket and a one-shot result file, the user has to model it as two processes.
- Bad, because cache/resume semantics are forced off — service tasks must always re-run.

### Option 2 — per-output `service(...)` marker

- Good, because it mirrors Snakemake exactly, easing migration.
- Good, because a single process can mix service and non-service outputs.
- Bad, because it requires changes to `OutParam` parsing in both DSL v1 and v2, a much larger surface area.
- Bad, because the semantic question "is this *task* a service?" gets split across multiple outputs, complicating the lifecycle observer.

### Option 3 — recipe-only, no language change

- Good, because no engine changes are needed.
- Bad, because the workarounds (background `nohup`, sentinel files, manual `kill -TERM`) are fragile, leak on abort, and break under `nextflow resume`.
- Bad, because it pushes complex lifecycle reasoning onto every user that hits the case.

## Solution or decision outcome

Adopt **Option 1**: a process-level `service true` directive. Service outputs are emitted as soon as the declared output paths appear in the work directory. A `ServiceLifecycleObserver` (a `TraceObserverV2`) tracks consumer processes via the workflow DAG and calls `TaskHandler.kill()` on the service task once all consumers have terminated. Initial executor coverage is **local only**; any non-local executor combined with `service true` raises a clear validation error.

## Rationale & discussion

### Output binding while the task is still running

Today `TaskProcessor.bindOutputs()` is invoked only on task completion. For a service this is too late — the consumer would wait for a script that never exits. We add a background watcher (`TaskProcessor.startServiceOutputWatcher`) kicked off when the task transitions to `RUNNING`. It polls the work directory using the existing `collectOutputs(task)` machinery, catches the "outputs not yet present" exceptions (`MissingFileException`, `IllegalArityException`), and once `collectOutputs` succeeds calls a private `bindServiceOutputs(task)` that:

1. sets a new `TaskRun.serviceOutputsBound` flag (so any later completion-time `bindOutputs` call is a no-op),
2. delegates to the existing `bindOutputs0` to emit values into the downstream channels.

Polling is a deliberate MVP choice over `WatchService`/inotify: it avoids initial-scan races (a fast `touch` before the watcher registers would otherwise be missed), works uniformly across local filesystems, and is bounded by `NXF_SERVICE_READY_POLL_MS` / `NXF_SERVICE_READY_TIMEOUT_MS` for tuning.

### Terminating the service

`TaskHandler.kill()` already exists and, on the local executor, sends SIGTERM via `LocalTaskHandler.killTask()`. We reuse it. The new `ServiceLifecycleObserver`:

- At `onFlowBegin`, walks `session.dag` and records each service process's downstream consumer processes (traversing through OPERATOR vertices, stopping at the first PROCESS vertex).
- At `onTaskStart`, captures the running service handler.
- At `onProcessTerminate(consumer)`, removes the consumer from each service's pending set; when a service's pending set empties, calls `handler.kill()`.
- At `onTaskComplete` for a service that exited unexpectedly (before we killed it), aborts the session with a clear error.
- At `onFlowComplete` / `onFlowError`, kills any still-tracked service as a safety net.

Choosing `onProcessTerminate` over per-task counting avoids fragile race conditions when a consumer process spawns more tasks dynamically: the event fires once the whole consumer process is finished, which is the correct trigger.

### Why force `cache false`

A service task's "output" is a transient resource — a Unix domain socket, a port-bound listener, a ramdisk symlink. Reusing yesterday's cached entry would point consumers at a dead socket. Making `isCacheable()` return `false` when `service == true` removes the foot-gun.

### Why local-executor-only for the MVP

A service must share a filesystem with its consumers and respond to SIGTERM from the same host as the Nextflow head process. The local executor satisfies both trivially. Grid executors (Slurm/LSF/SGE/PBS) often share a POSIX filesystem but require additional co-location enforcement so the consumer's compute node can reach the producer's socket; that is best handled in a follow-up that also defines a `service.colocate` strategy. Cloud and Kubernetes executors fundamentally violate the shared-FS assumption and require a different design (sidecar containers, headless services). Failing fast on these executors keeps the MVP honest.

### What gets validated up front

- `service true` + non-local executor → `IllegalArgumentException` at `TaskProcessor.run()` startup, before any task is submitted.
- Service processes are treated as singletons. Pluralized inputs that would otherwise produce multiple tasks per service process will be rejected in a follow-up validation pass (out of scope here; documented as a known limitation).

### Configuration surface

Two environment variables tune the watcher and are intentionally undocumented in the user-facing reference until we see real-world need:

- `NXF_SERVICE_READY_POLL_MS` (default `250`) — poll interval while waiting for declared outputs to appear.
- `NXF_SERVICE_READY_TIMEOUT_MS` (default `600000`) — abort if outputs don't appear within this window.

### Example

```groovy
process the_service {
    service true

    output:
    path 'foo.sock'

    script:
    """
    ln -s /dev/random foo.sock
    sleep 10000
    """
}

process consumer {
    input:
    path sock

    output:
    path 'out.txt'

    script:
    """
    head -n1 ${sock} > out.txt
    """
}

workflow {
    sock = the_service()
    consumer(sock)
    consumer2(sock)
}
```

`the_service` keeps running until both `consumer` and `consumer2` have finished, at which point Nextflow SIGTERMs it. The exit code (143) is treated as success for a service process.

## Links

- Snakemake service rules: <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#service-rules>
- Related directive: [hints](20260323-hints-process-directive.md) — same pattern of small, opt-in process directives
