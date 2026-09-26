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

- **GPU inference microservices**, e.g. an [NVIDIA NIM](https://docs.nvidia.com/nim/) container, [vLLM](https://docs.vllm.ai/), Triton, or Ollama, exposing a local HTTP endpoint. The model weights take tens of seconds to minutes to load and occupy several GB of GPU memory; loading them once and amortizing over many consumer tasks (e.g. one task per sample in a cohort) is dramatically cheaper than re-loading per task.
- **Embedded analytical databases run as a server**, e.g. [DuckDB's `quack` extension](https://duckdb.org/quack/) which turns DuckDB into a client-server system so that multiple consumer processes can run concurrent read/write queries against a shared in-memory or on-disk database.
- **Traditional databases as scratch state for a pipeline**: a `postgres`, `redis`, or `chromadb` instance brought up for the duration of a workflow run, populated by a producer step, queried by several downstream analyses, then torn down.
- **Named pipes, ramdisk-backed shared workspaces, or socket-based message buses** used to stream data between a producer and many consumers without going through the work-directory filesystem.
- **Local HTTP fixtures** for integration-style steps: `localstack`, `minio`, mock APIs, or a small test web server that consumer tasks hit over `http://127.0.0.1:<port>`.

In Nextflow today there is no first-class way to express this dependency. Process outputs are bound to downstream channels *only* when the task script exits, so a script that creates a socket then blocks (the usual shape of a server) never unblocks its consumers. Users currently work around this with `nohup` + manual `process.start.sh` files, with `beforeScript` hooks that leak processes across resume, or by serializing the whole pipeline into a single process — all brittle and non-portable.

Snakemake offers exactly this pattern via "service rules" (Snakemake 7+). Bringing the same semantic to Nextflow gives users a clean, portable way to model the co-running-resource case.

## Goals or Decision Drivers

- Allow a process to expose a resource (socket, port, file, ramdisk) while still running, so consumer processes can use it before the producer exits.
- Tie the producer's lifetime to its consumers — terminate it automatically when no consumer needs it anymore.
- Reuse existing Nextflow plumbing (directives, `TaskHandler.kill()`, `TraceObserverV2`) rather than introducing parallel lifecycles.
- Be opt-in and have zero impact on processes that don't use it.
- Keep the MVP small and reviewable; expand executor coverage in follow-ups.

## Non-goals (for the MVP)

The MVP is deliberately small. The following are explicit non-goals for the first release but are addressed in [Future work](#future-work) below:

- Cluster-wide or cross-node co-location of producer and consumers. The MVP runs on the **local executor only** and assumes producer and consumers share a POSIX work directory.
- Group-local service instances (one service per job group, as in Snakemake's `foo.{groupid}.socket` pattern). Single-instance services only for the MVP; multi-instance variants can be layered on later.
- Networked service discovery. The user is responsible for choosing a stable port / path and communicating it to consumers via the declared output.

Permanent non-goal:

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

### Examples

#### Minimal: socket-as-resource

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

workflow {
    sock = the_service()
    consumer1(sock)
    consumer2(sock)
}
```

`the_service` keeps running until both `consumer1` and `consumer2` have finished, at which point Nextflow SIGTERMs it. The exit code (143) is treated as success for a service process.

#### NVIDIA NIM inference server

A NIM container loads a 7B-parameter LLM into GPU memory once, then serves many consumer tasks (e.g. one variant-annotation call per sample) over HTTP. Without `service`, each sample task would re-load the model — minutes of GPU warm-up multiplied by N samples.

```groovy
process nim_server {
    service true
    accelerator 1, type: 'nvidia-h100'
    container 'nvcr.io/nim/meta/llama-3.1-8b-instruct:latest'

    output:
    path 'endpoint.txt'

    script:
    """
    # NIM listens on :8000 by default; write the endpoint file once
    # the server is ready, then block until SIGTERM'd.
    /opt/nim/start-server.sh &
    NIM_PID=\$!
    until curl -sf http://127.0.0.1:8000/v1/health/ready >/dev/null; do sleep 1; done
    echo "http://127.0.0.1:8000" > endpoint.txt
    wait \$NIM_PID
    """
}

process annotate {
    input:
    path endpoint
    path sample_vcf

    output:
    path "${sample_vcf.baseName}.annotated.json"

    script:
    """
    URL=\$(cat ${endpoint})
    annotate-variants --endpoint \$URL --in ${sample_vcf} \\
        --out ${sample_vcf.baseName}.annotated.json
    """
}

workflow {
    endpoint = nim_server()
    annotate(endpoint, channel.fromPath('samples/*.vcf.gz'))
}
```

The watcher unblocks consumers as soon as `endpoint.txt` is written (which the script writes only after `/v1/health/ready` reports 200), so consumers never race the model load. The single NIM stays up across all sample tasks and is SIGTERM'd once the last `annotate` finishes.

#### DuckDB-as-a-server (quack)

A producer step ingests raw data into a DuckDB database. Several downstream analyses query it concurrently — something a vanilla embedded DuckDB can't do because it requires single-writer/multi-reader file locking. The [`quack` extension](https://duckdb.org/quack/) exposes the database over a network socket, which is exactly the shape `service true` is for.

```groovy
process duckdb_server {
    service true

    input:
    path 'raw/*.parquet'

    output:
    path 'duckdb.endpoint'

    script:
    """
    duckdb analytics.db <<'SQL'
      INSTALL quack; LOAD quack;
      CREATE TABLE events AS SELECT * FROM read_parquet('raw/*.parquet');
      CALL quack_start(host := '127.0.0.1', port := 5432);
    SQL &
    DUCK_PID=\$!
    until nc -z 127.0.0.1 5432; do sleep 0.2; done
    echo "127.0.0.1:5432" > duckdb.endpoint
    wait \$DUCK_PID
    """
}

process daily_rollup { /* SELECT date_trunc('day', ts), count(*) ... */ }
process per_user_funnel { /* WITH steps AS (...) ... */ }
process anomaly_scan { /* SELECT ... WHERE z_score > 3 */ }

workflow {
    endpoint = duckdb_server(channel.fromPath('raw/*.parquet').collect())
    daily_rollup(endpoint)
    per_user_funnel(endpoint)
    anomaly_scan(endpoint)
}
```

All three analyses run in parallel against a shared in-process database, with the server torn down deterministically when they all finish. No long-lived state survives the pipeline run, which is the right default for reproducible workflows.

## Future work

The MVP intentionally restricts `service true` to the local executor. The directive is designed so that broader executor support can be added incrementally without breaking the user-facing surface (`service true` stays a single directive; the differences are in how the engine schedules, watches, and tears down the service task).

This section sketches the design for the three broader environments. None of it ships in the MVP, but the MVP code paths are factored so each can be layered on later.

### 1. Grid schedulers (Slurm, LSF, SGE, PBS)

These are the closest cousins of the local case because they typically have a shared POSIX work directory (NFS, Lustre). The three things that need to change versus the MVP:

- **Endpoint shape.** With a shared FS, a Unix domain socket no longer reaches a consumer on a different node. The declared output should be a connection string (e.g. `node05.cluster:5432`) written to a file in the work directory, and the service script should `bind` to `0.0.0.0` (or a specific interface) rather than `127.0.0.1`. This is a documentation/example change, not an engine change.
- **Co-location vs cluster-network access.** Two strategies:
    - *Sub-allocation*: schedule the service and its consumers inside a single multi-step job allocation (`salloc --nodes=N` + per-step `srun`). Lowest latency, fewest moving parts, but requires Nextflow to learn the concept of a "service allocation" and route the consumer submissions into it.
    - *Cluster-network access*: schedule the service as a normal job, let it publish `host:port` to its output file, let consumers reach it over the cluster network. Simpler to implement (no new scheduling primitive), but depends on the cluster permitting inter-node TCP and on the consumer knowing the network is reachable.

  The recommended first step is the cluster-network approach behind an opt-in `service.colocate` strategy directive (default `'network'`, future `'allocation'`).
- **Termination.** Grid `TaskHandler.killTask()` implementations already issue `scancel` / `bkill` / `qdel`. The `ServiceLifecycleObserver` works unchanged on top of these.

The watcher (`startServiceOutputWatcher`) is filesystem-only and works as-is over shared POSIX.

### 2. Kubernetes

Kubernetes is arguably the cleanest target after local because it natively models the service-and-consumers relationship.

- **Service task** becomes a `Pod` (or `Deployment` with replicas=1) plus a `Service` resource that gives it a stable in-cluster DNS name. The Pod's readiness probe handles "is the output ready" without polling — the engine watches the Pod's readiness condition via the K8s API instead of polling the work directory.
- **Endpoint** is the DNS name (`svc-<run-id>-<process>.<namespace>.svc.cluster.local:<port>`), written to the declared output file in the work directory so consumers can mount it.
- **Consumer tasks** are normal Pods that resolve that DNS name. No special scheduling — the Pod network handles connectivity.
- **Termination** deletes the Service and Pod via the K8s API.

This requires extensions to `nf-k8s`: a new `K8sServiceTaskHandler` that creates `Service` + `Pod`, and a new branch in `ServiceLifecycleObserver` (or executor-specific readiness signaling) that uses readiness probes instead of work-directory polling.

### 3. Cloud batch (AWS Batch, GCP Batch, Azure Batch)

The hardest case, because cloud batch jobs run in independent containers with no shared filesystem and no built-in cross-job discovery.

- **Endpoint propagation.** The service task's endpoint can't be written to a shared local path. Options:
    - Write `endpoint.txt` to the workflow's object-storage work directory; consumers read it from there. Works today but only signals readiness on the next consumer-poll cycle (high-ish latency).
    - Use a cloud-native service-discovery layer (AWS Cloud Map, GCP Service Directory, Azure DNS Private Zones). Cleaner but ties the executor implementation to provider-specific APIs.
- **Reachability.** Service container and consumer containers must be in the same VPC/subnet with the right security group / firewall rules. This is policy that the user (or a higher-level platform like Seqera Platform) controls; the directive can only document the requirement.
- **Readiness watcher.** Object-storage poll is the only portable option absent a callback channel. A future "callback URL" mechanism (the service POSTs to a Nextflow-provided URL when ready) would remove the poll latency.
- **Termination.** Batch `TaskHandler.killTask()` implementations already call `TerminateJob`/equivalent and work for service tasks unchanged.

For cloud batch, the recommended initial path is the object-storage `endpoint.txt` approach behind a documented "your service must be reachable in the workflow VPC" requirement, followed later by integration with Cloud Map / Service Directory and a callback channel.

### 4. Cross-cutting follow-ups

- **`service.colocate` directive** to choose the co-location strategy per executor where multiple strategies exist (`local`, `allocation`, `network`, `discovery`).
- **Group-local services** (one service instance per job group), once Nextflow has a stable job-group abstraction. Maps to Snakemake's `foo.{groupid}.socket` use case.
- **Readiness probes beyond "file appeared"**: a `serviceReady` sub-directive accepting an HTTP probe URL, a TCP `host:port`, or a shell command. Useful when the file is created before the server is actually accepting connections.
- **Multi-instance services** (sharded inference, read replicas) with a directive like `serviceReplicas N`. Each replica gets its own endpoint; consumer scheduling is round-robin or user-supplied.
- **Graceful drain**: on termination, send `SIGTERM` and wait up to a configurable timeout for the service to flush before sending `SIGKILL`.

## Links

- Snakemake service rules: <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#service-rules>
- Related directive: [hints](20260323-hints-process-directive.md) — same pattern of small, opt-in process directives
- NVIDIA NIM: <https://docs.nvidia.com/nim/>
- DuckDB quack extension: <https://duckdb.org/quack/>
