# Propagate the Seqera Intelligent Compute scheduler run identifier to Platform

- Authors: Jorge Ejarque
- Status: draft
- Deciders: Jorge Ejarque, Paolo Di Tommaso
- Date: 2026-06-09
- Tags: nf-seqera, nf-tower, intelligent-compute, scheduler, platform, trace, accounting

## Summary

When a workflow runs on the Seqera Intelligent Compute scheduler (the `seqera` executor),
the scheduler assigns a run identifier that Platform currently has no way to associate with
its workflow record. This ADR covers the **Nextflow side**: capturing that run identifier and
propagating it to Platform on trace requests, plus exposing a config-derived flag that marks a
run as scheduler-managed.

## Problem Statement

Platform needs to know, per workflow, (a) whether the run was executed via the Intelligent
Compute scheduler and (b) the scheduler's run identifier. These are required to retrieve
authoritative cost and resource-usage metrics (settled cost, VM types, etc.) from the scheduler
after completion. Today neither value reaches Platform:

- The scheduler `runId` is created inside `SeqeraExecutor` (lazily, on the first task
  submission) and never leaves the executor.
- There is no "scheduler enabled" marker analogous to the existing Wave/Fusion flags.

## Goals or Decision Drivers

- Make the run identifier available to Platform **as early as possible** — the first task
  submission both assigns the id and drives the trace sender loop, so it is reported (via the
  `PATCH /workflow/{id}` endpoint) at the earliest moment it exists, well before completion.
- Mark scheduler-managed runs with an `enabled` flag available at workflow create/begin time,
  mirroring the existing Wave and Fusion metadata.
- Keep `nextflow` core free of any dependency on the `nf-seqera` plugin.
- Reuse existing patterns (metadata value classes, the `platform`/`scheduler` shared metadata,
  the `TowerJsonGenerator` field-stripping hook) rather than introducing new mechanisms.

## Non-goals

- The Platform-side changes (the `scheduler_enabled` column, the `tw_workflow_ext` satellite table,
  the `PATCH /workflow/{workflowId}` endpoint, and the cost-retrieval cron). Those are the Platform
  half of this work.
- Per-task cost or VM attribution and workspace-level rollups (out of scope for the consumer
  feature's v1).

## Solution or decision outcome

Add a `SchedulerMetadata` value class to core, attached to `WorkflowMetadata` as `scheduler`
(peer of `wave`/`fusion`). It carries:

- `enabled` — config-derived (the run-level executor resolves to `seqera`), serialized on the
  workflow object of **begin/complete** requests via `WorkflowMetadata.toMap()`.
- `runId` — set by `SeqeraExecutor` on the first task submission. Delivered via a dedicated
  `PATCH /workflow/{workflowId}` request (the generic workflow-extension endpoint) once assigned,
  and — once assigned — also on the `workflow.scheduler` object of the **complete** request (via
  `toMap()`'s omit-when-null rule).

## Rationale & discussion

### Why the `enabled` flag is config-derived

The lifecycle is `onFlowCreate` → `onFlowBegin` → `callIgniters()` (workflow runs, executors
register, tasks submit). The `seqera` executor is therefore **not registered yet** at
create/begin time, so the flag cannot be read from the executor. It is instead derived from the
configuration (`process.executor`, then `executor.name`, then `NXF_EXECUTOR`), exactly the way
`WaveMetadata`/`FusionMetadata` read their `enabled` state from config. Limitation: per-process
executor overrides (`withName:… { executor='seqera' }`) are not reflected — the flag describes
the run-level executor selection. The literal `'seqera'` is used to avoid a core → plugin
dependency.

The Trace **create** request is a minimal "hello" that does not carry the workflow object, so —
like Wave/Fusion — the flag rides on **begin** (and complete) via `toMap()`, not on create.

### How the `runId` is carried

`runId` is conceptually scheduler metadata, so it is a (`volatile`) field of `SchedulerMetadata`,
written by `SeqeraExecutor.createRun()` on the executor thread and read by the `TowerObserver`
reporting thread. To make the field reliably mutable and non-null for the executor, `scheduler`
uses a lazy getter mirroring `getPlatform()`.

`SchedulerMetadata` owns a `toMap()` that exposes `enabled` and includes `runId` **only when
assigned**. Rather than coupling the serializer to the field name, nf-tower registers a converter in
`TowerJsonGenerator.create()` — `addConverter(SchedulerMetadata){ m -> m.toMap() }`, exactly like the
existing `NextflowMeta` converter. Because the id is assigned on the first task submission, this
"omit-when-null" rule means it naturally falls out of every place the workflow object is serialized
*before* the id exists, and appears only where it does:

- **begin** — id unset → `scheduler: {enabled}` only.
- **complete** — id set → `scheduler: {enabled, runId}`, giving Platform a clean entity-level
  persistence point (the embedded scheduler meta, alongside `wave`/`fusion`).

Independently — and **earliest** — the `TowerObserver` delivers the id via a dedicated
`PATCH /workflow/{workflowId}` request (`client.updateWorkflow`) with body `{ schedulerRunId }`,
the generic workflow-extension update endpoint Platform exposes. It is sent **exactly once**: the
sender thread (`sendTasks0`) calls `sendSchedulerRunId()` each loop iteration, which fires the PATCH
the first time the id is non-null and then latches. The id is assigned on the first task submission —
the same trigger that starts the sender loop — so it is propagated at the earliest moment it exists,
well before completion. (The run id is stable for the life of the run, so re-sending is unnecessary;
Platform also treats an identical re-supply as idempotent and rejects a *conflicting* value with
`400`.)

The id therefore travels in two shapes — a top-level `schedulerRunId` via `PATCH /workflow/{id}`
(early signal, stored in Platform's workflow-extension satellite table), and nested
`workflow.scheduler.runId` on complete (final, entity-persistable) — read by different Platform
handlers.

### Lineage consumes the same `toMap()`

`nf-lineage` also builds its workflow-run record from `WorkflowMetadata.toMap()`, and it collects
at `onFlowBegin` (before the first task), so `runId` is unset there — and the record is fed into the
lineage execution hash. `LinObserver.collectWorkflowMetadata()` converts the `scheduler` entry via
`SchedulerMetadata.toMap()`, mirroring its existing `nextflow → toJsonMap()` handling; the
omit-when-null rule keeps the (then-unset) `runId` out of the lineage provenance.

### Resulting wire contract

```
# begin request (workflow object, via toMap) — id not yet assigned
workflow: { …, wave: {enabled}, fusion: {enabled, version}, scheduler: {enabled} }

# PATCH /workflow/{workflowId} — earliest carrier, sent once when the id is assigned
{ schedulerRunId: "run-xyz" }      // not sent when not scheduler-managed / not yet assigned

# complete request (workflow object, via toMap) — id assigned by now
workflow: { …, scheduler: {enabled, runId: "run-xyz"} }
```

### Components changed (Nextflow side)

- `modules/nextflow` — new `SchedulerMetadata`; `WorkflowMetadata.scheduler` field + lazy
  `getScheduler()`.
- `nf-seqera` — `SeqeraExecutor.createRun()` publishes `runId` to
  `session.workflowMetadata.scheduler`.
- `nf-tower` — `TowerObserver` sends the run id once via `PATCH /workflow/{workflowId}`
  (`TowerClient.updateWorkflow`); `TowerClient` gains the `updateWorkflow` method, its URL builder,
  and `PATCH` support in `makeRequest`. `TowerJsonGenerator` registers a `SchedulerMetadata → toMap()`
  converter.
- `nf-lineage` — `LinObserver.collectWorkflowMetadata()` converts the `scheduler` entry via
  `SchedulerMetadata.toMap()`.

### Platform-side dependency / risk

Platform stores the `enabled` flag as a plain `scheduler_enabled` column on the `Workflow` entity
(read from the begin request) and the `runId` in a `tw_workflow_ext` workflow-extension satellite
table written by the new `PATCH /workflow/{workflowId}` endpoint. Until that endpoint ships, the
PATCH request would be rejected; the `enabled` flag on the begin/complete `workflow.scheduler` object
relies on the receiver ignoring unknown JSON properties (Micronaut's default). This Nextflow change
should therefore be coordinated with — and land alongside — the Platform endpoint.

## Links

- The Platform-side cost-alignment feature is the consumer of this change; it requires the
  scheduler run identifier on the Platform side.
