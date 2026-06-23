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
propagating it to Platform via a dedicated `PATCH /workflow/{workflowId}` request.

## Problem Statement

Platform needs the scheduler's run identifier, per workflow, to retrieve authoritative cost and
resource-usage metrics (settled cost, VM types, etc.) from the scheduler after completion. Today
the scheduler `runId` is created inside `SeqeraExecutor` (lazily, on the first task submission) and
never leaves the executor.

## Goals or Decision Drivers

- Make the run identifier available to Platform **as early as possible** — the first task
  submission both assigns the id and drives the trace sender loop, so it is reported (via the
  `PATCH /workflow/{id}` endpoint) at the earliest moment it exists, well before completion.
- Keep `nextflow` core free of any dependency on the `nf-seqera` plugin.
- Reuse existing patterns (`WorkflowMetadata`, the `TowerObserver` sender loop, `TowerClient`)
  rather than introducing new mechanisms.

## Non-goals

- The Platform-side changes (the `tw_workflow_ext` satellite table, the
  `PATCH /workflow/{workflowId}` endpoint, and the cost-retrieval cron). Those are the Platform
  half of this work.
- Per-task cost or VM attribution and workspace-level rollups (out of scope for the consumer
  feature's v1).
- A separate "scheduler enabled" flag. An earlier iteration carried a config-derived `enabled`
  boolean on the workflow object; it was dropped because (a) it is never consumed inside Nextflow,
  (b) deriving it from run-level config is unreliable (it misses per-process `withName:{ executor }`
  overrides, so it could report `false` for a run that does use the scheduler and sends a `runId`),
  and (c) the presence of `schedRunId` already tells Platform the run is scheduler-managed.

## Solution or decision outcome

The scheduler run id is held on `PlatformMetadata` as a single `volatile String schedRunId` field —
i.e. `workflow.platform.schedRunId`, grouped with the other Platform identifiers (`workflowId`,
`workflowUrl`) rather than promoted to a top-level `workflow.schedRunId` attribute. It is delivered to
Platform via a one-off `PATCH /workflow/{workflowId}` request with body `{ "schedRunId": "…" }`.

- `SeqeraExecutor.createRun()` assigns the run id on the first task submission and publishes it to
  `PlatformMetadata.schedRunId` (`session.workflowMetadata?.platform?.setSchedRunId(runId)`).
- `TowerObserver` reads it and sends the `PATCH` exactly once, at the earliest moment the id exists.

## Rationale & discussion

### Where the `runId` lives

The run id is held on `PlatformMetadata` (`workflow.platform.schedRunId`), next to the other
Platform identifiers, rather than as a top-level `WorkflowMetadata.schedRunId` field. This is a
deliberate choice:

- **It is not a general workflow property.** The run id is meaningless outside a Seqera Platform /
  Intelligent Compute run — it is a Platform-assigned identifier, conceptually a sibling of
  `platform.workflowId`. Promoting it to a top-level `workflow` attribute would advertise it to every
  pipeline as if it were portable run metadata (like `workflow.runName` or `workflow.sessionId`),
  which it is not.
- **It keeps the serialization behaviour coherent for free.** `TowerObserver.makeCompleteReq()`
  already strips the whole `platform` sub-object before sending, and at `onFlowBegin` (when the begin
  and lineage records capture `toMap()`) the id is always still null — it is assigned later, on the
  first task submission. So nesting under `platform` means the id is never carried as a real value on
  the begin/complete/lineage workflow object, **without** needing a special-case exclusion in
  `WorkflowMetadata.toMap()`. Platform receives the actual value only through the dedicated `PATCH`
  endpoint.

`schedRunId` is written by `SeqeraExecutor.createRun()` on the executor thread and read by the
`TowerObserver` reporting thread, so the field is `volatile` to publish that write across threads.

The `TowerObserver` delivers the id via `PATCH /workflow/{workflowId}` (`client.updateWorkflow`) with
body `{ schedRunId }`. It is sent **exactly once**: the sender thread (`sendTasks0`) calls
`sendSchedRunId()` each loop iteration, which fires the PATCH the first time the id is non-null and
then latches. The id is assigned on the first task submission — the same trigger that starts the
sender loop — so it is propagated at the earliest moment it exists, well before completion. (The run
id is stable for the life of the run, so re-sending is unnecessary; Platform also treats an identical
re-supply as idempotent and rejects a *conflicting* value with `400`.)

### Best-effort delivery

Propagating the run id is best-effort telemetry. `TowerClient.updateWorkflow` returns the HTTP
`Response` rather than throwing, and `sendSchedRunId()` logs a failed update and lets the run carry
on instead of aborting. Transient failures are already retried at the HTTP-client layer
(`retryConfig`), so the single attempt latches only once that retry policy is exhausted.

### Resulting wire contract

```
# PATCH /workflow/{workflowId} — sole carrier, sent once when the id is assigned
{ schedRunId: "run-xyz" }      // not sent when not scheduler-managed / not yet assigned
```

The begin/complete (and lineage) workflow object carries no scheduler run id: it is either stripped
with the `platform` sub-object (complete) or still null at capture time (begin/lineage).

### Components changed (Nextflow side)

- `modules/nextflow` — `PlatformMetadata.schedRunId` (`volatile`) field, grouped with the other
  Platform identifiers. No change to `WorkflowMetadata.toMap()`.
- `nf-seqera` — `SeqeraExecutor.createRun()` publishes `runId` to
  `session.workflowMetadata.platform.schedRunId`.
- `nf-tower` — `TowerObserver` sends the run id once via `PATCH /workflow/{workflowId}`
  (`TowerClient.updateWorkflow`, best-effort); `TowerClient` gains the `updateWorkflow` method, its
  URL builder, and `PATCH` support in `makeRequest`.

### Platform-side dependency / risk

Platform stores the `runId` in a `tw_workflow_ext` workflow-extension satellite table written by the
new `PATCH /workflow/{workflowId}` endpoint. Until that endpoint ships, the PATCH request would be
rejected — but because delivery is best-effort (the failure is logged and the run continues), this no
longer creates a hard Nextflow ↔ Platform version dependency. The change should still be coordinated
with the Platform endpoint so the id is actually persisted.

## Links

- The Platform-side cost-alignment feature is the consumer of this change; it requires the
  scheduler run identifier on the Platform side.
