# Automatic resource labels for any executor

- Authors: Paolo Di Tommaso
- Status: accepted
- Date: 2026-08-21
- Tags: executor, labels, platform, seqera

## Summary

The `seqera.executor.autoLabels` option derives resource labels from workflow metadata, but only the Seqera executor honors it. Move that mapping into the runtime so that every executor supporting `resourceLabels` (AWS Batch, Azure Batch, Google Batch, Kubernetes, Seqera) attaches the same labels, driven by a new `tower.autoLabels` option.

## Problem statement

The `resourceLabels` directive attaches name-value pairs to the compute resources backing a task, for operational purposes such as cloud cost attribution. Today the user must write those pairs by hand, and the values that matter most for attribution (run name, session id, project, Platform workspace, compute environment) are only known to Nextflow at runtime.

The Seqera executor already solves this. `seqera.executor.autoLabels` selects from thirteen workflow-metadata fields, maps them to `nextflow.io/*` and `seqera.io/platform/*` keys, and attaches them to the scheduler run. Each task sends only the delta against that run-level baseline.

That leaves a user running the same pipeline on AWS Batch, Azure Batch, Google Batch or Kubernetes with no equivalent, even though those executors all consume resource labels through one accessor:

| Executor | Call site | Target |
| --- | --- | --- |
| AWS Batch | `AwsBatchTaskHandler:833` | Batch job tags |
| Azure Batch | `AzBatchService:755` | auto-pool metadata |
| Google Batch | `GoogleBatchTaskHandler:538,595` | allocation policy and job labels |
| Kubernetes | `K8sTaskHandler:298` | pod labels |
| Seqera | `SeqeraTaskHandler:163` | scheduler task labels |

## Goals

- One implementation of the metadata-to-label mapping, shared by all executors.
- No behavior change for a user who does not opt in; the feature is off by default.
- User-declared labels are never altered or dropped.
- Labels emitted to a cloud API must be valid for that API.
- No new task-hash inputs. Enabling the feature must not invalidate a resumed run.
- Do not regress the Seqera executor's run-level baseline and per-task delta.

## Non-goals

- Platform *run* labels. A different concept, untouched here.
- Enforcing per-cloud label cardinality ceilings (AWS 50 tags, Google 64 labels). Thirteen auto-labels plus user labels stays well under.
- Making `resourceLabels` a repeatable directive.
- Recording resource labels in lineage metadata.

## Considered options

### A. Runtime computes, plugin supplies metadata

Both halves are already in the runtime. `PlatformMetadata` is a runtime class that nf-tower fills in at `onFlowCreate`, and `PlatformHelper.config()` already reads the `tower` config scope from the runtime. A runtime helper can therefore read `tower.autoLabels` and map `WorkflowMetadata` to labels with no new plugin machinery at all.

- Good, because it adds no service-provider interface for what has one producer.
- Good, because nf-tower's only change is declaring the config option.
- Good, because nf-seqera and every other executor consume one identical implementation.
- Good, because the core metadata labels still work when nf-tower is inactive.
- Bad, because the runtime reads a config key declared by a plugin scope. Mitigated by the precedent of `tower.accessToken` and friends, already read this way through `PlatformHelper`.

### B. New plugin extension point

- Good, because third-party plugins could contribute labels.
- Good, because layering is explicit: the runtime declares, the plugin implements.
- Bad, because it is discovery, ordering and conflict-resolution machinery for one implementation.
- Bad, because nf-seqera must either implement the interface or consume it, adding indirection to code that already has the answer.
- Bad, because labels then depend on plugin load state.

### C. Mutable session registry

- Good, because there is no discovery mechanism and the push point is an explicit lifecycle event.
- Bad, because label content becomes dependent on observer ordering.
- Bad, because paths that never fire `onFlowCreate` (preview, inspect) silently yield no labels.

## Solution

Adopt option A: a runtime `AutoLabels` helper plus a merge in `TaskConfig`, with per-executor sanitization applied to the auto-derived entries only. nf-seqera drops its private copy of the mapping and consumes the runtime one.

## Rationale and discussion

### Runtime: the mapping

A new `AutoLabels` class provides what is today private to nf-seqera:

- `VALID_NAMES`, the thirteen short names: `projectName`, `userName`, `runName`, `sessionId`, `resume`, `revision`, `commitId`, `repository`, `manifestName`, `runtimeVersion`, `workflowId`, `workspaceId`, `computeEnvId`.
- `parse(Object) -> Set<String>`, which accepts `true` (all names), `false` (none), a list, or a comma-separated string, and rejects unknown names with the existing error text.
- `labelsFor(WorkflowMetadata, Set<String>) -> Map<String,String>`, the canonical `nextflow.io/*` and `seqera.io/platform/*` mapping, including the `userName` fallback from the Platform user to the OS user.

`Session` gains a memoized `getAutoResourceLabels()` that loads the auto-labels config and maps the workflow metadata. The computation is lazy on first access, not eager in `init()`, because the Platform fields aren't populated until `notifyFlowCreate`. Entries with no source value are omitted, so a run without nf-tower still gets the metadata labels from the core runtime.

### Runtime: the merge

`TaskConfig` holds no session reference, and reaching for `Global.session` inside it would make the result depend on hidden global state. Instead `TaskConfig` gains an explicit transient field and setter, the same shape as its existing `cache` field.

Two accessors result:

- `getResourceLabels()`, auto labels merged under declared labels, where a declared label wins on key collision.
- `getResourceLabels(ResourceLabelPolicy)`, the same merge, with the *auto* entries sanitized by the given policy and the declared entries passed through untouched. Executors call this overload. The collision is decided *before* the sanitization, so that a label declared with the canonical key (`nextflow.io/runName`) overrides the auto one on a policy that mangles that key.

The overload exists because the two requirements meet here: sanitization is per-executor, but only auto-labels may be sanitized. Once the maps are merged an executor can no longer tell one from the other, so the distinction has to be preserved at the merge. The policy travels inward rather than the auto subset traveling outward.

### Per-executor sanitization

`ResourceLabelPolicy` describes a key/value charset, case folding, maximum length, and what to do with a value that sanitizes to empty. One policy per executor, applied at its existing call site:

| Executor | Policy | `nextflow.io/runName` becomes |
| --- | --- | --- |
| AWS Batch | permissive charset (`+ - = . _ : / @`), key <= 128, value <= 256 | unchanged |
| Google Batch | lowercase, `[a-z0-9_-]`, leading letter, <= 63 | `nextflow_io_runname` |
| Kubernetes | key prefix kept, any `/` after the first replaced with `_`; values stripped of scheme and slashes, <= 63 | unchanged (value fixed) |
| Azure Batch | near-identity, minus the reserved `microsoft` name prefix | unchanged |
| Seqera | identity | unchanged |

### nf-seqera convergence

The mapping logic is moved from nf-seqera to the core runtime. The Seqera executor builds the run label set from `session.autoResourceLabels` plus config-level `process.resourceLabels`. This way, the auto labels are not re-sent by the task delta.

### Config surface and precedence

`tower.autoLabels` accepts the same forms as the existing option: `true`, `false` (default), a list, or a comma-separated string.

`seqera.executor.autoLabels` is retained for backward compatibility and marked as deprecated. It takes precedence over `tower.autoLabels` when present.

Resolving a single session-wide set, rather than one per executor, guarantees the nf-seqera per-task delta collapses to empty. The run baseline and the merged task labels derive from the same map. The consequence is that `seqera.executor.autoLabels` now applies globally, which is acceptable for a deprecated option that is off by default.

### Risks

Azure pool churn. `AzBatchService.specFromAutoPool` derives the pool id from `CacheHelper.hasher([vmType.name, opts, metadata])`, where `metadata` is the resource label map. Any label that varies per run (`runName`, `sessionId`, `workflowId`) therefore produces a fresh auto-pool for every execution, and pools accumulate unless `azure.batch.deletePoolsOnCompletion` is set. This is inherent to Azure applying labels at pool rather than task granularity. Documented as a caveat, with the recommendation to select a stable subset (for example `projectName`, `workspaceId`, `computeEnvId`) on Azure. Gating volatile labels out of the Azure policy is a possible follow-up, deliberately not done here.

Cloud tag visibility. Enabling the option changes tags on cloud resources, which can affect cost reports and any tag-based IAM condition. Mitigated by the feature being off by default.

## Links

- Refines [scheduler run identifier propagation](20260609-scheduler-run-identifier-propagation.md)
