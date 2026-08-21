# Automatic resource labels for any executor

- Authors: Paolo Di Tommaso
- Status: draft
- Deciders: Paolo Di Tommaso
- Date: 2026-08-21
- Tags: executor, labels, platform, seqera

## Summary

The `seqera.executor.autoLabels` option derives resource labels from workflow metadata, but only the
Seqera executor honours it. Promote that mapping to the runtime so that every executor supporting
`resourceLabels` — AWS Batch, Azure Batch, Google Batch, Kubernetes, Seqera — attaches the same
labels, driven by a new `tower.autoLabels` option.

## Problem Statement

`resourceLabels` attaches name-value pairs to the compute resources backing a task, for operational
purposes such as cloud cost attribution. Today the user must write those pairs by hand, and the
values that matter most for attribution — run name, session id, project, Platform workspace and
compute environment — are only known to Nextflow at runtime.

The Seqera executor already solves this. `seqera.executor.autoLabels` selects from thirteen
workflow-metadata fields, which `Labels.withWorkflowMetadata` maps to `nextflow.io/*` and
`seqera.io/platform/*` keys, and `SeqeraExecutor.createRun` attaches to the scheduler run object.
Per task, `SeqeraTaskHandler` sends only the delta against that run-level baseline.

That leaves a user running the same pipeline on AWS Batch, Azure Batch, Google Batch or Kubernetes
with no equivalent, even though those executors all consume resource labels through one accessor:

| Executor | Call site | Target |
| --- | --- | --- |
| AWS Batch | `AwsBatchTaskHandler:833` | Batch job tags |
| Azure Batch | `AzBatchService:755` | auto-pool metadata |
| Google Batch | `GoogleBatchTaskHandler:538,595` | allocation policy and job labels |
| Kubernetes | `K8sTaskHandler:298` | pod labels |
| Seqera | `SeqeraTaskHandler:163` | scheduler task labels |

Every one of them reads `task.config.getResourceLabels()`, so a single merge point in
`TaskConfig` reaches all five.

Three properties of the existing code shape the design:

1. **Injecting into the config scope does not work.** `resourceLabels` is not in
   `REPEATABLE_DIRECTIVES`, so `ProcessConfigBuilder.applyConfigDefaults` skips it whenever the
   process declares its own value, and a `withName:`/`withLabel:` selector replaces it outright. A
   plugin writing into `session.config.process.resourceLabels` would therefore lose the labels for
   exactly the processes that customise them.

2. **Label syntax is not portable.** `nextflow.io/runName` is a legal AWS tag key and a legal
   Kubernetes label key, but an illegal Google Cloud label key (lowercase `[a-z0-9_-]` only,
   must start with a letter). Values are worse: a Kubernetes label value rejects
   `https://github.com/foo/bar`, so `repository`, `revision` and `workflowUrl` cannot be applied
   verbatim. None of the four executors sanitise today; each passes the map straight to its API.

3. **Timing is already correct.** `session.start()` fires `notifyFlowCreate()` before the script is
   parsed, and `TowerObserver.onFlowCreate` populates `workflowMetadata.platform` there. Every
   process definition, and therefore every `TaskConfig`, is built afterwards.

## Goals or Decision Drivers

- One implementation of the metadata-to-label mapping, shared by all executors.
- No behaviour change for a user who does not opt in; the feature is off by default.
- User-declared labels are never altered — not rewritten, not reordered away, not dropped.
- Labels emitted to a cloud API must be valid for that API.
- No new task-hash inputs: enabling the feature must not invalidate a resumed run.
- Do not regress the Seqera executor's run-level baseline and per-task delta.

## Non-goals

- Platform *run* labels (`PlatformMetadata.labels`, fetched by `TowerObserver` as a
  `List<String>`). A different concept, untouched here.
- Enforcing per-cloud label cardinality ceilings (AWS 50 tags, Google 64 labels). Thirteen
  auto-labels plus user labels stays well under; the limits are documented, not policed.
- Making `resourceLabels` a repeatable directive.
- Recording resource labels in lineage metadata. Verified absent from `nf-lineage`, and the
  directive documentation's claim that they are not recorded remains accurate.

## Considered Options

- **A. Runtime computes the labels; nf-tower only supplies metadata.**
- **B. New plugin extension point** (`ResourceLabelsProvider extends ExtensionPoint`).
- **C. Mutable session registry** pushed from `TowerObserver.onFlowCreate`.

## Pros and Cons of the Options

### A. Runtime computes, plugin supplies metadata

Both halves are already in the runtime: `PlatformMetadata` is a runtime class that nf-tower fills
in at `onFlowCreate`, and `PlatformHelper.config()` already reads the `tower` config scope from the
runtime. A runtime helper can therefore read `tower.autoLabels` and map `WorkflowMetadata` to
labels with no new plugin machinery at all.

- Good, because it adds no service-provider interface for what has one producer.
- Good, because nf-tower's only change is declaring the config option.
- Good, because nf-seqera and every other executor consume one identical implementation.
- Good, because the core metadata labels still work when nf-tower is inactive.
- Bad, because the runtime reads a config key declared by a plugin scope — mitigated by the
  precedent of `tower.accessToken` and friends, already read this way through `PlatformHelper`.

### B. New plugin extension point

- Good, because third-party plugins could contribute labels.
- Good, because layering is explicit: the runtime declares, the plugin implements.
- Bad, because it is discovery, ordering and conflict-resolution machinery for one implementation.
- Bad, because nf-seqera must either implement the interface or consume it, adding indirection to
  code that already has the answer.
- Bad, because labels then depend on plugin load state.

### C. Mutable session registry

- Good, because there is no discovery mechanism and the push point is an explicit lifecycle event.
- Bad, because label content becomes dependent on observer ordering.
- Bad, because paths that never fire `onFlowCreate` (preview, inspect) silently yield no labels.

## Solution or decision outcome

Adopt option A: a runtime `AutoLabels` helper plus a merge in `TaskConfig`, with per-executor
sanitisation applied to the auto-derived entries only. nf-seqera drops its private copy of the
mapping and consumes the runtime one.

## Rationale & discussion

### Runtime: the mapping

New `nextflow.platform.AutoLabels` in `modules/nextflow`, holding what is today private to
nf-seqera:

- `VALID_NAMES` — the thirteen short names: `projectName`, `userName`, `runName`, `sessionId`,
  `resume`, `revision`, `commitId`, `repository`, `manifestName`, `runtimeVersion`, `workflowId`,
  `workspaceId`, `computeEnvId`.
- `parse(Object) -> Set<String>` — accepts `true` (all names), `false` (none), a list, or a
  comma-separated string; rejects unknown names with the existing error text. Moved from
  `ExecutorOpts.parseAutoLabels`.
- `labelsFor(WorkflowMetadata, Set<String>) -> Map<String,String>` — the canonical
  `nextflow.io/*` and `seqera.io/platform/*` mapping, including the `userName` fallback from the
  Platform user to the OS user. Moved from `Labels.withWorkflowMetadata`.

`Session` gains a memoized `getAutoResourceLabels()` that reads
`PlatformHelper.config().autoLabels`, parses it, and maps the metadata; it returns an empty map when
the option is unset. The computation is lazy on first access, not eager in `init()`: the Platform
fields are only populated at `notifyFlowCreate`, and any earlier snapshot would capture nulls.
Entries whose source value is absent are omitted, so a run without nf-tower active still gets the
metadata labels the runtime knows on its own.

### Runtime: the merge

`TaskConfig` holds no session reference, and reaching for `Global.session` inside it would make the
result depend on hidden global state. Instead `TaskConfig` gains an explicit transient field and
setter — the same shape as its existing `cache` field — assigned at the two sites in
`TaskProcessor` that already initialise `config.process` and `config.executor`:
`createTaskRun` (`TaskProcessor:758`) and `createTaskPreview` (`TaskProcessor:399`).

Two accessors result:

- `getResourceLabels()` — auto labels merged under declared labels, **declared wins on key
  collision**. This is what `TaskBean`, trace records and the execution report observe, so the
  report reflects what was actually applied.
- `getResourceLabels(ResourceLabelPolicy)` — the same merge, with the *auto* entries sanitised by
  the given policy and the declared entries passed through untouched. Executors call this overload.

The overload exists because the two requirements meet here: sanitisation is per-executor, but only
auto-labels may be sanitised. Once the maps are merged an executor can no longer tell one from the
other, so the distinction has to be preserved at the merge — hence the policy travels inward rather
than the auto subset travelling outward.

`resourceLabels` is absent from `TaskHasher`, so nothing here alters a task hash or invalidates
resume. It is already listed in `TaskArrayCollector.ARRAY_DIRECTIVES`, so job arrays inherit the
merged value from the first task, as they do today.

### Per-executor sanitisation

`ResourceLabelPolicy` (runtime) describes a key/value charset, case folding, maximum length, and
what to do with a value that sanitises to empty. One policy per executor, applied at its existing
call site:

| Executor | Policy | `nextflow.io/runName` becomes |
| --- | --- | --- |
| AWS Batch | permissive charset (`+ - = . _ : / @`), key <= 128, value <= 256 | unchanged |
| Google Batch | lowercase, `[a-z0-9_-]`, leading letter, <= 63 | `nextflow_io_runname` |
| Kubernetes | key unchanged; values stripped of scheme and slashes, <= 63 | unchanged (value fixed) |
| Azure Batch | near-identity, minus the reserved `microsoft` name prefix | unchanged |
| Seqera | identity | unchanged |

### nf-seqera convergence

nf-seqera stops owning the mapping:

- Delete `Labels.withWorkflowMetadata` and `ExecutorOpts.parseAutoLabels`; both now live in the
  runtime. `Labels.toStringMap` and `Labels.delta` stay — they are scheduler-specific.
- `SeqeraExecutor.createRun` builds the run label set from `session.autoResourceLabels` plus
  config-level `process.resourceLabels`.
- `runResourceLabels` becomes the **complete** run baseline (auto plus config) rather than
  config-only. This is load-bearing: `Labels.delta` compares the task label set against that
  baseline, so leaving it config-only would make the runtime merge re-send all thirteen auto-labels
  on every task. With the full baseline the delta collapses to empty unless the task genuinely
  overrides something.

### Config surface and precedence

- `tower.autoLabels`, declared in `TowerConfig` with `@ConfigOption` and `@Description`, accepting
  the same forms as the existing option: `true`, `false` (default), a list, or a comma-separated
  string.
- `seqera.executor.autoLabels` is retained for backward compatibility and marked deprecated in its
  description.
- One set is resolved per session, with `seqera.executor.autoLabels` taking precedence when the key
  is present in the configuration (including when set to `false`) and `tower.autoLabels` otherwise.
  Off by default when neither is set.

Resolving a single session-wide set, rather than one per executor, is deliberate. It keeps
`Session.autoResourceLabels` the only source of truth, and it is what guarantees the nf-seqera
per-task delta collapses to empty: the run baseline and the merged task labels derive from the same
map. The cost is that a run mixing the Seqera executor with another one honours the legacy
Seqera-scoped option globally — acceptable for a deprecated option that defaults to off.

The runtime merge applies to the Seqera executor as well; it is not special-cased. Its per-task
delta absorbs the duplication, which is why the full-baseline change described above is a
prerequisite rather than an optimisation.

### Testing

- Runtime: `AutoLabelsTest` for parsing and mapping — largely the existing nf-seqera `LabelsTest`
  cases relocated; `TaskConfigTest` for merge precedence and the policy overload; `SessionTest`
  for lazy computation and for the no-Platform case.
- Per executor: policy unit tests covering the hostile inputs (a `repository` URL under the
  Kubernetes value rules, mixed-case and dotted keys under the Google rules), plus one handler test
  each asserting that auto-labels arrive sanitised and user labels arrive byte-identical.
- nf-seqera: existing tests updated for the moved code, plus one pinning that the per-task delta is
  empty when a task declares no labels of its own.

### Documentation

- `docs/reference/config/tower.mdx` — the new option.
- `docs/reference/config/seqera.mdx` — deprecation note on the existing option.
- `docs/reference/process/directives/resource-labels.mdx` — an auto-labels section and the
  per-executor key-mangling table above.

### Risks

**Azure pool churn.** `AzBatchService.specFromAutoPool` derives the pool id from
`CacheHelper.hasher([vmType.name, opts, metadata])`, where `metadata` is the resource label map.
Any label that varies per run — `runName`, `sessionId`, `workflowId` — therefore produces a fresh
auto-pool for every execution, and pools accumulate unless
`azure.batch.deletePoolsOnCompletion` is set. This is inherent to Azure applying labels at pool
rather than task granularity. Documented as a caveat, with the recommendation to select a stable
subset (for example `projectName`, `workspaceId`, `computeEnvId`) on Azure. Gating volatile labels
out of the Azure policy is a possible follow-up, deliberately not done here.

**Cloud tag visibility.** Enabling the option changes tags on cloud resources, which can affect
cost reports and any tag-based IAM condition. Mitigated by the feature being off by default.

## Links

- Refines [scheduler run identifier propagation](20260609-scheduler-run-identifier-propagation.md)
