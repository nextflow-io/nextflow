# Seqera executor: support `process.resourceLabels`

Date: 2026-04-17
Status: Approved

## Problem

The `nf-seqera` executor does not honour the `process.resourceLabels`
directive. `SeqeraTaskHandler.submit()` builds the scheduler `Task` with
`name`, `image`, `command`, `environment`, `resourceRequirement`,
`resourceLimit`, `machineRequirement`, and `nextflow(taskId/hash/workDir)` —
it never reads `task.config.getResourceLabels()`. The plugin's only label
path is at the run level (`SeqeraExecutor.createRun()`), where
`Labels.withUserLabels(seqeraConfig.labels)` and optional auto-labels are
attached to `CreateRunRequest`.

`AbstractComputePlatformProvider.addConfigResourceLabels()` emits
`process.resourceLabels = [...]` into the Nextflow config for every CE type
including the Seqera Compute default config. AWS Batch / GCP Batch / Azure /
K8s honour the directive; on the Seqera Compute path the directive is
effectively dead.

## Goal

Implement support for `process.resourceLabels` in the Seqera executor and
pass labels through to the `sched-client`, with cumulative semantics that
mirror Nextflow's existing label model.

## Label model

Nextflow labels are cumulative:

- `process.resourceLabels` at the top level of `nextflow.config` is the
  common baseline — it applies to every task across every process.
- Selector-scoped (`withName:`, `withLabel:`) and in-process-body
  `resourceLabels` directives merge on top, per process.
- `TaskConfig.getResourceLabels()` returns the final merged map for a given
  task.

The `sched-api` (≥ 0.51.0) exposes labels at two scopes:

- `CreateRunRequest.labels` — set once at run creation
- `Task.labels` — set per task

We map cumulative Nextflow labels onto these two scopes:

- **Run-level labels** = config-level `process.resourceLabels` (the common
  baseline) + `nextflow.io/*` auto-labels (when `seqera.executor.autoLabels`
  is enabled).
- **Per-task labels** = the *delta* between `task.config.getResourceLabels()`
  and the run-level baseline:
  - keys present on the task but absent from the run baseline
  - keys present in both where the task value differs from the run value
  - keys present in both with identical values are omitted

When the delta is empty, `Task.labels` is left unset.

The Sched scheduler is expected to merge run + task labels with task labels
overriding run labels on key collision; this preserves the Nextflow
semantic where a per-process `resourceLabels` directive overrides the
config-level default for the same key.

## Changes

### 1. Remove `seqera.executor.labels`

This config option becomes redundant once `process.resourceLabels` is the
canonical user-facing way to attach run-level labels.

- `ExecutorOpts` (`plugins/nf-seqera/src/main/io/seqera/config/ExecutorOpts.groovy`):
  remove the `labels` field, getter, and `@ConfigOption` declaration.
- `SeqeraExecutor.createRun()`: drop the
  `labels.withUserLabels(seqeraConfig.labels)` call.
- `Labels`: remove `withUserLabels(Map)` (no remaining callers).
- `ExecutorOptsTest`, `LabelsTest`, `SeqeraExecutorTest`: drop assertions
  for the removed option.
- Plugin `changelog.txt`: note the removal as a breaking change for plugin
  `nf-seqera` 0.18.0.

### 2. Add `withProcessResourceLabels` to `Labels`

`plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy`:

```groovy
Labels withProcessResourceLabels(Map<String,Object> map) {
    if( map )
        map.each { k, v -> entries.put(k.toString(), String.valueOf(v)) }
    return this
}
```

Values are coerced to `String` via `String.valueOf` to satisfy
`sched-api`'s `Map<String,String>` typing without rejecting non-string
values that Nextflow's `resourceLabels` directive accepts.

### 3. Wire run-level labels in `SeqeraExecutor.createRun()`

```groovy
final processLabels = (session.config.process as Map)?.resourceLabels as Map<String,Object>
final labels = new Labels()
if( seqeraConfig.autoLabels )
    labels.withWorkflowMetadata(session.workflowMetadata)
labels.withProcessResourceLabels(processLabels)
this.runResourceLabels = coerceToStringMap(processLabels)
```

The coerced map is cached on the executor as `runResourceLabels` so task
handlers can compute the delta without re-reading config or duplicating the
coercion logic. `coerceToStringMap` lives next to `Labels` (or as a
`static` helper on it) and applies `String.valueOf` to each value.

### 4. Compute and attach the per-task delta in `SeqeraTaskHandler.submit()`

```groovy
final taskLabels = coerceToStringMap(task.config.getResourceLabels())
final delta = deltaLabels(taskLabels, executor.runResourceLabels)
if( delta )
    schedTask.labels(delta)
```

`deltaLabels(task, run)` returns a `Map<String,String>` containing entries
in `task` that are missing in `run` or whose value differs from the value
in `run`. Empty map → return `null` so the caller can omit the field.

The helper lives alongside `Labels` (e.g. `Labels.delta(task, run)`).

### 5. `sched-client` dependency

- `plugins/nf-seqera/build.gradle`: bump `io.seqera:sched-client` to the
  released version exposing `Task.labels` (≥ 0.51.0 once published).
- `settings.gradle`: add an `includeBuild '../sched'` block matching the
  existing commented `includeBuild '../nextflow-plugin-gradle'` pattern, so
  development against an unreleased sched is opt-in via uncommenting.

### 6. Tests (Spock)

- `LabelsTest`:
  - `withProcessResourceLabels` merges entries, coerces non-String values,
    no-ops on null/empty.
  - `delta(task, run)` returns missing keys, returns differing keys, omits
    matching keys, returns empty/null when fully covered.
  - Existing `withUserLabels` assertions removed.
- `SeqeraExecutorTest`:
  - `createRun` populates `CreateRunRequest.labels` with config-level
    `process.resourceLabels` merged with auto-labels.
  - `runResourceLabels` accessor returns the coerced baseline map.
  - Removed assertions for `seqera.executor.labels`.
- `SeqeraTaskHandlerTest`:
  - `submit` attaches `Task.labels` containing the delta when the task adds
    new labels or overrides values.
  - `submit` leaves `Task.labels` unset when the task labels equal the run
    baseline.
- `ExecutorOptsTest`: remove `labels` parsing test.

### 7. Docs

- `docs/reference/process.md` (resourceLabels section, around line 1388):
  add `{ref}seqera-executor` to the list of executors that support
  `resourceLabels`.
- `plugins/nf-seqera/changelog.txt`: entry covering the new behaviour and
  the removal of `seqera.executor.labels`.

## Out of scope

- No deprecation warning shim for `seqera.executor.labels` — the plugin is
  early (0.17.0) and the user has approved removal.
- No changes to `seqera.executor.autoLabels` semantics or to the
  `nextflow.io/*` / `seqera:sched:*` label namespaces.
- No changes to the sched-api or sched-client itself; the `Task.labels`
  field is assumed already published in the version we depend on.
