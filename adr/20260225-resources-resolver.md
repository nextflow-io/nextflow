# Pre-flight Resource Optimization via ResourcesResolver

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2026-02-25
- Tags: scheduler, resources, plugin

## Summary

Add a plugin-based `ResourcesResolver` extension point that allows a scheduler optimization engine to modify task resource allocation (cpus, memory, time, disk) before execution. The Seqera scheduler plugin provides the concrete implementation calling the optimization API.

## Problem Statement

Tasks currently execute with user-declared resources or static `resourceLimits` caps. There is no mechanism for an external optimization engine to dynamically adjust resource allocation based on historical execution data, prediction models, or scheduling constraints before a task runs.

## Design

### Pattern

Follows the `ContainerResolver` plugin pattern:

- **`ResourcesResolver`** — interface extending `ExtensionPoint` in core (`modules/nextflow`)
- **`DefaultResourcesResolver`** — no-op implementation in core, returns `null`
- **`ResourcesResolverProvider`** — loads the active resolver via `Plugins.getPriorityExtensions()`
- **`SeqeraResourcesResolver`** — in `nf-seqera` plugin, calls the scheduler optimization API. Uses `@Priority(-10)` to override the default.

### Call site

`TaskProcessor.submitTask()` — after hash computation, before `executor.submit()`. This means:

- The cache hash reflects **declared** (original) resources, not optimized ones — caching stays deterministic based on user intent.
- All task metadata is fully resolved (context bound, closures evaluated, workDir set).
- Available to all executors, not just Seqera.

### Request/Response

The resolver receives a `ResourcesRequest` DTO containing:

- cpus, memory, disk, time
- accelerator (count, type)
- machineType
- resourceLimits (upper bounds)
- container image
- process name
- attempt number
- task index
- input files total size in bytes

The resolver returns a `ResourcesResponse` with the full optimized resource set (cpus, memory, disk, time). All fields always populated.

### Resource storage: overrides map on TaskConfig

`TaskConfig` gets an `overrides` map. When set, the resource getters (`getCpus()`, `getMemory()`, `getTime()`, `getDisk()`) check overrides first, falling back to the declared config values.

This preserves the original user-declared values in the config map (`config.get('cpus')` still returns the declared value) while downstream consumers (TraceRecord, TaskBean, executors) automatically see optimized values through the existing getters.

### Sync invocation

The resolver is called synchronously, matching the `ContainerResolver` pattern. Performance optimization (batching, caching by process name) is the responsibility of each implementation.

### Failure mode

On any exception, log a warning and fall back to declared resources. The optimizer is an enhancement, not a gate.

## Implementation

### 1. TaskConfig overrides

Add to `TaskConfig.groovy`:

- `private Map<String,Object> overrides` field
- `setOverrides(Map)` / `getOverrides()` methods
- Modify `getCpus()`, `getMemory()`, `getTime()`, `getDisk()` to check overrides first
- Copy overrides in `clone()`

### 2. Core extension point

Create in `modules/nextflow/src/main/groovy/nextflow/processor/`:

- `ResourcesResolver.groovy` — interface with `enabled()` and `resolve(ResourcesRequest)` methods
- `ResourcesRequest.groovy` — DTO with all task metadata fields listed above
- `ResourcesResponse.groovy` — DTO with cpus, memory, disk, time
- `DefaultResourcesResolver.groovy` — no-op, returns `null`
- `ResourcesResolverProvider.groovy` — loads via `Plugins.getPriorityExtensions()`

Register `DefaultResourcesResolver` in `modules/nextflow/src/main/resources/META-INF/extensions.idx`.

### 3. Wire into TaskProcessor

In `TaskProcessor.submitTask()`, after setting hash/workDir/name and before `executor.submit()`:

```groovy
resolveResources(task)
```

The `resolveResources(task)` method:
1. Loads the resolver via `ResourcesResolverProvider.load()`
2. Builds a `ResourcesRequest` from `task.config` and `task` metadata
3. Calls `resolver.resolve(request)`
4. On success, applies `response` as overrides on `task.config`
5. On failure, logs warning, continues with declared values

### 4. Seqera plugin implementation

Create `SeqeraResourcesResolver` in `plugins/nf-seqera/src/main/io/seqera/executor/`:

- Implements `ResourcesResolver` with `@Priority(-10)`
- Uses `SchedClient` to call the scheduler optimization API
- Maps `ResourcesRequest` to the scheduler API request format
- Maps the API response back to `ResourcesResponse`
- Returns `null` on API errors (fallback to declared)
- `enabled()` returns `true` only when the scheduler client is configured

Register in the nf-seqera plugin's `extensions.idx`.

### 5. Verify downstream consumers

No changes needed — `TaskHandler.getTraceRecord()`, `TaskBean` constructor, and executor submit methods all call `task.config.getCpus()` / `getMemory()` which will automatically return overridden values.

## File Summary

| File | Action |
|------|--------|
| `modules/nextflow/.../processor/TaskConfig.groovy` | Modify — add overrides |
| `modules/nextflow/.../processor/ResourcesResolver.groovy` | Create — interface |
| `modules/nextflow/.../processor/ResourcesRequest.groovy` | Create — request DTO |
| `modules/nextflow/.../processor/ResourcesResponse.groovy` | Create — response DTO |
| `modules/nextflow/.../processor/DefaultResourcesResolver.groovy` | Create — no-op default |
| `modules/nextflow/.../processor/ResourcesResolverProvider.groovy` | Create — provider |
| `modules/nextflow/.../processor/TaskProcessor.groovy` | Modify — call resolver in submitTask() |
| `modules/nextflow/.../processor/TaskRun.groovy` | Modify — add applyResourceOverrides() |
| `modules/nextflow/src/main/resources/META-INF/extensions.idx` | Modify — register default |
| `plugins/nf-seqera/.../executor/SeqeraResourcesResolver.groovy` | Create — Seqera impl |
| nf-seqera `extensions.idx` | Modify — register Seqera impl |
