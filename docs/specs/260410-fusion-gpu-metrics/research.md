# Research: Fusion GPU Metrics Collection

## R1: How to detect Fusion at trace collection time

**Decision**: Use `task.processor.executor.isFusionEnabled()` in `TaskHandler.getTraceRecord()`.

**Rationale**: TaskHandler already accesses the executor at line 222 (`task.processor.executor.getName()`), so this is a proven access path. The base `Executor.isFusionEnabled()` returns `false` by default, and Fusion-capable executors override it via `FusionHelper.isFusionEnabled(session)`. This works for all handler subclasses without requiring `instanceof` checks.

**Alternatives considered**:
- Checking `this instanceof FusionAwareTask`: Would miss custom executors that support Fusion but don't implement the trait. Also, `FusionAwareTask` is a trait on handler subclasses, not on the base `TaskHandler` where `getTraceRecord()` lives.
- Adding a Fusion flag to TaskRun/TaskConfig: Unnecessary complexity — Fusion is an executor-level property, not a per-task property.

## R2: Where to read `.fusion/trace.json`

**Decision**: Read it in `TaskHandler.getTraceRecord()`, right after the existing `.command.trace` parsing block (lines 244-253), gated by `task.processor.executor.isFusionEnabled()`.

**Rationale**: This is the single place where all task trace data is assembled, regardless of executor type. The existing `.command.trace` parsing already demonstrates the pattern: resolve a file in the work dir, parse it, handle `NoSuchFileException` and `IOException` gracefully.

**Alternatives considered**:
- Reading in each TaskHandler subclass: Would require changes across 7 handler subclasses in both core and plugins. Much higher blast radius.
- Reading in `TowerObserver`: Would couple Platform-specific code with file I/O. The observer should only transform data, not collect it.

## R3: Transient field pattern on TraceRecord

**Decision**: Add `transient private Map<String,Object> gpuMetrics` with getter/setter, following the exact `resourceAllocation` pattern.

**Rationale**: Transient fields on TraceRecord are the established mechanism for carrying executor-specific data to TowerObserver without persisting it in serialization (Kryo). The `resourceAllocation` field is the closest precedent — it's also a `Map<String,Object>` set during trace collection and consumed in `TowerObserver.makeTaskMap0()`.

**Implementation details**:
- Field: `transient private Map<String,Object> gpuMetrics`
- Getter: `Map<String,Object> getGpuMetrics()`
- Setter: `void setGpuMetrics(Map<String,Object> value)`
- In `makeTaskMap0()`: `record.gpuMetrics = trace.getGpuMetrics()`

## R4: JSON parsing approach

**Decision**: Use Groovy's `JsonSlurper` to parse `.fusion/trace.json` and extract the `gpu` key.

**Rationale**: `JsonSlurper` is already used throughout the Nextflow codebase (e.g., in tests and utilities). It parses JSON into native Groovy maps/lists, which is exactly what we need for the `Map<String,Object>` transient field. No additional dependencies required.

## R5: Test strategy

**Decision**: Three test locations following existing patterns.

1. **TraceRecordTest**: Verify `gpuMetrics` transient field is not persisted across serialization (follows `numSpotInterruptions` test pattern).
2. **TraceRecordTest**: Verify `parseFusionTraceFile()` correctly extracts GPU block from valid JSON, handles missing file, handles malformed JSON, handles missing GPU block.
3. **TowerClientTest**: Verify `gpuMetrics` is included in task map output (follows `resourceAllocation` test at lines 684-711).

**Rationale**: These three test locations mirror exactly how `resourceAllocation` and `numSpotInterruptions` are tested, ensuring consistency with project conventions.
