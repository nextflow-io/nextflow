# Implementation Plan: Fusion GPU Metrics Collection

**Branch**: `260410-fusion-gpu-metrics-v2` | **Date**: 2026-04-10 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/260410-fusion-gpu-metrics/spec.md`

## Summary

Collect GPU metrics from Fusion's `.fusion/trace.json` file on task completion and send them to Seqera Platform. The GPU block is carried as a transient `Map<String,Object>` field on `TraceRecord` (same pattern as `resourceAllocation`) and included in the task payload via `TowerObserver.makeTaskMap0()`.

## Technical Context

**Language/Version**: Groovy 4.0.29 / Java 17 target (Java 21 toolchain)
**Primary Dependencies**: Nextflow core (`modules/nextflow`), nf-tower plugin (`plugins/nf-tower`)
**Storage**: N/A (read-only file access to `.fusion/trace.json`)
**Testing**: Spock Framework (unit tests in both modules)
**Target Platform**: All Fusion-enabled executors (AWS Batch, Google Batch, Azure Batch, K8s, Seqera, SLURM)
**Project Type**: Multi-module Gradle project
**Performance Goals**: Negligible overhead — one small JSON file read per task completion
**Constraints**: Must not break existing trace pipeline; must be forward-compatible with evolving GPU block schema
**Scale/Scope**: 4 files modified, ~80 lines of production code, ~120 lines of test code

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Modular Architecture | PASS | Core trace logic in `modules/nextflow`, Platform integration in `plugins/nf-tower` — correct placement |
| II. Test-Driven Quality | PASS | Unit tests planned for both TraceRecord and TowerClient |
| III. Dataflow Programming | N/A | No changes to dataflow model |
| IV. Apache 2.0 License | PASS | All modified files already have headers |
| V. DCO Sign-off | PASS | Will use `git commit -s` |
| VI. Semantic Versioning | PASS | No version bump needed — feature addition within existing release cycle |
| VII. Groovy Idioms | PASS | Uses JsonSlurper, follows existing getter/setter patterns |

## Project Structure

### Files to Modify

```text
modules/nextflow/
├── src/main/groovy/nextflow/trace/TraceRecord.groovy        # Add transient field + parsing method
├── src/main/groovy/nextflow/processor/TaskHandler.groovy     # Read .fusion/trace.json on completion
└── src/test/groovy/nextflow/trace/TraceRecordTest.groovy     # Test transient field + parsing

plugins/nf-tower/
├── src/main/io/seqera/tower/plugin/TowerObserver.groovy      # Include gpuMetrics in task map
└── src/test/io/seqera/tower/plugin/TowerClientTest.groovy    # Test gpuMetrics in task map
```

## Implementation Tasks

### Task 1: Add transient `gpuMetrics` field to TraceRecord

**File**: `modules/nextflow/src/main/groovy/nextflow/trace/TraceRecord.groovy`

**Changes**:
1. Add field after `resourceAllocation` (line 128):
   ```groovy
   transient private Map<String,Object> gpuMetrics
   ```
2. Add getter/setter after existing `resourceAllocation` getter/setter (after line 649):
   ```groovy
   Map<String,Object> getGpuMetrics() {
       return gpuMetrics
   }

   void setGpuMetrics(Map<String,Object> value) {
       this.gpuMetrics = value
   }
   ```

### Task 2: Add Fusion trace file parsing method to TraceRecord

**File**: `modules/nextflow/src/main/groovy/nextflow/trace/TraceRecord.groovy`

**Changes**:
Add a static method to parse `.fusion/trace.json` and extract the `gpu` block:
```groovy
static Map<String,Object> parseFusionTraceFile(Path file) {
    final text = file.text
    final json = (Map) new JsonSlurper().parseText(text)
    return (Map<String,Object>) json.get('gpu')
}
```

This keeps parsing logic on TraceRecord (consistent with `parseTraceFile()` for `.command.trace`).

### Task 3: Read `.fusion/trace.json` in TaskHandler.getTraceRecord()

**File**: `modules/nextflow/src/main/groovy/nextflow/processor/TaskHandler.groovy`

**Changes**:
After the existing `.command.trace` parsing block (after line 253), add:
```groovy
// collect Fusion GPU metrics
if( task.processor.executor.isFusionEnabled() ) {
    final fusionTrace = task.workDir?.resolve('.fusion/trace.json')
    try {
        if( fusionTrace ) {
            final gpu = TraceRecord.parseFusionTraceFile(fusionTrace)
            if( gpu )
                record.gpuMetrics = gpu
        }
    }
    catch( NoSuchFileException e ) {
        // ignore - Fusion trace may not exist
    }
    catch( Exception e ) {
        log.debug "[WARN] Cannot read Fusion trace file: $fusionTrace -- Cause: ${e.message}"
    }
}
```

**Key design decisions**:
- Gated by `task.processor.executor.isFusionEnabled()` — no file access when Fusion is not enabled (FR-007)
- Placed inside `isCompleted()` block but NOT gated by task status — runs for both success and failure (FR-005)
- Same error handling pattern as `.command.trace` parsing above it (FR-006)

### Task 4: Include `gpuMetrics` in TowerObserver task payload

**File**: `plugins/nf-tower/src/main/io/seqera/tower/plugin/TowerObserver.groovy`

**Changes**:
In `makeTaskMap0()` method, add after `record.resourceAllocation = trace.getResourceAllocation()` (after line 476):
```groovy
record.gpuMetrics = trace.getGpuMetrics()
```

### Task 5: Unit tests for TraceRecord

**File**: `modules/nextflow/src/test/groovy/nextflow/trace/TraceRecordTest.groovy`

**Tests to add**:

1. **Transient field serialization test** (follows `numSpotInterruptions` pattern):
   - Set `gpuMetrics` on a TraceRecord
   - Serialize and deserialize
   - Verify deserialized record has `null` for `gpuMetrics`

2. **parseFusionTraceFile with GPU block**:
   - Create a temp file with valid trace.json content including a `gpu` block
   - Verify the returned map contains all GPU fields with correct values

3. **parseFusionTraceFile without GPU block**:
   - Create a temp file with valid trace.json content without a `gpu` key
   - Verify `null` is returned

4. **parseFusionTraceFile with malformed JSON**:
   - Create a temp file with invalid JSON
   - Verify an exception is thrown (caller handles it)

### Task 6: Unit tests for TowerClient/TowerObserver

**File**: `plugins/nf-tower/src/test/io/seqera/tower/plugin/TowerClientTest.groovy`

**Test to add** (follows `resourceAllocation` test at lines 684-711):
- Create a TraceRecord with `gpuMetrics` set to a GPU metrics map
- Call `makeTasksReq([trace])`
- Verify `req.tasks[0].gpuMetrics` contains the GPU data

## Implementation Order

1. **Task 1 + Task 2** (TraceRecord changes) — no dependencies
2. **Task 3** (TaskHandler) — depends on Task 1+2
3. **Task 4** (TowerObserver) — depends on Task 1
4. **Task 5** (TraceRecord tests) — depends on Task 1+2
5. **Task 6** (TowerClient tests) — depends on Task 4

Tasks 1+2 and 5 can be done in parallel with Tasks 4 and 6.

## Verification

After implementation, run:
```bash
./gradlew :nextflow:test --tests "TraceRecordTest"
./gradlew :nf-tower:test --tests "TowerClientTest"
make smoke  # verify no regressions
```
