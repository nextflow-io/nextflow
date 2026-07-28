# Tasks: Fusion GPU Metrics Collection

**Input**: Design documents from `/specs/260410-fusion-gpu-metrics/`
**Prerequisites**: plan.md, spec.md, research.md

**Tests**: Included — the spec requires unit tests for both TraceRecord and TowerClient.

**Organization**: Tasks grouped by user story for independent implementation and testing.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2)
- Exact file paths included in descriptions

## Phase 1: Foundational (TraceRecord transient field)

**Purpose**: Add the `gpuMetrics` transient field to TraceRecord — all subsequent tasks depend on this.

- [ ] T001 Add transient `gpuMetrics` field with getter/setter to `modules/nextflow/src/main/groovy/nextflow/trace/TraceRecord.groovy` (after `resourceAllocation` field at line 128, getter/setter after line 649)
- [ ] T002 Add static `parseFusionTraceFile(Path)` method to `modules/nextflow/src/main/groovy/nextflow/trace/TraceRecord.groovy` — parse `.fusion/trace.json` and return the `gpu` block as `Map<String,Object>`

**Checkpoint**: TraceRecord can hold and parse GPU metrics. No behavior change yet.

---

## Phase 2: User Story 1 - GPU metrics sent to Platform (Priority: P1)

**Goal**: Read `.fusion/trace.json` on task completion, extract GPU block, send to Platform via TowerObserver.

**Independent Test**: Run a Fusion-enabled task with `.fusion/trace.json` containing a `gpu` block, verify GPU data appears in the Platform task payload.

### Implementation

- [ ] T003 [US1] Read `.fusion/trace.json` in `TaskHandler.getTraceRecord()` at `modules/nextflow/src/main/groovy/nextflow/processor/TaskHandler.groovy` — add after `.command.trace` parsing block (after line 253), gated by `task.processor.executor.isFusionEnabled()`
- [ ] T004 [US1] Include `gpuMetrics` in task payload in `TowerObserver.makeTaskMap0()` at `plugins/nf-tower/src/main/io/seqera/tower/plugin/TowerObserver.groovy` — add `record.gpuMetrics = trace.getGpuMetrics()` after `resourceAllocation` line (line 476)

### Tests

- [ ] T005 [P] [US1] Test `parseFusionTraceFile` with valid GPU block in `modules/nextflow/src/test/groovy/nextflow/trace/TraceRecordTest.groovy` — create temp file with full trace.json, verify returned map has all GPU fields
- [ ] T006 [P] [US1] Test `gpuMetrics` transient field is not persisted across serialization in `modules/nextflow/src/test/groovy/nextflow/trace/TraceRecordTest.groovy` — set field, serialize/deserialize, verify null
- [ ] T007 [US1] Test `gpuMetrics` included in task map in `plugins/nf-tower/src/test/io/seqera/tower/plugin/TowerClientTest.groovy` — create TraceRecord with gpuMetrics set, call `makeTasksReq()`, verify output contains GPU data

**Checkpoint**: GPU metrics flow end-to-end from `.fusion/trace.json` to Platform payload. Run:
```bash
./gradlew :nextflow:test --tests "TraceRecordTest"
./gradlew :nf-tower:test --tests "TowerClientTest"
```

---

## Phase 3: User Story 2 - Graceful error handling (Priority: P2)

**Goal**: Ensure missing, malformed, or GPU-less trace files don't break task reporting.

**Independent Test**: Simulate tasks with missing/malformed `.fusion/trace.json`, verify task trace is sent without GPU data and no errors.

### Tests

- [ ] T008 [P] [US2] Test `parseFusionTraceFile` without GPU block in `modules/nextflow/src/test/groovy/nextflow/trace/TraceRecordTest.groovy` — create temp file with valid JSON but no `gpu` key, verify null returned
- [ ] T009 [P] [US2] Test `parseFusionTraceFile` with malformed JSON in `modules/nextflow/src/test/groovy/nextflow/trace/TraceRecordTest.groovy` — create temp file with invalid JSON, verify exception is thrown

**Checkpoint**: Error handling verified. The implementation in T003 already handles these cases via try/catch — these tests confirm the behavior.

---

## Phase 4: Verification

**Purpose**: End-to-end validation across both modules.

- [ ] T010 Run smoke tests to verify no regressions: `make smoke`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1** (T001, T002): No dependencies — start immediately
- **Phase 2** (T003-T007): Depends on Phase 1 completion
- **Phase 3** (T008-T009): Depends on Phase 1 (T002 specifically)
- **Phase 4** (T010): Depends on all previous phases

### Parallel Opportunities

- T001 and T002 modify the same file but different sections — execute sequentially
- T005, T006 are [P] — can run in parallel (same file but independent test methods)
- T008, T009 are [P] — can run in parallel
- T004 and T005/T006 are in different modules — can run in parallel after T001

### Within Each Phase

```
Phase 1:  T001 → T002
Phase 2:  T003 → T004 (sequential: different modules but T004 depends on field from T001)
          T005, T006 (parallel, after T002)
          T007 (after T004)
Phase 3:  T008, T009 (parallel, after T002)
Phase 4:  T010 (after all)
```

---

## Implementation Strategy

### MVP (User Story 1 Only)

1. Complete Phase 1: TraceRecord field + parser (T001-T002)
2. Complete Phase 2: TaskHandler + TowerObserver + tests (T003-T007)
3. **STOP and VALIDATE**: Run unit tests for both modules
4. GPU metrics flow to Platform

### Full Feature

1. MVP above
2. Add Phase 3: Error handling tests (T008-T009)
3. Phase 4: Smoke tests (T010)

---

## Summary

| Metric | Value |
|--------|-------|
| Total tasks | 10 |
| US1 tasks | 5 (T003-T007) |
| US2 tasks | 2 (T008-T009) |
| Foundational | 2 (T001-T002) |
| Verification | 1 (T010) |
| Files modified | 5 |
| Parallel opportunities | T005+T006, T008+T009 |
