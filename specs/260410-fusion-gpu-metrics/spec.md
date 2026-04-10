# Feature Specification: Fusion GPU Metrics Collection

**Feature Branch**: `260410-fusion-gpu-metrics`  
**Created**: 2026-04-10  
**Status**: Draft  
**Input**: User description: "Collect GPU metrics from Fusion trace.json and send to Seqera Platform via TowerClient"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - GPU metrics sent to Platform on task completion (Priority: P1)

A user runs a Nextflow pipeline with Fusion enabled on a GPU-equipped executor (e.g., AWS Batch, Google Batch, Kubernetes). When each task completes, Nextflow reads the Fusion-generated `.fusion/trace.json` file from the task work directory, extracts the `gpu` block, and includes it in the task trace data sent to Seqera Platform. The user can then view GPU utilization metrics (compute %, memory %, active time, etc.) for each task in the Platform UI.

**Why this priority**: This is the core feature. Without it, GPU usage is invisible to Platform users running Fusion-enabled pipelines.

**Independent Test**: Can be tested by running a Fusion-enabled task that produces a `.fusion/trace.json` with a `gpu` block, then verifying the GPU data appears in the task payload sent to Platform.

**Acceptance Scenarios**:

1. **Given** a completed task with Fusion enabled and a valid `.fusion/trace.json` containing a `gpu` block, **When** the task trace is collected, **Then** all GPU metrics from the `gpu` block are included in the task data sent to Platform.
2. **Given** a completed task with Fusion enabled and a valid `.fusion/trace.json` without a `gpu` block (CPU-only task), **When** the task trace is collected, **Then** no GPU metrics are sent and no error occurs.
3. **Given** a failed task with Fusion enabled and a valid `.fusion/trace.json` containing a `gpu` block, **When** the task trace is collected, **Then** GPU metrics are still sent (metrics are collected irrespective of task status).

---

### User Story 2 - Graceful handling when trace.json is missing or malformed (Priority: P2)

When Fusion's `.fusion/trace.json` file is missing (e.g., task was killed before Fusion wrote it) or contains invalid JSON, the system logs a debug-level warning and proceeds without GPU metrics. The task trace is still sent to Platform with all other fields intact.

**Why this priority**: Robustness is essential — GPU metrics are supplementary data and must never cause task reporting to fail.

**Independent Test**: Can be tested by simulating a completed task where `.fusion/trace.json` is absent or contains malformed JSON, and verifying the task trace is still sent successfully without GPU data.

**Acceptance Scenarios**:

1. **Given** a completed Fusion-enabled task where `.fusion/trace.json` does not exist, **When** the task trace is collected, **Then** no GPU metrics are included and no error is raised.
2. **Given** a completed Fusion-enabled task where `.fusion/trace.json` contains invalid JSON, **When** the task trace is collected, **Then** the file is skipped with a debug log message and the task trace is sent without GPU data.
3. **Given** a completed Fusion-enabled task where `.fusion/trace.json` exists but the `gpu` block is null/absent, **When** the task trace is collected, **Then** no GPU metrics are included and no error is raised.

---

### Edge Cases

- What happens when the `gpu` block contains unexpected or extra fields not in the known schema? They are included as-is (forward compatibility).
- What happens when Fusion is not enabled for a task? No attempt is made to read `.fusion/trace.json`.
- What happens when the task work directory is inaccessible at trace collection time (e.g., remote storage timeout)? The same error handling as existing `.command.trace` parsing applies — log and continue.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST read the file `.fusion/trace.json` from the task work directory on task completion when the executor has Fusion enabled.
- **FR-002**: System MUST extract the entire `gpu` block from the parsed `trace.json` as a map.
- **FR-003**: System MUST store the GPU metrics as a transient field on `TraceRecord` (following the same pattern as `resourceAllocation`).
- **FR-004**: System MUST include the GPU metrics map in the task payload sent to Seqera Platform via the Tower observer.
- **FR-005**: System MUST collect GPU metrics irrespective of task completion status (success or failure).
- **FR-006**: System MUST NOT fail or disrupt task trace reporting if `.fusion/trace.json` is missing, unreadable, or malformed.
- **FR-007**: System MUST only attempt to read `.fusion/trace.json` when Fusion is enabled for the executor.

### Key Entities

- **Fusion Trace File**: JSON file at `.fusion/trace.json` in the task work directory, produced by the Fusion client. Contains `proc`, `gpu`, and `cgroup` blocks with runtime metrics.
- **GPU Metrics Block**: The `gpu` object within `trace.json`, containing fields: `name`, `mem`, `driver`, `active_time`, `pct`, `peak`, `pct_mem`, `peak_mem`, `avg_mem`, `peak_mem_used`, `avg_mem_bw_util`, `peak_mem_bw_util`.
- **TraceRecord GPU field**: New transient field on `TraceRecord` that carries the GPU metrics map through the existing trace pipeline to the Tower observer, following the `resourceAllocation` pattern.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: GPU metrics from Fusion trace files are visible in Seqera Platform for all Fusion-enabled tasks that ran on GPU hardware.
- **SC-002**: Tasks without GPU usage or without Fusion enabled report successfully with no GPU data and no errors.
- **SC-003**: A missing or malformed `.fusion/trace.json` does not cause any task to fail reporting — 100% of tasks still have their standard metrics delivered.
- **SC-004**: GPU metrics collection adds negligible overhead — reading and parsing a single small JSON file per task completion.

## Assumptions

- The Fusion client is responsible for creating `.fusion/trace.json` in the task work directory. Nextflow only reads it.
- The `gpu` block schema may evolve over time. The implementation forwards the entire block as a map rather than mapping to fixed fields, ensuring forward compatibility.
- Seqera Platform API already accepts or will be updated to accept the GPU metrics payload alongside existing task trace data.
- The file path `.fusion/trace.json` is stable and defined by the Fusion client contract.
- All executors that support Fusion (AWS Batch, Google Batch, Azure Batch, Kubernetes, Seqera, SLURM) benefit from this feature without executor-specific code — the detection is based on whether Fusion is enabled, not on the executor type.
