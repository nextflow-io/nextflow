# Implementation Plan: AWS ECS Managed Instances Executor

**Branch**: `001-ecs-executor` | **Date**: 2025-12-30 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/001-ecs-executor/spec.md`

## Summary

Implement a new `awsecs` executor in the nf-amazon plugin that runs Nextflow tasks on AWS ECS using the Managed Instances capacity provider. The executor will support configurable CPU, memory, GPU, and disk resources, use S3/Fusion for storage, and follow existing patterns from the AWS Batch executor while remaining completely independent.

## Technical Context

**Language/Version**: Groovy 4.0.29, Java 17 target (Java 21 toolchain)
**Primary Dependencies**: AWS SDK v2 (ECS, EC2, CloudWatch Logs, S3), GPars 1.2.1
**Storage**: AWS S3 via Seqera Fusion filesystem
**Testing**: Spock Framework, JaCoCo for coverage
**Target Platform**: Nextflow runtime on any platform with AWS credentials
**Project Type**: Plugin module within nf-amazon
**Performance Goals**: Task submission <60s, status polling efficient via batch aggregation
**Constraints**: ECS Managed Instances 14-day lifecycle limit, Wave + Fusion required for S3 access
**Scale/Scope**: Support equivalent scale to AWS Batch executor

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Modular Architecture | ✅ PASS | Plugin feature in `plugins/nf-amazon`, independent from Batch executor |
| II. Test-Driven Quality | ✅ PASS | Unit tests (Spock), integration tests, smoke tests planned |
| III. Dataflow Programming Model | ✅ PASS | Executor is infrastructure layer, doesn't affect dataflow semantics |
| IV. Apache 2.0 License | ✅ PASS | All new files will include Apache 2.0 headers |
| V. DCO Sign-off | ✅ PASS | All commits will use `-s` flag |
| VI. Semantic Versioning | ✅ PASS | Plugin version bump required in `plugins/nf-amazon/VERSION` |
| VII. Groovy Idioms | ✅ PASS | Follow existing Batch executor patterns and conventions |

**Gate Status**: PASSED - No violations requiring justification.

## Project Structure

### Documentation (this feature)

```text
specs/001-ecs-executor/
├── spec.md              # Feature specification (complete)
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
├── contracts/           # Phase 1 output (ECS API patterns)
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code (repository root)

```text
plugins/nf-amazon/src/main/nextflow/cloud/aws/
├── ecs/                           # NEW - ECS executor package
│   ├── AwsEcsExecutor.groovy      # Main executor class
│   ├── AwsEcsTaskHandler.groovy   # Task lifecycle management
│   ├── AwsEcsHelper.groovy        # CloudWatch logs, metadata
│   ├── AwsEcsOptions.groovy       # Configuration wrapper
│   └── model/                     # ECS-specific model classes
│       ├── RegisterTaskDefinitionModel.groovy
│       └── ContainerDefinitionModel.groovy
├── config/
│   ├── AwsConfig.groovy           # MODIFY - Add ecs config field
│   └── AwsEcsConfig.groovy        # NEW - ECS configuration
└── AmazonPlugin.groovy            # VERIFY - Extension point may auto-register

# DEFERRED (optimization phase):
# ├── ecs/AwsEcsProxy.groovy       # Throttled ECS client wrapper

plugins/nf-amazon/src/test/nextflow/cloud/aws/
├── ecs/                           # NEW - ECS executor tests
│   ├── AwsEcsExecutorTest.groovy
│   ├── AwsEcsTaskHandlerTest.groovy
│   └── AwsEcsConfigTest.groovy
```

**Structure Decision**: Follow the existing `batch/` package structure exactly, creating a parallel `ecs/` package. This maintains consistency with the modular architecture and keeps Batch and ECS code independent.

## Architecture Design

### Component Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                         Nextflow Runtime                        │
├─────────────────────────────────────────────────────────────────┤
│  ┌──────────────────┐    ┌──────────────────┐                   │
│  │   TaskProcessor  │────│  TaskMonitor     │                   │
│  └────────┬─────────┘    │  (Polling Loop)  │                   │
│           │              └────────┬─────────┘                   │
├───────────┼───────────────────────┼─-───────────────────────────┤
│           ▼                       ▼                             │
│  ┌──────────────────┐    ┌──────────────────┐                   │
│  │  AwsEcsExecutor  │────│ AwsEcsTaskHandler│                   │
│  │  (@ServiceName   │    │  (per task)      │                   │
│  │   'awsecs')      │    └────────┬─────────┘                   │
│  └────────┬─────────┘             │                             │
│           │                       │                             │
│  ┌────────┴─────────┐    ┌────────┴─────────┐                   │
│  │  AwsEcsProxy     │    │  BatchContext    │                   │
│  │  (throttling)    │    │  (aggregation)   │                   │
│  └────────┬─────────┘    └────────┬─────────┘                   │
├───────────┼───────────────────────┼─-───────────────────────────┤
│           ▼                       ▼                             │
│  ┌─────────────────────────────────────────────────────────────┐│
│  │                    AWS ECS API                              ││
│  │  RunTask │ DescribeTasks │ StopTask │ RegisterTaskDefinition││
│  └─────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────┘
```

### Key Classes

| Class | Responsibility | Key Methods |
|-------|---------------|-------------|
| `AwsEcsExecutor` | Executor lifecycle, client init | `register()`, `createTaskHandler()`, `createTaskMonitor()` |
| `AwsEcsTaskHandler` | Task submission, polling, cancel | `submit()`, `checkIfRunning()`, `checkIfCompleted()`, `killTask()` |
| `AwsEcsConfig` | Configuration validation | `@ConfigOption` annotated fields |
| `AwsEcsHelper` | CloudWatch logs, metadata | `getTaskLogStream()`, `describeCluster()` |
| `AwsEcsProxy` | Request throttling (DEFERRED) | Wraps `EcsClient` with `ThrottlingExecutor` |

### ECS API Mapping

| Nextflow Operation | AWS Batch | AWS ECS |
|--------------------|-----------|---------|
| Create definition | `RegisterJobDefinition` | `RegisterTaskDefinition` |
| Submit task | `SubmitJob` | `RunTask` |
| Query status | `DescribeJobs` | `DescribeTasks` |
| Cancel task | `TerminateJob` | `StopTask` |
| Get logs | CloudWatch Logs | CloudWatch Logs |

### Configuration Schema

**Minimal configuration** (only 2 required settings):

```groovy
aws {
    region = 'us-east-1'
    ecs {
        cluster = 'my-nextflow-cluster'           // REQUIRED
        executionRole = 'arn:aws:iam::...'        // REQUIRED
    }
}

wave {
    enabled = true                                // Required for Fusion
}

fusion {
    enabled = true
}
```

**Full configuration** (with optional settings):

```groovy
aws {
    region = 'us-east-1'
    ecs {
        cluster = 'my-nextflow-cluster'           // REQUIRED
        executionRole = 'arn:aws:iam::...'        // REQUIRED
        taskRole = 'arn:aws:iam::...'             // Optional: for S3 access
        subnets = ['subnet-abc', 'subnet-def']    // Optional: auto-discovered from default VPC
        securityGroups = ['sg-xyz']               // Optional: auto-discovered from default VPC
        logsGroup = '/aws/ecs/nextflow'           // Default: /aws/ecs/nextflow
        maxSpotAttempts = 5                       // Default: 5
        assignPublicIp = true                     // Default: true
    }
}

wave {
    enabled = true
}

fusion {
    enabled = true
}
```

### Resource Mapping

| Nextflow Directive | ECS Task Definition Field | Notes |
|--------------------|--------------------------|-------|
| `cpus` | `cpu` (CPU units) | 1 vCPU = 1024 units |
| `memory` | `memory` (MiB) | Direct mapping |
| `accelerator` | `resourceRequirements[GPU]` | Triggers GPU instance selection |
| `disk` | `storageConfiguration.storageSizeinGiB` | 30-16384 GiB |
| `container` | `containerDefinitions[0].image` | Docker image |
| `machineType` | Capacity provider attributes | Instance type constraint |

### Task Definition Caching Strategy

```groovy
// Global cache keyed by container:hash
static final Map<String, String> taskDefinitions = [:]

String resolveTaskDefinition(TaskRun task) {
    def key = "${task.container}:${computeHash(task)}"

    // Check cache
    if (taskDefinitions.containsKey(key))
        return taskDefinitions[key]

    synchronized(taskDefinitions) {
        // Double-check after lock
        if (taskDefinitions.containsKey(key))
            return taskDefinitions[key]

        // Find existing or create new
        def taskDefArn = findOrCreateTaskDefinition(task)
        taskDefinitions[key] = taskDefArn
        return taskDefArn
    }
}
```

### Spot Retry Strategy

```groovy
// Detect spot interruption from task stop reason
boolean isSpotInterruption(Task task) {
    def reason = task.stoppedReason()
    return reason?.contains('Host EC2') ||
           reason?.contains('spot interruption') ||
           task.stopCode() == TaskStopCode.SPOT_INTERRUPTION
}

// Retry logic in checkIfCompleted()
if (isSpotInterruption(task) && spotAttempts < maxSpotAttempts) {
    spotAttempts++
    log.debug "Spot interruption detected, retry $spotAttempts/$maxSpotAttempts"
    resubmit()
    return false  // Not completed yet
}
```

## Complexity Tracking

> No constitution violations requiring justification.

## Implementation Phases

### Phase 1: Core Infrastructure (P1 - MVP)
- AwsEcsExecutor skeleton with lifecycle methods
- AwsEcsConfig with required options (cluster, executionRole)
- AwsEcsTaskHandler with basic submit/poll/cancel
- Task definition registration (no caching initially)
- Basic Fusion integration

### Phase 2: Resource Management (P2)
- CPU/memory mapping with validation
- GPU support via attribute-based instance selection
- Disk size configuration
- Instance type specification (optional)
- Task definition caching

### Phase 3: Reliability & Observability (P2/P3)
- Spot retry logic (maxSpotAttempts)
- CloudWatch Logs integration
- Cluster validation at startup
- BatchContext for efficient status polling

### Phase 4: Optimization (Deferred)
- Rate limiting via AwsEcsProxy (throttled ECS client wrapper)

### Phase 5: Documentation & Testing
- Unit tests for all components
- Integration tests with mock ECS
- Smoke tests
- User documentation for cluster setup
- Configuration reference documentation

## Dependencies

### Internal (nf-amazon plugin)
- `AwsClientFactory` - Reuse for ECS client creation
- `AwsConfig` - Extend with `ecs` field
- `AwsFusionEnv` - Reuse for S3 credentials
- `ThrottlingExecutor` - Reuse for rate limiting

### External (AWS SDK)
- `software.amazon.awssdk:ecs:2.33.2` - Already included in plugin
- `software.amazon.awssdk:cloudwatchlogs:2.33.2` - Already included

### Core Nextflow
- `BatchHandler` trait - Implement for batch aggregation
- `FusionAwareTask` trait - Implement for Fusion support
- `TaskHandler` - Extend for task lifecycle

## Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| ECS API rate limits | Task submission delays | ThrottlingExecutor with configurable rate |
| 14-day lifecycle limit | Long tasks fail | Clear error message, documentation |
| Managed Instances capacity | Task scheduling delays | Document cluster setup, capacity provider config |
| Spot interruptions | Task failures | Automatic retry with configurable max attempts |
