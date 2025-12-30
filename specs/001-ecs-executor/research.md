# Research: AWS ECS Managed Instances Executor

**Date**: 2025-12-30 | **Branch**: `001-ecs-executor`

## 1. ECS Managed Instances API Patterns

### Decision: Use ECS RunTask API with Managed Instances Capacity Provider

**Rationale**: ECS Managed Instances uses `RunTask` with capacity provider strategy pointing to a Managed Instances capacity provider. This is different from AWS Batch which abstracts away instance management.

**Alternatives Considered**:
- ECS Services: Rejected - Services are for long-running tasks, not batch jobs
- Direct EC2 with ECS agent: Rejected - Managed Instances handles this automatically

### Key API Operations

| Operation | API | Notes |
|-----------|-----|-------|
| Register task definition | `RegisterTaskDefinition` | One-time per container/resource combo |
| Run task | `RunTask` | With `capacityProviderStrategy` |
| Check status | `DescribeTasks` | Poll task ARN for status |
| Stop task | `StopTask` | For cancellation |
| Describe cluster | `DescribeClusters` | For validation at startup |

### RunTask Request Structure

```java
RunTaskRequest.builder()
    .cluster(clusterArn)
    .taskDefinition(taskDefArn)
    .capacityProviderStrategy(
        CapacityProviderStrategyItem.builder()
            .capacityProvider("managed-instances-cp")  // User-configured
            .weight(1)
            .build()
    )
    .networkConfiguration(
        NetworkConfiguration.builder()
            .awsvpcConfiguration(
                AwsVpcConfiguration.builder()
                    .subnets(subnets)
                    .securityGroups(securityGroups)
                    .assignPublicIp(AssignPublicIp.ENABLED)
                    .build()
            )
            .build()
    )
    .overrides(
        TaskOverride.builder()
            .containerOverrides(
                ContainerOverride.builder()
                    .name("main")
                    .command(command)
                    .environment(envVars)
                    .build()
            )
            .build()
    )
    .build()
```

## 2. Task Definition Management

### Decision: Cache task definitions by container+resource hash

**Rationale**: ECS task definitions are versioned (family:revision). Creating a new revision for each task would be inefficient. Cache by hash of image + CPU + memory + GPU + disk configuration.

**Alternatives Considered**:
- Single task definition with overrides: Rejected - Resource requirements (CPU/memory at task level) cannot be overridden
- No caching: Rejected - Would create excessive task definition revisions

### Task Definition Structure

```java
RegisterTaskDefinitionRequest.builder()
    .family("nf-" + normalizedImageName)
    .requiresCompatibilities(Compatibility.MANAGED_INSTANCES)
    .networkMode(NetworkMode.AWSVPC)
    .cpu(cpuUnits)           // Task-level CPU (cannot override)
    .memory(memoryMiB)       // Task-level memory (cannot override)
    .executionRoleArn(executionRole)
    .taskRoleArn(taskRole)
    .containerDefinitions(
        ContainerDefinition.builder()
            .name("main")
            .image(containerImage)
            .command("true")  // Overridden at RunTask
            .essential(true)
            .logConfiguration(cloudWatchLogs)
            .resourceRequirements(gpuRequirements)
            .build()
    )
    .managedInstancesProvider(
        ManagedInstancesProvider.builder()
            .instanceLaunchTemplate(
                InstanceLaunchTemplate.builder()
                    .storageConfiguration(
                        StorageConfiguration.builder()
                            .storageSizeinGiB(diskSizeGiB)
                            .build()
                    )
                    .build()
            )
            .build()
    )
    .build()
```

## 3. Resource Mapping

### Decision: Map Nextflow directives to ECS task definition fields

**CPU Mapping**:
- Nextflow `cpus = 2` → ECS `cpu = "2048"` (1 vCPU = 1024 units)
- ECS Managed Instances supports flexible CPU allocation

**Memory Mapping**:
- Nextflow `memory = '4 GB'` → ECS `memory = "4096"` (MiB)
- Direct conversion from MemoryUnit

**GPU Mapping**:
- Nextflow `accelerator 1, type: 'nvidia-tesla-t4'`
- ECS uses `resourceRequirements` with `GPU` type
- Triggers attribute-based instance selection for GPU instances

**Disk Mapping**:
- Nextflow `disk = '500 GB'`
- ECS `storageConfiguration.storageSizeinGiB = 500`
- Range: 30 GiB - 16,384 GiB

## 4. Task Status Polling

### Decision: Use BatchContext pattern for efficient status polling

**Rationale**: Following AWS Batch executor pattern, aggregate multiple task status queries into single `DescribeTasks` call (up to 100 tasks per call).

**Task Status Mapping**:

| ECS Task Status | Nextflow TaskStatus | Action |
|-----------------|---------------------|--------|
| PROVISIONING | SUBMITTED | Continue polling |
| PENDING | SUBMITTED | Continue polling |
| ACTIVATING | SUBMITTED | Continue polling |
| RUNNING | RUNNING | Continue polling |
| DEACTIVATING | RUNNING | Continue polling |
| STOPPING | RUNNING | Continue polling |
| STOPPED | COMPLETED | Extract exit code |

**Exit Code Extraction**:
```groovy
def exitCode = task.containers()[0].exitCode()
def stoppedReason = task.stoppedReason()
def stopCode = task.stopCode()

// Handle special cases
if (stopCode == TaskStopCode.SPOT_INTERRUPTION) {
    // Trigger spot retry if configured
}
if (exitCode == null && stopCode == TaskStopCode.ESSENTIAL_CONTAINER_EXITED) {
    exitCode = 1  // Default to failure
}
```

## 5. Spot Instance Retry

### Decision: Implement spot retry similar to aws.batch.maxSpotAttempts

**Rationale**: ECS Managed Instances can use spot capacity. When spot interruption occurs, the task should automatically retry up to a configurable limit.

**Detection**:
```groovy
boolean isSpotInterruption(Task task) {
    return task.stopCode() == TaskStopCode.SPOT_INTERRUPTION ||
           task.stoppedReason()?.contains('spot') ||
           task.stoppedReason()?.contains('Host EC2')
}
```

**Retry Logic**:
- Track `spotAttempts` counter per task handler
- On spot interruption, increment counter and resubmit if under limit
- Log warning on each retry
- Fail task if max attempts exceeded

## 6. CloudWatch Logs Integration

### Decision: Use awslogs driver with configurable log group

**Rationale**: Standard ECS logging pattern. Log group configurable via `aws.ecs.logsGroup`.

**Log Configuration**:
```java
LogConfiguration.builder()
    .logDriver(LogDriver.AWSLOGS)
    .options(Map.of(
        "awslogs-group", logsGroup,
        "awslogs-region", region,
        "awslogs-stream-prefix", "nf"
    ))
    .build()
```

**Log Retrieval**:
- Log stream name format: `nf/<container-name>/<task-id>`
- Use CloudWatch Logs `GetLogEvents` API
- Implement in AwsEcsHelper class

## 7. Cluster Validation

### Decision: Validate cluster and capacity provider at startup

**Rationale**: Fail fast if cluster doesn't exist or lacks Managed Instances capacity provider.

**Validation Steps**:
1. Call `DescribeClusters` with cluster name/ARN
2. Check cluster status is ACTIVE
3. Verify capacity providers include a Managed Instances provider
4. If validation fails, throw `ProcessUnrecoverableException` with clear message

## 8. Fusion Integration

### Decision: Reuse existing AwsFusionEnv implementation

**Rationale**: Fusion environment setup for S3 is identical between Batch and ECS. Both need:
- AWS credentials in environment
- S3 endpoint configuration
- Fusion-specific environment variables

**Implementation**:
- Implement `FusionAwareTask` trait in AwsEcsTaskHandler
- Use `FusionScriptLauncher` when Fusion enabled
- Pass Fusion environment variables to container

## 9. Configuration Namespace

### Decision: Use `aws.ecs.*` namespace with minimal required settings

**Rationale**: Minimize user configuration burden by requiring only cluster and executionRole. Network settings (subnets, security groups) are auto-discovered from the default VPC when not specified.

**Required Configuration** (only 2 settings):

| Option | Type | Description |
|--------|------|-------------|
| `cluster` | String | ECS cluster name or ARN with Managed Instances capacity provider |
| `executionRole` | String | Task execution IAM role ARN (for image pull, CloudWatch logs) |

**Optional Configuration** (with defaults/auto-discovery):

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `taskRole` | String | - | Task IAM role ARN (for S3 access from containers) |
| `subnets` | List | Auto-discovered | VPC subnets (discovers from default VPC) |
| `securityGroups` | List | Auto-discovered | Security groups (discovers from default VPC) |
| `logsGroup` | String | `/aws/ecs/nextflow` | CloudWatch log group |
| `maxSpotAttempts` | Integer | 5 | Max spot retry attempts |
| `assignPublicIp` | Boolean | true | Assign public IP (enables internet without NAT) |

### Network Auto-Discovery

When `subnets` or `securityGroups` are not specified, the executor auto-discovers them from the default VPC using EC2 API:

```groovy
// Auto-discover default VPC subnets
List<String> discoverDefaultSubnets() {
    def response = ec2Client.describeSubnets(
        DescribeSubnetsRequest.builder()
            .filters(
                Filter.builder()
                    .name("default-for-az")
                    .values("true")
                    .build()
            )
            .build()
    )
    return response.subnets().collect { it.subnetId() }
}

// Auto-discover default security group
List<String> discoverDefaultSecurityGroup() {
    def response = ec2Client.describeSecurityGroups(
        DescribeSecurityGroupsRequest.builder()
            .filters(
                Filter.builder()
                    .name("group-name")
                    .values("default")
                    .build()
            )
            .build()
    )
    return response.securityGroups().collect { it.groupId() }
}
```

## 10. Extension Point Registration

### Decision: Use @ServiceName annotation for auto-discovery

**Rationale**: Following existing executor pattern, the executor class annotated with `@ServiceName('awsecs')` and implementing `ExtensionPoint` is automatically discovered by Nextflow's plugin system.

**Registration**:
```groovy
@ServiceName('awsecs')
@CompileStatic
class AwsEcsExecutor extends Executor implements ExtensionPoint {
    // ...
}
```

**Build Configuration** (build.gradle):
```groovy
extensionPoints = [
    'nextflow.cloud.aws.batch.AwsBatchExecutor',
    'nextflow.cloud.aws.ecs.AwsEcsExecutor',  // ADD
    // ...
]
```

## Summary

All technical decisions are aligned with:
- Existing AWS Batch executor patterns in nf-amazon plugin
- AWS ECS Managed Instances best practices
- Nextflow executor architecture requirements

No clarifications remain outstanding.
