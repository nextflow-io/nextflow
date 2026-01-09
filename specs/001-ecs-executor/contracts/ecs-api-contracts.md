# ECS API Contracts

**Branch**: `001-ecs-executor` | **Date**: 2025-12-30

This document defines the AWS ECS API contracts used by the `awsecs` executor.

## RegisterTaskDefinition

Creates a new task definition for running Nextflow tasks.

### Request Contract

```java
RegisterTaskDefinitionRequest.builder()
    // Required fields
    .family(String)                           // "nf-{normalized-image-name}"
    .requiresCompatibilities(Compatibility.MANAGED_INSTANCES)
    .networkMode(NetworkMode.AWSVPC)
    .executionRoleArn(String)                 // From aws.ecs.executionRole

    // Resource configuration (task level)
    .cpu(String)                              // CPU units: "1024" = 1 vCPU
    .memory(String)                           // Memory MiB: "4096" = 4 GB

    // Optional task role
    .taskRoleArn(String)                      // From aws.ecs.taskRole

    // Container definition
    .containerDefinitions(
        ContainerDefinition.builder()
            .name("main")                     // Always "main"
            .image(String)                    // Container image from process
            .command("true")                  // Placeholder, overridden at RunTask
            .essential(true)
            .logConfiguration(LogConfiguration.builder()
                .logDriver(LogDriver.AWSLOGS)
                .options(Map.of(
                    "awslogs-group", String,      // From aws.ecs.logsGroup
                    "awslogs-region", String,     // From aws.region
                    "awslogs-stream-prefix", "nf"
                ))
                .build())
            // GPU resources (optional)
            .resourceRequirements(List.of(
                ResourceRequirement.builder()
                    .type(ResourceType.GPU)
                    .value(String)            // Number of GPUs
                    .build()
            ))
            .build()
    )

    // Storage configuration for Managed Instances
    .managedInstancesProvider(
        ManagedInstancesProvider.builder()
            .instanceLaunchTemplate(
                InstanceLaunchTemplate.builder()
                    .storageConfiguration(
                        StorageConfiguration.builder()
                            .storageSizeinGiB(Integer)  // 30-16384 GiB
                            .build()
                    )
                    .build()
            )
            .build()
    )
    .build()
```

### Response Contract

```java
RegisterTaskDefinitionResponse {
    TaskDefinition taskDefinition() {
        String taskDefinitionArn()    // "arn:aws:ecs:region:account:task-definition/family:revision"
        String family()               // Task definition family
        Integer revision()            // Revision number
        String status()               // "ACTIVE" or "INACTIVE"
    }
}
```

### Error Cases

| Error Code | Condition | Handling |
|------------|-----------|----------|
| `ClientException` | Invalid parameters | Fail task with error message |
| `ServerException` | AWS service error | Retry with backoff |
| `InvalidParameterException` | Resource limits exceeded | Fail task with clear message |

---

## RunTask

Submits a task for execution on the ECS cluster.

### Request Contract

```java
RunTaskRequest.builder()
    // Required fields
    .cluster(String)                          // From aws.ecs.cluster
    .taskDefinition(String)                   // ARN from cache/registration

    // Capacity provider strategy (for Managed Instances)
    .capacityProviderStrategy(
        CapacityProviderStrategyItem.builder()
            .capacityProvider(String)         // User-configured capacity provider
            .weight(1)
            .build()
    )

    // Network configuration (required for awsvpc)
    .networkConfiguration(
        NetworkConfiguration.builder()
            .awsvpcConfiguration(
                AwsVpcConfiguration.builder()
                    .subnets(List<String>)        // From aws.ecs.subnets
                    .securityGroups(List<String>) // From aws.ecs.securityGroups
                    .assignPublicIp(AssignPublicIp.ENABLED | DISABLED)
                    .build()
            )
            .build()
    )

    // Task overrides (command, environment)
    .overrides(
        TaskOverride.builder()
            .containerOverrides(
                ContainerOverride.builder()
                    .name("main")
                    .command(List<String>)        // ["bash", "-c", "..."]
                    .environment(List.of(
                        KeyValuePair.builder()
                            .name(String)
                            .value(String)
                            .build()
                    ))
                    .build()
            )
            .build()
    )

    // Optional: Tags for tracking
    .tags(List.of(
        Tag.builder()
            .key("nextflow:taskId")
            .value(String)
            .build()
    ))
    .build()
```

### Response Contract

```java
RunTaskResponse {
    List<Task> tasks() {
        Task {
            String taskArn()              // ARN of started task
            String clusterArn()
            String lastStatus()           // Initially "PROVISIONING" or "PENDING"
            String desiredStatus()        // "RUNNING"
            Instant createdAt()
        }
    }
    List<Failure> failures() {
        Failure {
            String arn()                  // Resource that failed
            String reason()               // Failure reason
            String detail()               // Additional details
        }
    }
}
```

### Error Cases

| Error Code | Condition | Handling |
|------------|-----------|----------|
| `ClusterNotFoundException` | Invalid cluster | Fail with config error |
| `InvalidParameterException` | Bad request parameters | Fail task with message |
| `AccessDeniedException` | Missing IAM permissions | Fail with permissions error |
| `PlatformTaskDefinitionIncompatibilityException` | Wrong launch type | Fail with config error |
| `BlockedException` | Account blocked | Fail with clear message |

---

## DescribeTasks

Polls task status for monitoring execution.

### Request Contract

```java
DescribeTasksRequest.builder()
    .cluster(String)                          // From aws.ecs.cluster
    .tasks(List<String>)                      // Task ARNs (up to 100)
    .include(TaskField.TAGS)                  // Include tags in response
    .build()
```

### Response Contract

```java
DescribeTasksResponse {
    List<Task> tasks() {
        Task {
            String taskArn()
            String taskDefinitionArn()
            String clusterArn()
            String lastStatus()               // PROVISIONING, PENDING, ACTIVATING,
                                              // RUNNING, DEACTIVATING, STOPPING, STOPPED
            String desiredStatus()            // RUNNING or STOPPED
            TaskStopCode stopCode()           // TASK_FAILED_TO_START, SPOT_INTERRUPTION,
                                              // ESSENTIAL_CONTAINER_EXITED, USER_INITIATED, etc.
            String stoppedReason()            // Human-readable reason
            Instant startedAt()
            Instant stoppedAt()
            List<Container> containers() {
                Container {
                    String name()             // "main"
                    String containerArn()
                    Integer exitCode()        // 0-255, null if still running
                    String lastStatus()
                    String reason()           // Failure reason if any
                }
            }
        }
    }
    List<Failure> failures() {
        Failure {
            String arn()
            String reason()                   // Usually "MISSING" for deleted tasks
        }
    }
}
```

### Status Mapping

| ECS lastStatus | Nextflow TaskStatus | Action |
|----------------|---------------------|--------|
| PROVISIONING | SUBMITTED | Continue polling |
| PENDING | SUBMITTED | Continue polling |
| ACTIVATING | SUBMITTED | Continue polling |
| RUNNING | RUNNING | Continue polling |
| DEACTIVATING | RUNNING | Continue polling |
| STOPPING | RUNNING | Continue polling |
| STOPPED | COMPLETED | Extract exit code, complete |

### Exit Code Extraction Logic

```groovy
Integer extractExitCode(Task task) {
    def container = task.containers().find { it.name() == "main" }
    def exitCode = container?.exitCode()

    // Handle special cases
    if (exitCode == null) {
        def stopCode = task.stopCode()
        if (stopCode == TaskStopCode.TASK_FAILED_TO_START) {
            return 1  // Failed to start
        }
        if (stopCode == TaskStopCode.ESSENTIAL_CONTAINER_EXITED) {
            return 1  // Container crashed
        }
        // Unknown failure
        return 1
    }

    return exitCode
}
```

---

## StopTask

Cancels a running task.

### Request Contract

```java
StopTaskRequest.builder()
    .cluster(String)                          // From aws.ecs.cluster
    .task(String)                             // Task ARN
    .reason("Cancelled by Nextflow")          // Stop reason
    .build()
```

### Response Contract

```java
StopTaskResponse {
    Task task() {
        String taskArn()
        String lastStatus()                   // May still be RUNNING briefly
        String desiredStatus()                // STOPPED
    }
}
```

### Error Cases

| Error Code | Condition | Handling |
|------------|-----------|----------|
| `InvalidParameterException` | Task not found | Log warning, ignore (already stopped) |
| `ClusterNotFoundException` | Invalid cluster | Log error |

---

## DescribeClusters

Validates cluster configuration at startup.

### Request Contract

```java
DescribeClustersRequest.builder()
    .clusters(String)                         // Cluster name or ARN
    .include(ClusterField.ATTACHMENTS,
             ClusterField.SETTINGS,
             ClusterField.STATISTICS,
             ClusterField.TAGS)
    .build()
```

### Response Contract

```java
DescribeClustersResponse {
    List<Cluster> clusters() {
        Cluster {
            String clusterArn()
            String clusterName()
            String status()                   // ACTIVE, PROVISIONING, DEPROVISIONING,
                                              // FAILED, INACTIVE
            List<String> capacityProviders()  // Attached capacity providers
            Integer registeredContainerInstancesCount()
            Integer runningTasksCount()
            Integer pendingTasksCount()
        }
    }
    List<Failure> failures() {
        Failure {
            String arn()
            String reason()                   // "MISSING" if not found
        }
    }
}
```

### Validation Logic

```groovy
void validateCluster(String clusterName) {
    def response = ecsClient.describeClusters(
        DescribeClustersRequest.builder()
            .clusters(clusterName)
            .build()
    )

    if (response.failures()) {
        throw new ProcessUnrecoverableException(
            "ECS cluster not found: ${clusterName}")
    }

    def cluster = response.clusters()[0]

    if (cluster.status() != "ACTIVE") {
        throw new ProcessUnrecoverableException(
            "ECS cluster is not active: ${cluster.status()}")
    }

    // Check for Managed Instances capacity provider
    def hasManaged = cluster.capacityProviders().any {
        // User must configure capacity provider name
        it == options.capacityProvider
    }

    if (!hasManaged) {
        throw new ProcessUnrecoverableException(
            "ECS cluster lacks configured Managed Instances capacity provider")
    }
}
```

---

## CloudWatch Logs Integration

### GetLogEvents Request

```java
GetLogEventsRequest.builder()
    .logGroupName(String)                     // From aws.ecs.logsGroup
    .logStreamName(String)                    // "nf/main/{task-id}"
    .startFromHead(true)
    .build()
```

### Log Stream Name Format

```
{prefix}/{container-name}/{task-id}
nf/main/1a2b3c4d-5e6f-7g8h-9i0j-1k2l3m4n5o6p
```

Where `task-id` is extracted from the task ARN:
```
arn:aws:ecs:us-east-1:123456789:task/cluster-name/1a2b3c4d-5e6f-7g8h-9i0j-1k2l3m4n5o6p
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                   task-id
```

---

## Rate Limiting

ECS API has the following rate limits (requests per second):

| API Operation | Burst | Sustained |
|---------------|-------|-----------|
| RunTask | 100 | 20 |
| DescribeTasks | 100 | 40 |
| StopTask | 100 | 20 |
| RegisterTaskDefinition | 100 | 1 |
| DescribeClusters | 100 | 20 |

### Throttling Implementation

The `AwsEcsProxy` class wraps the ECS client with a `ThrottlingExecutor`:

```groovy
class AwsEcsProxy {
    private final EcsClient client
    private final ThrottlingExecutor throttle

    RunTaskResponse runTask(RunTaskRequest request) {
        throttle.execute { client.runTask(request) }
    }

    DescribeTasksResponse describeTasks(DescribeTasksRequest request) {
        throttle.execute { client.describeTasks(request) }
    }

    // ... other methods
}
```
