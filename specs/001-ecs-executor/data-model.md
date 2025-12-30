# Data Model: AWS ECS Managed Instances Executor

**Branch**: `001-ecs-executor` | **Date**: 2025-12-30

## Entity Relationship Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              AWS Infrastructure                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────────┐       ┌──────────────────────────┐                    │
│  │   ECS Cluster    │◄──────│  Capacity Provider       │                    │
│  │                  │       │  (Managed Instances)     │                    │
│  │  - clusterArn    │       │                          │                    │
│  │  - clusterName   │       │  - capacityProviderName  │                    │
│  │  - status        │       │  - autoScalingGroupArn   │                    │
│  └────────┬─────────┘       └──────────────────────────┘                    │
│           │                                                                  │
│           │ runs                                                             │
│           ▼                                                                  │
│  ┌──────────────────┐       ┌──────────────────────────┐                    │
│  │   ECS Task       │◄──────│  Task Definition         │                    │
│  │                  │       │                          │                    │
│  │  - taskArn       │       │  - taskDefinitionArn     │                    │
│  │  - lastStatus    │       │  - family                │                    │
│  │  - desiredStatus │       │  - revision              │                    │
│  │  - stoppedReason │       │  - cpu                   │                    │
│  │  - stopCode      │       │  - memory                │                    │
│  │  - containers[]  │       │  - containerDefinitions  │                    │
│  └────────┬─────────┘       │  - executionRoleArn      │                    │
│           │                 │  - taskRoleArn           │                    │
│           │ writes to       │  - storageConfiguration  │                    │
│           ▼                 └──────────────────────────┘                    │
│  ┌──────────────────┐                                                       │
│  │ CloudWatch Logs  │                                                       │
│  │                  │                                                       │
│  │  - logGroupName  │                                                       │
│  │  - logStreamName │                                                       │
│  └──────────────────┘                                                       │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                              Nextflow Runtime                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────────┐       ┌──────────────────────────┐                    │
│  │   TaskRun        │──────►│  AwsEcsTaskHandler       │                    │
│  │                  │       │                          │                    │
│  │  - id            │       │  - taskArn               │                    │
│  │  - name          │       │  - taskDefArn            │                    │
│  │  - container     │       │  - spotAttempts          │                    │
│  │  - config        │       │  - status                │                    │
│  │  - workDir       │       │  - exitCode              │                    │
│  └──────────────────┘       └──────────────────────────┘                    │
│           │                            │                                     │
│           │                            │ uses                                │
│           │                            ▼                                     │
│           │                 ┌──────────────────────────┐                    │
│           │                 │  AwsEcsExecutor          │                    │
│           │                 │                          │                    │
│           │                 │  - ecsClient             │                    │
│           │                 │  - taskDefCache          │                    │
│           │                 │  - options               │                    │
│           │                 └──────────────────────────┘                    │
│           │                            │                                     │
│           │ work directory             │ configured by                       │
│           ▼                            ▼                                     │
│  ┌──────────────────┐       ┌──────────────────────────┐                    │
│  │  S3 Work Dir     │       │  AwsEcsConfig            │                    │
│  │  (via Fusion)    │       │                          │                    │
│  │                  │       │  - cluster               │                    │
│  │  - bucket        │       │  - executionRole         │                    │
│  │  - prefix        │       │  - taskRole              │                    │
│  └──────────────────┘       │  - logsGroup             │                    │
│                             │  - maxSpotAttempts       │                    │
│                             │  - subnets               │                    │
│                             │  - securityGroups        │                    │
│                             └──────────────────────────┘                    │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Entity Definitions

### ECS Cluster

The pre-configured AWS ECS cluster where tasks are executed.

| Attribute | Type | Description |
|-----------|------|-------------|
| clusterArn | String | Full ARN of the ECS cluster |
| clusterName | String | Short name of the cluster |
| status | String | ACTIVE, PROVISIONING, DEPROVISIONING, FAILED, INACTIVE |
| capacityProviders | List<String> | Attached capacity provider names |
| defaultCapacityProviderStrategy | List | Default strategy for tasks |

**Lifecycle**: Pre-created by user, validated at executor startup.

### Capacity Provider

The Managed Instances capacity provider that handles automatic instance provisioning.

| Attribute | Type | Description |
|-----------|------|-------------|
| capacityProviderArn | String | Full ARN of the capacity provider |
| name | String | Capacity provider name (user-defined) |
| autoScalingGroupArn | String | Associated Auto Scaling group |
| status | String | ACTIVE, INACTIVE |
| managedScaling | Object | Scaling configuration |

**Lifecycle**: Pre-created by user as part of cluster setup.

### Task Definition

ECS task definition created from Nextflow process requirements.

| Attribute | Type | Description |
|-----------|------|-------------|
| taskDefinitionArn | String | Full ARN with revision |
| family | String | Task definition family name |
| revision | Integer | Revision number (auto-incremented) |
| cpu | String | CPU units (1024 = 1 vCPU) |
| memory | String | Memory in MiB |
| containerDefinitions | List | Container configurations |
| executionRoleArn | String | IAM role for ECS agent |
| taskRoleArn | String | IAM role for task (S3 access) |
| networkMode | String | Always AWSVPC for Managed Instances |
| requiresCompatibilities | List | Always MANAGED_INSTANCES |
| managedInstancesProvider | Object | Storage configuration |

**Lifecycle**: Auto-created by executor, cached by container+resource hash.

### ECS Task

Running ECS task instance corresponding to a Nextflow task.

| Attribute | Type | Description |
|-----------|------|-------------|
| taskArn | String | Full ARN of the running task |
| clusterArn | String | Cluster where task runs |
| taskDefinitionArn | String | Associated task definition |
| lastStatus | String | Current status (see Status Mapping) |
| desiredStatus | String | Target status |
| containers | List | Container instances with exit codes |
| stoppedReason | String | Reason for task stop |
| stopCode | String | ECS stop code |
| startedAt | Timestamp | Task start time |
| stoppedAt | Timestamp | Task stop time |

**Lifecycle**: Created by RunTask, transitions through states, deleted after polling.

### Container Instance (within Task)

| Attribute | Type | Description |
|-----------|------|-------------|
| containerArn | String | Container instance ARN |
| name | String | Container name (always "main") |
| exitCode | Integer | Process exit code (0-255) |
| lastStatus | String | Container status |
| reason | String | Failure reason if any |

### Work Directory

S3 bucket path used as Nextflow work directory, accessed via Fusion.

| Attribute | Type | Description |
|-----------|------|-------------|
| bucket | String | S3 bucket name |
| prefix | String | Path prefix within bucket |
| uri | String | Full S3 URI (s3://bucket/prefix) |

**Lifecycle**: User-configured, used throughout workflow execution.

## State Diagrams

### ECS Task Status Transitions

```
                    ┌──────────────┐
                    │ PROVISIONING │
                    └──────┬───────┘
                           │
                           ▼
                    ┌──────────────┐
                    │   PENDING    │
                    └──────┬───────┘
                           │
                           ▼
                    ┌──────────────┐
                    │  ACTIVATING  │
                    └──────┬───────┘
                           │
                           ▼
                    ┌──────────────┐
         ┌─────────│   RUNNING    │─────────┐
         │         └──────────────┘         │
         │                                  │
         │ (normal completion)     (cancel/error)
         │                                  │
         ▼                                  ▼
┌──────────────┐                    ┌──────────────┐
│ DEACTIVATING │                    │   STOPPING   │
└──────┬───────┘                    └──────┬───────┘
       │                                   │
       └───────────────┬───────────────────┘
                       │
                       ▼
                ┌──────────────┐
                │   STOPPED    │
                └──────────────┘
```

### Nextflow Task Handler State Transitions

```
┌───────────────┐
│    NEW        │  (TaskHandler created)
└───────┬───────┘
        │ submit()
        ▼
┌───────────────┐
│   SUBMITTED   │  (RunTask called, taskArn obtained)
└───────┬───────┘
        │ checkIfRunning() → ECS status RUNNING
        ▼
┌───────────────┐
│    RUNNING    │  (Task executing)
└───────┬───────┘
        │ checkIfCompleted() → ECS status STOPPED
        │
        ├─────────────────────────────────────────┐
        │ (spot interruption && attempts < max)   │
        │                                         │
        │     resubmit()                          │
        ▼                                         │
┌───────────────┐                                 │
│   SUBMITTED   │◄────────────────────────────────┘
└───────┬───────┘
        │
        │ (normal completion or max attempts)
        ▼
┌───────────────┐
│   COMPLETED   │  (exitCode extracted, status reported)
└───────────────┘
```

## Caching Strategy

### Task Definition Cache

```groovy
// Key structure
String cacheKey = "${containerImage}:${cpuUnits}:${memoryMiB}:${gpuCount}:${diskGiB}"

// Example keys
"ubuntu:latest:2048:4096:0:50"           // 2 vCPU, 4GB RAM, no GPU, 50GB disk
"nvidia/cuda:12:4096:16384:1:100"        // 4 vCPU, 16GB RAM, 1 GPU, 100GB disk
```

### Cache Implementation

```
┌─────────────────────────────────────────────────────────────┐
│                   Task Definition Cache                      │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│   ConcurrentHashMap<String, String>                         │
│                                                              │
│   Key (hash)                    Value (taskDefArn)          │
│   ──────────────────────────    ──────────────────────────  │
│   "ubuntu:2048:4096:0:50"   →   "arn:aws:ecs:...:nf-       │
│                                  ubuntu:1"                   │
│   "nginx:1024:2048:0:30"    →   "arn:aws:ecs:...:nf-       │
│                                  nginx:1"                    │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

## Configuration Data Model

### AwsEcsConfig Fields

| Field | Type | Required | Default | Validation |
|-------|------|----------|---------|------------|
| cluster | String | **Yes** | - | Non-empty, valid cluster name/ARN |
| executionRole | String | **Yes** | - | Valid IAM role ARN |
| taskRole | String | No | - | Valid IAM role ARN if provided |
| subnets | List<String> | No | Auto-discovered | Valid subnet IDs (discovers from default VPC) |
| securityGroups | List<String> | No | Auto-discovered | Valid security group IDs (discovers from default VPC) |
| logsGroup | String | No | `/aws/ecs/nextflow` | Valid log group name |
| maxSpotAttempts | Integer | No | 5 | 1-100 |
| assignPublicIp | Boolean | No | true | - |

### Configuration Hierarchy

**Minimal** (only required settings):
```
aws {
    region = 'us-east-1'              // AwsConfig.region
    ecs {
        cluster = '...'               // AwsEcsConfig.cluster (REQUIRED)
        executionRole = '...'         // AwsEcsConfig.executionRole (REQUIRED)
    }
}
wave { enabled = true }
fusion { enabled = true }
```

**Full** (with all optional settings):
```
aws {
    region = 'us-east-1'              // AwsConfig.region
    ecs {
        cluster = '...'               // AwsEcsConfig.cluster (REQUIRED)
        executionRole = '...'         // AwsEcsConfig.executionRole (REQUIRED)
        taskRole = '...'              // AwsEcsConfig.taskRole (optional)
        subnets = ['...']             // AwsEcsConfig.subnets (auto-discovered)
        securityGroups = ['...']      // AwsEcsConfig.securityGroups (auto-discovered)
        logsGroup = '...'             // AwsEcsConfig.logsGroup (default: /aws/ecs/nextflow)
        maxSpotAttempts = 5           // AwsEcsConfig.maxSpotAttempts (default: 5)
        assignPublicIp = true         // AwsEcsConfig.assignPublicIp (default: true)
    }
}
wave { enabled = true }
fusion { enabled = true }
```
