# Quickstart: AWS ECS Managed Instances Executor

**Branch**: `001-ecs-executor` | **Date**: 2025-12-30

This guide helps you get started with the `awsecs` executor to run Nextflow workflows on AWS ECS Managed Instances.

## Prerequisites

Before using the `awsecs` executor, ensure you have:

1. **AWS Account** with appropriate permissions for ECS, EC2, S3, CloudWatch, and IAM
2. **Nextflow** installed with the nf-amazon plugin
3. **Wave containers** enabled (required for Fusion)
4. **Fusion** filesystem configured (required for S3 work directory)
5. **ECS Cluster** with Managed Instances capacity provider (see Setup section)

## Quick Example

### Minimal Configuration

The executor requires only two AWS settings. Network configuration (subnets, security groups) is auto-discovered from the default VPC.

```groovy
// nextflow.config
process {
    executor = 'awsecs'
}

aws {
    region = 'us-east-1'
    ecs {
        cluster = 'nextflow-cluster'                                        // REQUIRED
        executionRole = 'arn:aws:iam::123456789:role/ecsTaskExecutionRole'  // REQUIRED
    }
}

wave {
    enabled = true
}

fusion {
    enabled = true
}
```

### Full Configuration (with explicit networking)

For production deployments or custom VPC setups, you can explicitly specify network configuration:

```groovy
// nextflow.config
process {
    executor = 'awsecs'
}

aws {
    region = 'us-east-1'
    ecs {
        cluster = 'nextflow-cluster'
        executionRole = 'arn:aws:iam::123456789:role/ecsTaskExecutionRole'
        taskRole = 'arn:aws:iam::123456789:role/ecsTaskRole'  // Optional: for S3 access
        subnets = ['subnet-abc123', 'subnet-def456']          // Optional: auto-discovered
        securityGroups = ['sg-xyz789']                         // Optional: auto-discovered
        logsGroup = '/aws/ecs/nextflow'                        // Optional: default shown
        maxSpotAttempts = 5                                    // Optional: default shown
        assignPublicIp = true                                  // Optional: default shown
    }
}

wave {
    enabled = true
}

fusion {
    enabled = true
}
```

### Sample Workflow

```groovy
// main.nf
process HELLO {
    container 'ubuntu:latest'
    cpus 2
    memory '4 GB'

    output:
    stdout

    script:
    '''
    echo "Hello from ECS Managed Instances!"
    hostname
    '''
}

workflow {
    HELLO()
}
```

### Run

```bash
nextflow run main.nf -work-dir s3://my-bucket/work
```

## AWS Infrastructure Setup

### 1. Create ECS Cluster with Managed Instances

```bash
# Create cluster
aws ecs create-cluster \
    --cluster-name nextflow-cluster \
    --capacity-providers MANAGED_INSTANCES \
    --default-capacity-provider-strategy \
        capacityProvider=MANAGED_INSTANCES,weight=1

# Or via CloudFormation/Terraform (recommended for production)
```

### 2. Create IAM Roles

**Task Execution Role** (required - allows ECS to pull images and write logs):

```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "ecr:GetAuthorizationToken",
                "ecr:BatchCheckLayerAvailability",
                "ecr:GetDownloadUrlForLayer",
                "ecr:BatchGetImage",
                "logs:CreateLogStream",
                "logs:PutLogEvents"
            ],
            "Resource": "*"
        }
    ]
}
```

Trust policy:
```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Principal": {
                "Service": "ecs-tasks.amazonaws.com"
            },
            "Action": "sts:AssumeRole"
        }
    ]
}
```

**Task Role** (optional - for S3 access from within tasks):

```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:DeleteObject",
                "s3:ListBucket"
            ],
            "Resource": [
                "arn:aws:s3:::my-bucket",
                "arn:aws:s3:::my-bucket/*"
            ]
        }
    ]
}
```

### 3. Configure VPC Networking

Ensure your subnets:
- Have internet access (NAT gateway for private subnets, or use public subnets)
- Allow outbound traffic to AWS services (ECR, S3, CloudWatch)
- Allow any inter-task communication if needed

Security group should allow:
- Outbound HTTPS (443) for AWS API calls
- Outbound for any external services your workflows need

### 4. Create CloudWatch Log Group

```bash
aws logs create-log-group --log-group-name /aws/ecs/nextflow
```

## Configuration Reference

### Required Options

| Option | Description |
|--------|-------------|
| `aws.ecs.cluster` | ECS cluster name or ARN with Managed Instances capacity provider |
| `aws.ecs.executionRole` | Task execution IAM role ARN (for image pull, CloudWatch logs) |

### Optional Options

| Option | Default | Description |
|--------|---------|-------------|
| `aws.ecs.taskRole` | - | IAM role for task (S3 access from within containers) |
| `aws.ecs.subnets` | Auto-discovered | VPC subnet IDs (uses default VPC subnets if not specified) |
| `aws.ecs.securityGroups` | Auto-discovered | Security group IDs (uses default VPC security group if not specified) |
| `aws.ecs.logsGroup` | `/aws/ecs/nextflow` | CloudWatch log group |
| `aws.ecs.maxSpotAttempts` | `5` | Max retries for spot interruptions |
| `aws.ecs.assignPublicIp` | `true` | Assign public IP to tasks (enables internet without NAT gateway) |

### Process Directives

| Directive | ECS Mapping | Notes |
|-----------|-------------|-------|
| `cpus` | Task CPU units | 1 vCPU = 1024 units |
| `memory` | Task memory (MiB) | Direct mapping |
| `disk` | Storage size (GiB) | 30-16384 GiB |
| `accelerator` | GPU resources | Triggers GPU instance selection |
| `container` | Container image | Required |
| `machineType` | Instance type | Optional constraint |

## Usage Examples

### Basic CPU/Memory Task

```groovy
process CPU_TASK {
    container 'ubuntu:latest'
    cpus 4
    memory '8 GB'

    script:
    '''
    # Your computation here
    '''
}
```

### GPU-Accelerated Task

```groovy
process GPU_TRAINING {
    container 'nvidia/cuda:12.0-runtime-ubuntu22.04'
    cpus 4
    memory '16 GB'
    accelerator 1, type: 'nvidia-tesla-t4'

    script:
    '''
    nvidia-smi
    python train_model.py
    '''
}
```

### Large Storage Task

```groovy
process BIG_DATA {
    container 'ubuntu:latest'
    cpus 2
    memory '4 GB'
    disk '500 GB'

    script:
    '''
    # Process large files
    '''
}
```

### Specific Instance Type

```groovy
process MEMORY_INTENSIVE {
    container 'ubuntu:latest'
    cpus 4
    memory '64 GB'
    machineType 'r6i.2xlarge'

    script:
    '''
    # Memory-intensive computation
    '''
}
```

## Monitoring and Debugging

### View Task Logs

Logs are written to CloudWatch:

```bash
# Via AWS CLI
aws logs tail /aws/ecs/nextflow --follow

# Or in AWS Console: CloudWatch > Log groups > /aws/ecs/nextflow
```

### Check Task Status

```bash
# List running tasks
aws ecs list-tasks --cluster nextflow-cluster

# Describe specific task
aws ecs describe-tasks \
    --cluster nextflow-cluster \
    --tasks <task-arn>
```

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Task stuck in PROVISIONING | No capacity | Check capacity provider configuration |
| Task fails to start | Image pull error | Verify ECR permissions, image exists |
| Network timeout | No internet access | Check NAT gateway, security groups |
| Permission denied (S3) | Missing task role | Add S3 permissions to task role |

## Limitations

1. **14-Day Lifecycle**: ECS Managed Instances tasks have a 14-day maximum runtime. Tasks exceeding this limit will fail.

2. **No Task Arrays**: Unlike AWS Batch, ECS doesn't support native task arrays. Each task is submitted individually.

3. **Wave + Fusion Required**: The executor requires Wave containers and Fusion filesystem for S3 access. Local work directories are not supported.

## Comparison with AWS Batch

| Feature | awsecs | awsbatch |
|---------|--------|----------|
| Infrastructure | ECS Managed Instances | AWS Batch compute environments |
| Task arrays | Not supported | Supported |
| Max runtime | 14 days | Unlimited |
| GPU support | Yes (via EC2) | Yes |
| Spot support | Yes (via capacity provider) | Yes |
| Pricing | ECS + EC2 | Batch + EC2 |

Choose `awsecs` when:
- You prefer direct ECS control
- You have existing ECS infrastructure
- You need faster task startup (no Batch queue overhead)

Choose `awsbatch` when:
- You need task arrays for large parallel jobs
- You need tasks running longer than 14 days
- You want AWS Batch's managed compute environments
