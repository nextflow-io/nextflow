# Feature Specification: AWS ECS Managed Instances Executor

**Feature Branch**: `001-ecs-executor`
**Created**: 2025-12-30
**Status**: Draft
**Input**: User description: "Add an executor to offload task executions via AWS ECS Managed instances, see https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ManagedInstances.html. The executor should make part of the nf-amazon plugin, however it should be independent from the AWS Batch implementation. The executor should allow running tasks as containers and should allow definition the following compute resources: CPUs, memory, GPUs, Disk size and type and the EC2 instance type (optionally). The storage strategy should be based on AWS S3 via Seqera Fusion file system."

## Clarifications

### Session 2025-12-30

- Q: Should the executor auto-create ECS clusters or require pre-configured infrastructure? → A: Pre-configured cluster required; executor validates cluster exists and has capacity provider; documentation provided for cluster setup.
- Q: How should task failure retries be handled? → A: Executor retries spot interruptions automatically (configurable max attempts, similar to aws.batch.maxSpotAttempts); all other errors delegate to Nextflow's errorStrategy/maxRetries.
- Q: How should tasks exceeding the 14-day ECS Managed Instances lifecycle be handled? → A: Fail task with clear error message; document this limitation prominently in executor documentation.
- Q: What configuration namespace should be used for ECS executor settings? → A: `aws.ecs.*` (e.g., `aws.ecs.cluster`, `aws.ecs.maxSpotAttempts`) to parallel existing `aws.batch` pattern.
- Q: What should the executor name be? → A: `awsecs`, following the `awsbatch` naming pattern.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Run Basic Nextflow Task on ECS Managed Instances (Priority: P1)

A pipeline developer wants to run a containerized Nextflow task on AWS ECS Managed Instances with basic compute resources (CPUs and memory). The system should automatically provision appropriate EC2 instances, execute the task in the specified container, and make results available via S3/Fusion.

**Why this priority**: This is the core functionality - without basic task execution, no other features matter. This enables the MVP where users can offload compute to ECS Managed Instances.

**Independent Test**: Can be fully tested by running a simple Nextflow workflow with `executor = 'awsecs'` and verifying task completes successfully with output files available in the S3 work directory.

**Acceptance Scenarios**:

1. **Given** a Nextflow workflow with process specifying `executor = 'awsecs'`, `cpus = 2`, `memory = '4 GB'`, and a container image, **When** the workflow is executed, **Then** the task runs on an ECS Managed Instance with at least the requested resources and completes successfully.

2. **Given** a Nextflow workflow with ECS executor configured, **When** the task fails (non-zero exit code), **Then** the error is properly reported back to Nextflow with the exit code and error logs accessible.

3. **Given** multiple concurrent tasks in a workflow, **When** the workflow executes, **Then** tasks are submitted to ECS and run in parallel according to available capacity.

---

### User Story 2 - Configure GPU-Accelerated Tasks (Priority: P2)

A machine learning pipeline developer wants to run GPU-accelerated tasks (e.g., model training, inference) on ECS Managed Instances with NVIDIA GPUs. They need to specify the number and type of GPUs required.

**Why this priority**: GPU support is essential for ML/AI workloads which are a growing use case. AWS ECS Managed Instances supports GPU instance types through attribute-based instance selection.

**Independent Test**: Can be tested by running a workflow with `accelerator` directive and verifying the task runs on a GPU-enabled instance with correct GPU visibility in the container.

**Acceptance Scenarios**:

1. **Given** a process with `accelerator 1, type: 'nvidia-tesla-t4'`, **When** the task executes, **Then** ECS provisions a GPU-enabled instance and the container has access to the specified GPU.

2. **Given** a process requiring GPU but no GPU instances available, **When** the task is submitted, **Then** the system waits for capacity and eventually provisions appropriate resources (or fails with clear error after timeout).

---

### User Story 3 - Configure Custom Storage for Large Data Tasks (Priority: P2)

A bioinformatics pipeline developer needs to run tasks that process large files requiring significant local disk space. They want to specify disk size to ensure sufficient ephemeral storage for temporary files during task execution.

**Why this priority**: Many scientific workflows require substantial temporary storage for intermediate files. ECS Managed Instances supports 30 GiB to 16 TiB of EBS storage.

**Independent Test**: Can be tested by running a workflow with `disk` directive and verifying the task has access to the requested storage capacity.

**Acceptance Scenarios**:

1. **Given** a process with `disk '500 GB'`, **When** the task executes, **Then** the ECS task has at least 500 GB of ephemeral storage available.

2. **Given** a process specifying disk type (e.g., `disk '100 GB' type: 'gp3'`), **When** the task executes, **Then** the underlying EBS volume uses the specified volume type.

---

### User Story 4 - Specify EC2 Instance Type (Priority: P3)

An advanced user wants to pin tasks to specific EC2 instance types for performance predictability, cost optimization, or specific hardware requirements (e.g., high-memory instances, compute-optimized instances).

**Why this priority**: While automatic instance selection works for most cases, power users need control over instance types for specific workloads or cost management.

**Independent Test**: Can be tested by specifying an instance type in configuration and verifying the task runs on that exact instance type.

**Acceptance Scenarios**:

1. **Given** a process with `machineType 'c6i.2xlarge'`, **When** the task executes, **Then** ECS provisions that specific instance type.

2. **Given** a process without explicit instance type, **When** the task executes, **Then** ECS automatically selects a cost-optimized instance matching the CPU/memory requirements.

---

### User Story 5 - Monitor and Debug Task Execution (Priority: P3)

A pipeline developer needs to monitor task progress and debug failures. They need access to task logs and status information during and after execution.

**Why this priority**: Observability is critical for production workflows but builds on top of core execution functionality.

**Independent Test**: Can be tested by running a workflow and verifying logs appear in CloudWatch and are accessible through standard Nextflow mechanisms.

**Acceptance Scenarios**:

1. **Given** a running task on ECS, **When** the user checks Nextflow logs, **Then** they see the task status (pending, running, completed/failed).

2. **Given** a completed or failed task, **When** the user requests logs, **Then** container stdout/stderr logs are available via CloudWatch Logs integration.

---

### Edge Cases

- What happens when requested resources exceed ECS Managed Instances limits (e.g., >16 vCPUs, >120 GB memory per instance)?
- How does the system handle ECS capacity unavailable scenarios (instance type not available in region)?
- What happens when a task runs longer than 14 days (ECS Managed Instances lifecycle limit)? → Task fails with clear error; limitation documented.
- How does the system handle S3/Fusion connectivity failures during task execution?
- What happens when the container image pull fails (image not found, authentication error)?
- How does the system behave when AWS API rate limits are encountered?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a new executor named `awsecs` that submits Nextflow tasks to AWS ECS using the Managed Instances capacity provider.

- **FR-002**: System MUST allow users to specify CPU requirements via the standard `cpus` process directive, mapped to ECS task definition CPU units.

- **FR-003**: System MUST allow users to specify memory requirements via the standard `memory` process directive, mapped to ECS task definition memory allocation.

- **FR-004**: System MUST allow users to specify GPU requirements via the standard `accelerator` process directive, triggering attribute-based instance selection for GPU-enabled instances.

- **FR-005**: System MUST allow users to specify disk size via the standard `disk` process directive, mapped to ECS storage configuration (30 GiB - 16,384 GiB range).

- **FR-006**: System MUST allow users to optionally specify EC2 instance type for precise control over compute resources.

- **FR-007**: System MUST use AWS S3 as the work directory with Seqera Fusion filesystem for efficient data access.

- **FR-008**: System MUST automatically create ECS task definitions based on process requirements (container image, resources, environment variables).

- **FR-009**: System MUST poll and track task status (PENDING, RUNNING, STOPPED) and report back to Nextflow runtime.

- **FR-010**: System MUST extract and report task exit codes from completed containers.

- **FR-011**: System MUST support task cancellation when the user terminates the workflow.

- **FR-012**: System MUST be implemented as part of the nf-amazon plugin but independent from the existing AWS Batch executor code.

- **FR-013**: System MUST support container image specification via the standard `container` process directive.

- **FR-014**: System MUST integrate with CloudWatch Logs for task output logging.

- **FR-015**: System MUST support ECS task role configuration for AWS service access from within tasks.

- **FR-016**: System MUST validate at startup that the configured ECS cluster exists and has a Managed Instances capacity provider, failing fast with a clear error message if not.

- **FR-017**: Documentation MUST be provided describing how to set up the required ECS cluster with Managed Instances capacity provider.

- **FR-018**: System MUST automatically retry tasks that fail due to spot instance interruption, with configurable maximum retry attempts (similar to `aws.batch.maxSpotAttempts`).

- **FR-019**: System MUST delegate all non-spot failures to Nextflow's standard error handling (`errorStrategy`, `maxRetries` directives).

- **FR-020**: System MUST report a clear error when a task is terminated due to the 14-day ECS Managed Instances lifecycle limit.

- **FR-021**: System MUST use configuration namespace `aws.ecs.*` for all ECS executor settings (e.g., `aws.ecs.cluster`, `aws.ecs.maxSpotAttempts`), consistent with the existing `aws.batch.*` pattern.

### Infrastructure Configuration Requirements

The executor MUST minimize required user configuration by using sensible defaults and auto-discovery.

**Required Settings** (no defaults possible):

- **FR-022**: System MUST require `aws.ecs.cluster` setting specifying the ECS cluster name or ARN with Managed Instances capacity provider.

- **FR-023**: System MUST require `aws.ecs.executionRole` setting specifying the IAM role ARN for ECS task execution (image pull, CloudWatch logs).

**Auto-Discovery** (use defaults when not specified):

- **FR-024**: System MUST auto-discover VPC subnets from the default VPC when `aws.ecs.subnets` is not specified. If no default VPC exists, system MUST fail with a clear error message instructing the user to explicitly configure `aws.ecs.subnets`.

- **FR-025**: System MUST auto-discover the default security group from the default VPC when `aws.ecs.securityGroups` is not specified. If no default VPC exists, system MUST fail with a clear error message instructing the user to explicitly configure `aws.ecs.securityGroups`.

**Default Values**:

- **FR-026**: System MUST use `/aws/ecs/nextflow` as the default CloudWatch log group when `aws.ecs.logsGroup` is not specified.

- **FR-027**: System MUST default `aws.ecs.assignPublicIp` to `true` to enable internet access without requiring NAT gateway configuration.

- **FR-028**: System MUST default `aws.ecs.maxSpotAttempts` to `5` for automatic spot interruption retry.

### Key Entities

- **ECS Cluster**: The AWS ECS cluster with Managed Instances capacity provider where tasks are executed.

- **Task Definition**: ECS task definition created from Nextflow process requirements (container, CPU, memory, GPU, storage).

- **Task**: Running ECS task instance corresponding to a Nextflow task execution.

- **Capacity Provider**: The Managed Instances capacity provider configuration that handles automatic instance provisioning.

- **Work Directory**: S3 bucket path used as the Nextflow work directory, accessed via Fusion filesystem.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users can successfully execute Nextflow workflows using ECS Managed Instances executor with equivalent functionality to AWS Batch executor for container-based tasks.

- **SC-002**: Task submission to first status update takes less than 60 seconds for available capacity scenarios.

- **SC-003**: 99% of successfully completed tasks report correct exit codes back to Nextflow.

- **SC-004**: GPU-accelerated tasks can access the requested number of GPUs within the container environment.

- **SC-005**: Tasks can utilize up to 16 TiB of ephemeral storage when configured via disk directive.

- **SC-006**: Task logs are available in CloudWatch Logs within 30 seconds of task completion.

- **SC-007**: The executor operates independently from AWS Batch code with no shared mutable state.

## Assumptions

- Users have an AWS account with appropriate IAM permissions for ECS, EC2, S3, and CloudWatch.
- An ECS cluster with Managed Instances capacity provider MUST be pre-configured by the user; the executor validates but does not create infrastructure.
- Wave containers service is enabled (required for Fusion filesystem).
- Fusion filesystem is properly configured and licensed for S3 access.
- Container images are accessible from the ECS execution environment (ECR, Docker Hub, or other registry).
- Network configuration (VPC, subnets, security groups) allows ECS tasks to access S3 and pull container images.
- Task definitions use `MANAGED_INSTANCES` as the required compatibility mode.
- Tasks are expected to complete within the 14-day ECS Managed Instances lifecycle limit.
