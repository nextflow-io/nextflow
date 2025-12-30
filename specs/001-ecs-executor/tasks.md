# Tasks: AWS ECS Managed Instances Executor

**Input**: Design documents from `/specs/001-ecs-executor/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/

**Organization**: Tasks are grouped by user story to enable independent implementation and testing.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1-US5)
- Paths are relative to `plugins/nf-amazon/src/`

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Package structure and configuration foundation

- [x] T001 Create ECS executor package directory structure at `main/nextflow/cloud/aws/ecs/`
- [x] T002 Create ECS model subpackage at `main/nextflow/cloud/aws/ecs/model/`
- [x] T003 [P] Create test package at `test/nextflow/cloud/aws/ecs/`
- [x] T004 Register extension point in `plugins/nf-amazon/build.gradle` for AwsEcsExecutor

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core infrastructure that MUST be complete before ANY user story can be implemented

**CRITICAL**: No user story work can begin until this phase is complete

- [x] T005 Create AwsEcsConfig class with @ConfigScope annotation at `main/nextflow/cloud/aws/config/AwsEcsConfig.groovy`
- [x] T006 Add `ecs` field to AwsConfig class at `main/nextflow/cloud/aws/config/AwsConfig.groovy`
- [x] T007 [P] Create AwsEcsOptions wrapper class at `main/nextflow/cloud/aws/ecs/AwsEcsOptions.groovy`
- [x] T008 Create AwsEcsExecutor skeleton with @ServiceName('awsecs') at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [x] T009 Create AwsEcsTaskHandler skeleton extending TaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [x] T010 [P] Create AwsEcsConfigTest at `test/nextflow/cloud/aws/ecs/AwsEcsConfigTest.groovy`

**Checkpoint**: Foundation ready - user story implementation can now begin

---

## Phase 3: User Story 1 - Run Basic Nextflow Task (Priority: P1) MVP

**Goal**: Execute containerized tasks on ECS Managed Instances with basic CPU/memory resources

**Independent Test**: Run simple workflow with `executor = 'awsecs'`, verify task completes with output in S3

### Implementation for User Story 1

- [ ] T011 [US1] Implement AwsEcsConfig required fields (cluster, executionRole) with validation at `main/nextflow/cloud/aws/config/AwsEcsConfig.groovy`
- [ ] T012 [US1] Implement AwsEcsConfig optional fields with defaults (logsGroup, maxSpotAttempts, assignPublicIp) at `main/nextflow/cloud/aws/config/AwsEcsConfig.groovy`
- [ ] T013 [US1] Implement VPC auto-discovery for subnets/securityGroups in AwsEcsConfig (fail with clear error if no default VPC) at `main/nextflow/cloud/aws/config/AwsEcsConfig.groovy`
- [ ] T014 [US1] Implement AwsEcsExecutor.register() with ECS and EC2 client initialization (EC2 for VPC auto-discovery) at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T015 [US1] Implement AwsEcsExecutor.createTaskHandler() at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T016 [US1] Implement AwsEcsExecutor.createTaskMonitor() with ParallelPollingMonitor at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T017 [P] [US1] Create RegisterTaskDefinitionModel at `main/nextflow/cloud/aws/ecs/model/RegisterTaskDefinitionModel.groovy`
- [ ] T018 [P] [US1] Create ContainerDefinitionModel at `main/nextflow/cloud/aws/ecs/model/ContainerDefinitionModel.groovy`
- [ ] T019 [US1] Implement task definition registration with RegisterTaskDefinition API in AwsEcsTaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T020 [US1] Implement CPU/memory mapping (cpus directive -> CPU units, memory -> MiB) in AwsEcsTaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T021 [US1] Implement AwsEcsTaskHandler.submit() with RunTask API at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T022 [US1] Implement AwsEcsTaskHandler.checkIfRunning() with status polling at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T023 [US1] Implement AwsEcsTaskHandler.checkIfCompleted() with exit code extraction at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T024 [US1] Implement AwsEcsTaskHandler.killTask() with StopTask API at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T025 [US1] Implement FusionAwareTask trait for S3/Fusion integration in AwsEcsTaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T026 [US1] Implement task definition caching by container+resource hash in AwsEcsExecutor at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T027 [P] [US1] Create AwsEcsExecutorTest with mock ECS client at `test/nextflow/cloud/aws/ecs/AwsEcsExecutorTest.groovy`
- [ ] T028 [P] [US1] Create AwsEcsTaskHandlerTest at `test/nextflow/cloud/aws/ecs/AwsEcsTaskHandlerTest.groovy`

**Checkpoint**: Basic task execution works - MVP complete

---

## Phase 4: User Story 2 - GPU-Accelerated Tasks (Priority: P2)

**Goal**: Run tasks on GPU-enabled instances using the accelerator directive

**Independent Test**: Run workflow with `accelerator 1, type: 'nvidia-tesla-t4'`, verify GPU access in container

### Implementation for User Story 2

- [ ] T029 [US2] Add GPU resourceRequirements to RegisterTaskDefinitionModel at `main/nextflow/cloud/aws/ecs/model/RegisterTaskDefinitionModel.groovy`
- [ ] T030 [US2] Implement accelerator directive parsing in AwsEcsTaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T031 [US2] Map GPU requirements to ECS resourceRequirements in task definition at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T032 [US2] Add GPU to task definition cache key computation at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T033 [P] [US2] Add GPU tests to AwsEcsTaskHandlerTest at `test/nextflow/cloud/aws/ecs/AwsEcsTaskHandlerTest.groovy`

**Checkpoint**: GPU tasks work independently

---

## Phase 5: User Story 3 - Custom Storage (Priority: P2)

**Goal**: Configure ephemeral disk storage (30 GiB - 16 TiB) via disk directive

**Independent Test**: Run workflow with `disk '500 GB'`, verify storage capacity available

### Implementation for User Story 3

- [ ] T034 [US3] Add storageConfiguration to RegisterTaskDefinitionModel at `main/nextflow/cloud/aws/ecs/model/RegisterTaskDefinitionModel.groovy`
- [ ] T035 [US3] Implement disk directive parsing with validation (30-16384 GiB) in AwsEcsTaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T036 [US3] Map disk size to managedInstancesProvider.storageConfiguration in task definition at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T037 [US3] Add disk size to task definition cache key computation at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T038 [P] [US3] Add disk storage tests to AwsEcsTaskHandlerTest at `test/nextflow/cloud/aws/ecs/AwsEcsTaskHandlerTest.groovy`

**Checkpoint**: Custom storage works independently

---

## Phase 6: User Story 4 - EC2 Instance Type (Priority: P3)

**Goal**: Allow pinning tasks to specific EC2 instance types via machineType directive

**Independent Test**: Run workflow with `machineType 'c6i.2xlarge'`, verify instance type used

### Implementation for User Story 4

- [ ] T039 [US4] Add instance type constraint to RegisterTaskDefinitionModel at `main/nextflow/cloud/aws/ecs/model/RegisterTaskDefinitionModel.groovy`
- [ ] T040 [US4] Implement machineType directive parsing in AwsEcsTaskHandler at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T041 [US4] Map machineType to capacity provider attributes in RunTask request at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T042 [P] [US4] Add instance type tests to AwsEcsTaskHandlerTest at `test/nextflow/cloud/aws/ecs/AwsEcsTaskHandlerTest.groovy`

**Checkpoint**: Instance type selection works independently

---

## Phase 7: User Story 5 - Monitoring & Debugging (Priority: P3)

**Goal**: CloudWatch Logs integration and task status observability

**Independent Test**: Run workflow, verify logs appear in CloudWatch and are accessible via Nextflow

### Implementation for User Story 5

- [ ] T043 [US5] Create AwsEcsHelper class for CloudWatch logs operations at `main/nextflow/cloud/aws/ecs/AwsEcsHelper.groovy`
- [ ] T044 [US5] Implement getTaskLogStream() for log retrieval in AwsEcsHelper at `main/nextflow/cloud/aws/ecs/AwsEcsHelper.groovy`
- [ ] T045 [US5] Implement CloudWatch log configuration in task definition (awslogs driver) at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T046 [US5] Implement describeCluster() for cluster validation in AwsEcsHelper at `main/nextflow/cloud/aws/ecs/AwsEcsHelper.groovy`
- [ ] T047 [US5] Add cluster validation at executor startup using AwsEcsHelper at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T048 [US5] Implement spot interruption detection in checkIfCompleted() at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T049 [US5] Implement spot retry logic with maxSpotAttempts counter at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T050 [US5] Implement 14-day lifecycle error detection and clear error message at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T051 [US5] Implement BatchHandler trait for efficient status polling via BatchContext at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T052 [US5] Implement resource limits validation with clear error for exceeded limits (>16 vCPUs, >120 GB memory) at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T053 [US5] Implement capacity unavailable error handling (instance type not available) at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T054 [US5] Implement container image pull failure detection and clear error message at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`
- [ ] T055 [US5] Implement AWS API rate limit handling with appropriate backoff at `main/nextflow/cloud/aws/ecs/AwsEcsExecutor.groovy`
- [ ] T056 [P] [US5] Add error handling tests to AwsEcsTaskHandlerTest at `test/nextflow/cloud/aws/ecs/AwsEcsTaskHandlerTest.groovy`
- [ ] T057 [P] [US5] Add CloudWatch logs tests to AwsEcsExecutorTest at `test/nextflow/cloud/aws/ecs/AwsEcsExecutorTest.groovy`

**Checkpoint**: Monitoring and debugging works independently

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, integration validation, cleanup

- [ ] T058 [P] Add Apache 2.0 license headers to all new source files
- [ ] T059 [P] Update plugins/nf-amazon/VERSION with appropriate version bump
- [ ] T060 Run quickstart.md validation scenarios
- [ ] T061 Code cleanup and Groovy idiom consistency review
- [ ] T062 [P] Add integration test workflow in tests/ directory
- [ ] T063 [P] [Low Priority] Implement disk type support (e.g., `type: 'gp3'`) in storage configuration at `main/nextflow/cloud/aws/ecs/AwsEcsTaskHandler.groovy`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3-7)**: All depend on Foundational phase completion
  - User stories can proceed in parallel (if staffed)
  - Or sequentially in priority order (P1 -> P2 -> P3)
- **Polish (Phase 8)**: Depends on all desired user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Can start after Foundational (Phase 2) - No dependencies on other stories
- **User Story 2 (P2)**: Can start after Foundational - Extends US1 task definition model
- **User Story 3 (P2)**: Can start after Foundational - Extends US1 task definition model
- **User Story 4 (P3)**: Can start after Foundational - Extends US1 RunTask handling
- **User Story 5 (P3)**: Can start after Foundational - Adds observability to US1 infrastructure

### Within Each User Story

- Models before services/handlers
- Core implementation before enhancements
- Tests can be written in parallel with implementation

### Parallel Opportunities

- All Setup tasks marked [P] can run in parallel
- All Foundational tasks marked [P] can run in parallel
- Once Foundational phase completes, all user stories can start in parallel
- Model classes (T017, T018) can be created in parallel
- Test files can be created in parallel with implementation

---

## Parallel Example: User Story 1

```bash
# Launch model classes in parallel:
Task: "Create RegisterTaskDefinitionModel at main/nextflow/cloud/aws/ecs/model/RegisterTaskDefinitionModel.groovy"
Task: "Create ContainerDefinitionModel at main/nextflow/cloud/aws/ecs/model/ContainerDefinitionModel.groovy"

# Launch tests in parallel with implementation:
Task: "Create AwsEcsExecutorTest at test/nextflow/cloud/aws/ecs/AwsEcsExecutorTest.groovy"
Task: "Create AwsEcsTaskHandlerTest at test/nextflow/cloud/aws/ecs/AwsEcsTaskHandlerTest.groovy"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup
2. Complete Phase 2: Foundational (CRITICAL - blocks all stories)
3. Complete Phase 3: User Story 1
4. **STOP and VALIDATE**: Test basic task execution independently
5. Deploy/demo if ready - users can run basic workflows

### Incremental Delivery

1. Complete Setup + Foundational -> Foundation ready
2. Add User Story 1 -> Test independently -> Deploy (MVP!)
3. Add User Story 2 (GPU) -> Test independently -> Deploy
4. Add User Story 3 (Storage) -> Test independently -> Deploy
5. Add User Story 4 (Instance Type) -> Test independently -> Deploy
6. Add User Story 5 (Monitoring) -> Test independently -> Deploy
7. Each story adds value without breaking previous stories

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together
2. Once Foundational is done:
   - Developer A: User Story 1 (MVP)
   - Developer B: User Story 2 + 3 (after US1 models exist)
3. Stories complete and integrate independently

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story should be independently completable and testable
- Commit after each task or logical group
- Stop at any checkpoint to validate story independently
- Avoid: vague tasks, same file conflicts, cross-story dependencies that break independence
- AwsEcsProxy (throttling) is DEFERRED to optimization phase - not included in these tasks
