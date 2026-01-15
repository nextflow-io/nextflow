# Feature Specification: Nextflow Module System Client

**Feature Branch**: `251117-module-system`
**Created**: 2026-01-15
**Status**: Draft
**Input**: User description: "Implement Nextflow module system client based on ADR 20251114-module-system.md. Focus on client-side implementation only - CLI commands, DSL parser extensions, dependency resolution, and local storage. Registry backend is assumed to be already implemented."

## Overview

This specification covers the **Nextflow client-side implementation** of the module system, enabling pipeline developers to:
- Include remote modules from the Nextflow registry using `@scope/name` syntax
- Manage module versions through `nextflow.config`
- Use CLI commands to install, search, list, remove, freeze, publish, and run modules
- Configure tool arguments through structured `meta.yaml` definitions

**Out of Scope**: Registry backend implementation (assumed already available at `registry.nextflow.io`)

## User Scenarios & Testing

### User Story 1 - Install and Use Registry Module (Priority: P1)

A pipeline developer wants to use a pre-built module from the Nextflow registry in their workflow without manually downloading or managing module files.

**Why this priority**: This is the core value proposition - enabling code reuse from the ecosystem. Without this, the module system provides no benefit.

**Independent Test**: Can be fully tested by running `nextflow module install nf-core/fastqc` and then executing a workflow that includes the module. Delivers immediate value by enabling module consumption.

**Acceptance Scenarios**:

1. **Given** a new Nextflow project with no modules installed, **When** user runs `nextflow module install nf-core/fastqc`, **Then** the module is downloaded to `modules/@nf-core/fastqc/`, a `.checksum` file is created, and `nextflow.config` is updated with the version
2. **Given** a workflow file with `include { FASTQC } from '@nf-core/fastqc'`, **When** user runs `nextflow run main.nf`, **Then** Nextflow resolves the module from local storage and executes the process
3. **Given** a module version declared in `nextflow.config`, **When** user includes the module, **Then** the declared version is used (not latest)
4. **Given** a module with transitive dependencies in `meta.yaml`, **When** user installs the module, **Then** all transitive dependencies are also installed

---

### User Story 2 - Run Module Directly (Priority: P1)

A user wants to run a module directly from the command line without writing a wrapper workflow.

**Why this priority**: Enables immediate productivity - users can test and execute modules without boilerplate code, essential for AI agents and quick experimentation.

**Independent Test**: Can be tested by running `nextflow module run nf-core/fastqc --input 'data/*.fq'` and verifying the process executes.

**Acceptance Scenarios**:

1. **Given** a module is available (locally or in registry), **When** user runs `nextflow module run nf-core/fastqc --input 'data/*.fastq'`, **Then** the module is executed with the provided inputs mapped to process parameters
2. **Given** a module with tool arguments defined in `meta.yaml`, **When** user runs `nextflow module run nf-core/bwa-align --tools:bwa:K 100000`, **Then** the tool argument is validated and passed to the process
3. **Given** a module is not installed locally, **When** user runs `nextflow module run nf-core/salmon`, **Then** the module is automatically downloaded before execution

---

### User Story 3 - Structured Tool Arguments (Priority: P1)

A module author wants to define typed, documented tool arguments that replace the legacy `ext.args` pattern.

**Why this priority**: Critical for module usability - provides type-safe, documented arguments that enable IDE autocompletion and validation, replacing the opaque `ext.args` pattern.

**Independent Test**: Can be tested by configuring `tools.bwa.args.K = 100000` in config and verifying the argument is applied in the script.

**Acceptance Scenarios**:

1. **Given** a module with `tools.*.args` defined in `meta.yaml`, **When** user configures `tools.bwa.args.K = 100000` in config, **Then** the argument is accessible in scripts as `tools.bwa.args.K` returning `-K 100000`
2. **Given** all tool arguments are configured, **When** script uses `${tools.bwa.args}`, **Then** all configured arguments are concatenated in the output
3. **Given** an argument with enum validation, **When** user provides an invalid value, **Then** a validation error is displayed

---

### User Story 4 - Module Version Management (Priority: P2)

A pipeline developer wants to pin and manage module versions to ensure reproducible workflow executions.

**Why this priority**: Reproducibility is important for scientific workflows - version pinning ensures consistent results.

**Independent Test**: Can be tested by modifying `nextflow.config` module versions and verifying the correct version is used on workflow run.

**Acceptance Scenarios**:

1. **Given** a module is installed at version 1.0.0, **When** user changes `nextflow.config` to specify version 1.1.0 and runs the workflow, **Then** version 1.1.0 is automatically downloaded and replaces the local copy
2. **Given** modules with transitive dependencies, **When** user runs `nextflow module freeze`, **Then** all transitive dependencies are pinned to exact versions in `nextflow.config`
3. **Given** modules installed locally, **When** user runs `nextflow module list`, **Then** configured version, installed version, latest available version, and status are displayed for each module

---

### User Story 5 - Module Integrity Protection (Priority: P2)

A pipeline developer who has locally modified a module (for debugging or customization) wants to be protected from accidentally losing those changes.

**Why this priority**: Protects user work - important for developer experience but not blocking core functionality.

**Independent Test**: Can be tested by modifying a module's `main.nf` locally, then attempting to install a different version and verifying the warning appears.

**Acceptance Scenarios**:

1. **Given** a locally modified module (checksum mismatch with `.checksum`), **When** user tries to install a different version, **Then** Nextflow warns about local modifications and does NOT override
2. **Given** a locally modified module, **When** user runs `nextflow module install -force`, **Then** the local module is replaced with the registry version
3. **Given** a locally modified module, **When** user runs the workflow, **Then** a warning is displayed about checksum mismatch but execution continues

---

### User Story 6 - Remove Module (Priority: P3)

A pipeline developer wants to remove a module they no longer need.

**Why this priority**: Housekeeping feature - useful but not blocking core workflows.

**Independent Test**: Can be tested by running `nextflow module remove nf-core/fastqc` and verifying files are deleted and config is updated.

**Acceptance Scenarios**:

1. **Given** a module is installed, **When** user runs `nextflow module remove nf-core/fastqc`, **Then** the module directory is deleted and the entry is removed from `nextflow.config`
2. **Given** a module is referenced in workflow files, **When** user runs `nextflow module remove`, **Then** a warning is displayed about the reference but removal proceeds
3. **Given** orphaned transitive dependencies exist, **When** user removes a module, **Then** orphaned dependencies are identified and user is prompted to remove them

---

### User Story 7 - Search and Discover Modules (Priority: P3)

A pipeline developer wants to find available modules in the registry that match their analysis needs.

**Why this priority**: Discovery feature - useful but users can find modules through documentation or registry web UI.

**Independent Test**: Can be tested by running `nextflow module search bwa` and verifying results are displayed with name, version, and description.

**Acceptance Scenarios**:

1. **Given** modules exist in the registry, **When** user runs `nextflow module search alignment`, **Then** matching modules are displayed with name, latest version, description, and download count
2. **Given** user wants JSON output for scripting, **When** user runs `nextflow module search fastqc -json`, **Then** results are returned in parseable JSON format
3. **Given** many results exist, **When** user runs `nextflow module search quality -limit 5`, **Then** only 5 results are returned

---

### User Story 8 - Publish Module to Registry (Priority: P3)

A module author wants to publish their module to the Nextflow registry for others to use.

**Why this priority**: Ecosystem contribution feature - important for growth but users can consume modules without publishing capability.

**Independent Test**: Can be tested by creating a valid module structure and running `nextflow module publish -dry-run` to validate.

**Acceptance Scenarios**:

1. **Given** a valid module with `main.nf`, `meta.yaml`, and `README.md`, **When** user runs `nextflow module publish myorg/my-module`, **Then** the module is uploaded to the registry and becomes available for installation
2. **Given** an invalid module (missing required fields), **When** user runs `nextflow module publish`, **Then** validation errors are displayed listing the missing requirements
3. **Given** no authentication configured, **When** user runs `nextflow module publish`, **Then** a clear error message indicates authentication is required

---

### Edge Cases

- What happens when the registry is unreachable during module resolution?
  - Nextflow uses locally cached modules if available, otherwise fails with a clear network error
- How does the system handle circular module dependencies?
  - Dependency resolver detects cycles and fails with an error listing the cycle
- What happens when two modules require incompatible versions of the same dependency?
  - Version conflict is reported with the conflicting requirements
- How are modules resolved when multiple registries are configured?
  - Registries are tried in order; first match wins
- What happens when `meta.yaml` is missing from a module?
  - Module is treated as having no dependencies; basic functionality works
- What happens when local module directory is corrupted or incomplete?
  - Checksum mismatch triggers warning; `-force` allows re-download

## Requirements

### Functional Requirements

#### DSL Parser Extension

- **FR-001**: System MUST recognize `@scope/name` syntax in `include` statements as registry module references
- **FR-002**: System MUST distinguish between local file paths (starting with `.` or `/`) and registry modules (starting with `@`)
- **FR-003**: System MUST resolve module versions from `nextflow.config` `modules {}` block before downloading
- **FR-004**: System MUST parse and validate `meta.yaml` files for module metadata and dependencies

#### Module Resolution

- **FR-005**: System MUST resolve modules at workflow parse time (after plugin resolution)
- **FR-006**: System MUST check local `modules/@scope/name/` directory before querying registry
- **FR-007**: System MUST verify module integrity using `.checksum` file on every run
- **FR-008**: System MUST download modules from registry when not present locally or when version differs
- **FR-009**: System MUST NOT override locally modified modules (checksum mismatch) unless `-force` is used
- **FR-010**: System MUST recursively resolve transitive dependencies from `meta.yaml`
- **FR-011**: System MUST warn when transitive dependencies are not pinned in `nextflow.config`

#### Local Storage

- **FR-012**: System MUST store modules in `modules/@scope/name/` directory structure (single version per module)
- **FR-013**: System MUST create `.checksum` file from registry's X-Checksum header on download
- **FR-014**: System MUST store module's `main.nf`, `meta.yaml`, and supporting files in the module directory

#### CLI Commands

- **FR-015**: System MUST provide `nextflow module install [scope/name]` command to download modules
- **FR-016**: System MUST provide `nextflow module search <query>` command to search the registry
- **FR-017**: System MUST provide `nextflow module list` command to show installed vs configured modules
- **FR-018**: System MUST provide `nextflow module remove scope/name` command to delete modules
- **FR-019**: System MUST provide `nextflow module freeze` command to pin all transitive dependencies
- **FR-020**: System MUST provide `nextflow module publish scope/name` command to upload modules to registry
- **FR-021**: System MUST provide `nextflow module run scope/name` command to execute modules directly

#### Configuration

- **FR-022**: System MUST read module versions from `modules {}` block in `nextflow.config`
- **FR-023**: System MUST support `registry {}` block for configuring registry URL and authentication
- **FR-024**: System MUST support `NXF_REGISTRY_TOKEN` environment variable for authentication
- **FR-025**: System MUST support multiple registry URLs with fallback ordering

#### Tool Arguments

- **FR-026**: System MUST provide `tools.<toolname>.args.<argname>` implicit variable in module scripts
- **FR-027**: System MUST validate tool arguments against `meta.yaml` schema (type, enum)
- **FR-028**: System MUST support boolean, integer, float, string, file, and path argument types
- **FR-029**: System MUST concatenate all tool arguments when `tools.<toolname>.args` is accessed

#### Registry Communication

- **FR-030**: System MUST communicate with registry via documented Module API endpoints
- **FR-031**: System MUST handle authentication using Bearer token in Authorization header
- **FR-032**: System MUST verify SHA-256 checksum on module download

### Key Entities

- **Module**: A reusable Nextflow process definition with `main.nf` entry point, optional `meta.yaml` manifest, and README documentation
- **Module Reference**: A scoped identifier (`@scope/name`) pointing to a registry module
- **Module Manifest (meta.yaml)**: YAML file containing module metadata, version, dependencies, tool arguments schema
- **Checksum File (.checksum)**: Local cache of registry checksum for integrity verification
- **Registry Configuration**: Settings for registry URL, authentication, and fallback ordering

## Success Criteria

### Measurable Outcomes

- **SC-001**: Pipeline developers can install and use a registry module within 5 minutes of starting a new project
- **SC-002**: Module resolution adds less than 2 seconds to workflow startup time when modules are cached locally
- **SC-003**: Users can successfully search, install, and run any module from the registry without reading documentation
- **SC-004**: 100% of module version changes in `nextflow.config` result in automatic module updates without manual intervention
- **SC-005**: Users receive clear, actionable error messages for all failure scenarios (network, validation, authentication)
- **SC-006**: Module authors can publish a new module version within 3 minutes using the CLI
- **SC-007**: Locally modified modules are never accidentally overwritten during normal operations
- **SC-008**: All transitive dependencies can be pinned with a single `nextflow module freeze` command

## Assumptions

- Registry backend is fully implemented and available at `registry.nextflow.io` with the Module API as documented in the ADR
- Existing plugin authentication system can be reused for module registry authentication
- Module bundle size limit of 1MB (uncompressed) is enforced by the registry
- Network connectivity is available for initial module downloads; offline operation uses local cache only
- The `modules/` directory is intended to be committed to the pipeline's git repository
- Version constraints in `meta.yaml` follow the same syntax as existing Nextflow plugin version constraints
- SHA-256 is used for all checksum operations
- Tool arguments CLI syntax uses colon-separated format: `--tools:<tool>:<arg>`

## Dependencies

- Registry backend API (Module API endpoints as specified in ADR)
- Existing Nextflow plugin system (for authentication reuse)
- Existing DSL parser infrastructure (for `include` statement extension)
- Existing config parser (for `modules {}` and `registry {}` blocks)