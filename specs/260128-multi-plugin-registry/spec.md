# Feature Specification: Multiple Nextflow Registry Support

**Feature Branch**: `260128-multi-plugin-registry`
**Created**: 2026-01-28
**Status**: Draft
**Input**: User description: "Organizations running Nextflow in enterprise or air-gapped environments need to host private Nextflow registries while optionally falling back to the public registry. Extend the Nextflow registry configuration to support an ordered list of registries with deterministic resolution order."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Configure Multiple Registries via Config File (Priority: P1)

An enterprise administrator configures Nextflow to check their private Nextflow registry first, then fall back to the public registry when a plugin is not found internally.

**Why this priority**: This is the core use case - enabling organizations to maintain private plugins while still accessing the public ecosystem.

**Independent Test**: Can be fully tested by creating a `nextflow.config` with multiple registry URLs and verifying plugin resolution order.

**Acceptance Scenarios**:

1. **Given** a `nextflow.config` with `plugins.registry = ['https://internal.company.com/plugins', 'https://registry.nextflow.io/api']`, **When** a user runs a workflow requiring plugin `nf-custom@1.0.0` that exists only on the internal registry, **Then** the plugin is downloaded from the internal registry.

2. **Given** a `nextflow.config` with multiple registries configured, **When** a user runs a workflow requiring plugin `nf-amazon@3.0.0` that exists only on the public registry, **Then** the plugin is downloaded from the public registry after the internal registry returns "not found".

3. **Given** a `nextflow.config` with multiple registries where the plugin exists in both, **When** resolving the plugin, **Then** the plugin is downloaded from the first registry in the list (internal) and subsequent registries are not queried.

---

### User Story 2 - Configure Multiple Registries via Environment Variable (Priority: P1)

A CI/CD pipeline operator sets the `NXF_PLUGINS_REGISTRY_URL` environment variable to configure multiple registries without modifying config files.

**Why this priority**: Environment variable configuration is essential for containerized and automated environments.

**Independent Test**: Can be fully tested by setting the environment variable with comma-separated URLs and verifying resolution behavior.

**Acceptance Scenarios**:

1. **Given** `NXF_PLUGINS_REGISTRY_URL=https://internal.company.com/plugins,https://registry.nextflow.io/api`, **When** running a workflow, **Then** registries are queried in the specified order.

2. **Given** both `NXF_PLUGINS_REGISTRY_URL` environment variable and `plugins.registry` config are set, **When** running a workflow, **Then** the environment variable takes precedence and completely overrides the config setting.

---

### User Story 3 - Backward Compatible Single Registry (Priority: P1)

An existing user with a single registry URL configured continues to have their setup work without modification.

**Why this priority**: Backward compatibility ensures zero disruption for existing users.

**Independent Test**: Can be fully tested by using existing single-URL syntax and verifying it works identically to current behavior.

**Acceptance Scenarios**:

1. **Given** a `nextflow.config` with `plugins.registry = 'https://registry.nextflow.io/api'` (single string), **When** running a workflow, **Then** plugin resolution works exactly as before.

2. **Given** `NXF_PLUGINS_REGISTRY_URL=https://registry.nextflow.io/api` (single URL), **When** running a workflow, **Then** plugin resolution works exactly as before.

---

### User Story 4 - Air-Gapped Environment (Priority: P2)

An organization in an air-gapped environment configures only their internal registry without any fallback to external registries.

**Why this priority**: Air-gapped environments are a key enterprise use case but build on the core multi-registry capability.

**Independent Test**: Can be fully tested by configuring only internal registry URLs and verifying no external network requests are made.

**Acceptance Scenarios**:

1. **Given** a `nextflow.config` with `plugins.registry = ['https://internal.company.com/plugins']` (no public registry), **When** a plugin is not found, **Then** resolution fails with a clear error message and no attempt is made to contact external registries.

2. **Given** only internal registries are configured, **When** running in offline mode (`NXF_OFFLINE=true`), **Then** plugins are resolved from local cache only, consistent with current behavior.

---

### User Story 5 - Plugin Resolution Failure Reporting (Priority: P2)

A user receives clear diagnostic information when a plugin cannot be found in any configured registry.

**Why this priority**: Good error messages are essential for troubleshooting but not part of core functionality.

**Independent Test**: Can be fully tested by requesting a non-existent plugin and verifying error message content.

**Acceptance Scenarios**:

1. **Given** multiple registries configured, **When** a plugin is not found in any registry, **Then** the error message lists all registries that were queried.

2. **Given** a registry is unreachable (network error), **When** resolution proceeds to the next registry, **Then** a warning is logged indicating the unreachable registry.

---

### Edge Cases

- What happens when the registry list is empty? System should fail with a clear configuration error.
- What happens when all registries are unreachable? System should fail with an error listing all attempted registries and their failure reasons.
- What happens when a registry returns a malformed response? System should log a warning and proceed to the next registry.
- How does the system handle duplicate URLs in the registry list? Duplicates should be silently deduplicated while preserving order.
- What happens when mixing legacy `.json` format URLs with API URLs? Both formats should be supported in the same registry list.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST accept `plugins.registry` configuration as either a single string (backward compatible) or a list of strings (new functionality).
- **FR-002**: System MUST accept `NXF_PLUGINS_REGISTRY_URL` environment variable as either a single URL or comma-separated list of URLs.
- **FR-003**: System MUST query registries in the order specified when resolving plugins.
- **FR-004**: System MUST stop querying registries once a plugin is found (first match wins).
- **FR-005**: System MUST fall back to the next registry when a plugin is not found (HTTP 404 or equivalent).
- **FR-006**: System MUST fall back to the next registry when a registry is unreachable (network error, timeout).
- **FR-007**: System MUST fail with an error listing all queried registries when a plugin is not found in any registry.
- **FR-008**: System MUST allow environment variable to completely override config file settings (Nextflow convention).
- **FR-009**: System MUST NOT implicitly add the default public registry - users must explicitly include it if desired.
- **FR-010**: System MUST support both HTTP API endpoints (`/api`) and legacy JSON file URLs (`.json`) in the same registry list.
- **FR-011**: System MUST log which registry was used when successfully resolving a plugin (at debug level).
- **FR-012**: System MUST log warnings when a registry is skipped due to errors (unreachable, malformed response).
- **FR-013**: System MUST deduplicate registry URLs while preserving the order of first occurrence.
- **FR-014**: System MUST validate registry URLs (must start with `http://` or `https://`).

### Key Entities

- **Nextflow Registry**: A Nextflow plugin repository endpoint identified by URL, supporting either HTTP API format or legacy JSON format.
- **Registry List**: An ordered collection of registries defining the resolution priority (first registry has highest priority).
- **Plugin Resolution**: The process of locating and downloading a plugin, now iterating through the registry list until found or exhausted.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users can configure multiple registries and plugins resolve from the first registry containing them.
- **SC-002**: Existing single-registry configurations continue to work without modification.
- **SC-003**: Plugin resolution fails gracefully when all registries are exhausted, with actionable error messages.
- **SC-004**: Environment variable configuration overrides config file settings completely.
- **SC-005**: Air-gapped environments can operate with only internal registries, with no attempts to contact external networks.

## Assumptions

- Registry URLs follow existing validation rules (must be HTTP/HTTPS).
- The existing retry logic (3 attempts on connection errors) applies per-registry before falling back to the next.
- The `NXF_PLUGINS_TEST_REPOSITORY` environment variable continues to function as before (appended to the registry list for testing purposes).
- Plugin version resolution logic remains unchanged - only the registry iteration is new.
- The prefetch optimization (single HTTP request for multiple plugins) works per-registry.

## Out of Scope

- Authentication/credentials per registry (may be future enhancement).
- Registry health checks or automatic failover based on latency.
- Caching of "plugin not found" responses to avoid repeated queries.
- GUI or web interface for registry management.
- Registry priority weights (only ordered list is supported).
