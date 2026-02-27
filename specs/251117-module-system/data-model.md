# Data Model: Nextflow Module System Client

**Date**: 2026-01-19
**Feature**: 251117-module-system
**Last Updated**: 2026-02-27 (reflects final implementation)

## Overview

This document defines the data entities, their attributes, relationships, and state transitions for the Nextflow module system client implementation.

---

## Entity Definitions

### 1. ModuleReference

Represents a reference to a module in DSL `include` statements.

```groovy
@CompileStatic
class ModuleReference {
    String scope           // e.g., "nf-core"
    String name            // e.g., "fastqc"
    String fullName        // e.g., "@nf-core/fastqc"

    static ModuleReference parse(String source) {
        // Parses "@scope/name" format
    }

    boolean isRegistryModule() {
        return fullName.startsWith('@')
    }
}
```

**Validation Rules**:
- `scope`: lowercase alphanumeric with hyphens, pattern `[a-z0-9][a-z0-9-]*`
- `name`: lowercase alphanumeric with underscores/hyphens, pattern `[a-z][a-z0-9_-]*`
- `fullName`: must match `^@[a-z0-9][a-z0-9-]*/[a-z][a-z0-9_-]*$`

---

### 2. ModuleSpec

Parsed representation of `meta.yaml` file. Class: `nextflow.module.ModuleSpec`.

```groovy
@CompileStatic
class ModuleSpec {
    String name            // e.g., "nf-core/fastqc" (without @)
    String version         // e.g., "1.0.0"
    String description     // Module description
    List<String> keywords  // Discovery keywords
    List<String> authors   // GitHub handles
    String license         // SPDX identifier
    Map<String, String> requires  // dependency -> version constraint

    static ModuleSpec load(Path metaYamlPath) { ... }
    List<String> validate() { ... }  // Returns list of validation errors
    boolean isValid() { ... }
}
```

**Validation Rules**:
- `name`: Must match `scope/name` or `scope/path/to/name` pattern
- `version`: Must be valid SemVer (`MAJOR.MINOR.PATCH[-prerelease]`)
- `description`, `license`: Required fields (validate() reports missing)

**Note**: Tool/argument definitions were removed from the ADR and are not part of `ModuleSpec`.

---

### 3. InstalledModule

Represents a module in the local `modules/` directory.

```groovy
@CompileStatic
class InstalledModule {
    ModuleReference reference
    Path directory         // e.g., /project/modules/@nf-core/fastqc
    Path mainFile          // e.g., /project/modules/@nf-core/fastqc/main.nf
    Path manifestFile      // e.g., /project/modules/@nf-core/fastqc/meta.yaml
    Path checksumFile      // e.g., /project/modules/@nf-core/fastqc/.checksum
    String installedVersion
    String expectedChecksum

    ModuleIntegrity getIntegrity() {
        // Compute and compare checksum
    }
}

enum ModuleIntegrity {
    VALID,              // Checksum matches
    MODIFIED,           // Checksum mismatch (local changes)
    MISSING_CHECKSUM,   // No .checksum file
    CORRUPTED           // Missing required files
}
```

**State Transitions**:
```
[NOT_INSTALLED] --install--> [VALID]
[VALID] --user edits--> [MODIFIED]
[MODIFIED] --install -force--> [VALID]
[VALID] --version change in config--> [VALID] (replaced)
[MODIFIED] --version change in config--> [MODIFIED] (blocked, warn)
```

---

### 4. ModulesConfig and RegistryConfig

Modules configuration loaded from `nextflow_spec.json` ( or the `modules {}` block in `nextflow.config` as alternative). Registry settings from the `registry {}` block in `nextflow.config`.

```groovy
@ScopeName("modules")
@CompileStatic
class ModulesConfig implements ConfigScope {
    Map<String, String> modules = [:]  // module fullName -> version

    String getVersion(String moduleName) { ... }
    boolean hasVersion(String moduleName) { ... }
}

@ScopeName("registry")
@CompileStatic
class RegistryConfig implements ConfigScope {
    static final String DEFAULT_REGISTRY_URL = 'https://registry.nextflow.io/api'

    Collection<String> url    // Registry URL(s) in priority order
    String apiKey             // API key (falls back to NXF_REGISTRY_TOKEN env var)

    String getUrl()           // Returns primary (first) URL
    Collection<String> getAllUrls()
    String getApiKey()        // Returns apiKey or NXF_REGISTRY_TOKEN
}
```

**Config Syntax**:
```nextflow
// nextflow_spec.json (current approach)
{
  "modules": {
    "@nf-core/fastqc": "1.0.0",
    "@nf-core/bwa-align": "1.2.0"
  }
}

// nextflow.config (alternative not currently used)
modules {
    '@nf-core/fastqc' = '1.0.0'
    '@nf-core/bwa-align' = '1.2.0'
}

registry {
    url = [
        'https://private.registry.myorg.com',
        'https://registry.nextflow.io/api'
    ]
    apiKey = '${MYORG_TOKEN}'  // Only applied to the primary registry
}
```

---

### 5. PipelineSpec

Reads and writes `nextflow_spec.json` in the project root. Class: `nextflow.pipeline.PipelineSpec`.

```groovy
class PipelineSpec {
    PipelineSpec(Path baseDir)
    Map<String, String> getModules()
    void addModuleEntry(String name, String version)
    boolean removeModuleEntry(String name)
}
```

---

### 6. ModuleResolutionResult

Result of module resolution process.

```groovy
@CompileStatic
class ModuleResolutionResult {
    ModuleReference reference
    Path resolvedPath          // Absolute path to main.nf
    ResolutionAction action
    String message             // Warning/info message if any
}

enum ResolutionAction {
    USE_LOCAL,          // Used existing local module
    DOWNLOADED,         // Downloaded from registry
    REPLACED,           // Replaced local with different version
    BLOCKED_MODIFIED,   // Local modified, not replaced (warning issued)
    FAILED              // Resolution failed (error)
}
```

---

## Relationships

```
PipelineSpec (1) -----> (*) ModuleReference  (nextflow_spec.json)
ModulesConfig (1) -----> (*) ModuleReference (nextflow.config alternative)
RegistryConfig (1) -----> (*) Registry URLs

ModuleReference (1) -----> (0..1) InstalledModule
     |
     v (via registry)
ModuleSpec (1) <----- InstalledModule (from meta.yaml)
```

---

## Storage Layout

```
project-root/
├── nextflow.config              # registry{} block; optional modules{} block
├── nextflow_spec.json           # auto-managed module version pins
├── main.nf                      # include { X } from '@scope/name'
└── modules/
    └── @scope/
        └── name/
            ├── .checksum        # SHA-256 from registry (download integrity)
            ├── main.nf          # Entry point (required)
            ├── meta.yaml        # Manifest (required for publishing)
            ├── README.md        # Documentation (required for publishing)
            └── [other files]    # Supporting files
```

---

## Validation Summary

| Entity | Field | Validation |
|--------|-------|------------|
| ModuleReference | fullName | Pattern: `^@[a-z0-9][a-z0-9-]*/[a-z][a-z0-9_-]*$` |
| ModuleSpec | name | Pattern: `scope/name` or `scope/path/to/name` |
| ModuleSpec | version | SemVer: `MAJOR.MINOR.PATCH[-prerelease]` |
| ModuleSpec | description, license | Required (non-empty) |
| InstalledModule | directory | Must contain main.nf |
| ModulesConfig | modules keys | Must be valid module fullName |
| RegistryConfig | url | Valid HTTPS URL |