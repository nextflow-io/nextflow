# Data Model: Nextflow Module System Client

**Date**: 2026-01-19
**Feature**: 251117-module-system

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

### 2. ModuleManifest

Parsed representation of `meta.yaml` file.

```groovy
@CompileStatic
class ModuleManifest {
    String name            // e.g., "nf-core/fastqc" (without @)
    String version         // e.g., "1.0.0"
    String description     // Module description
    List<String> keywords  // Discovery keywords
    List<String> authors   // GitHub handles
    List<String> maintainers
    String license         // SPDX identifier

    ModuleRequirements requires
    List<ToolDefinition> tools
    List<InputDefinition> input
    Map<String, OutputDefinition> output
}

@CompileStatic
class ModuleRequirements {
    String nextflow        // Version constraint, e.g., ">=24.04.0"
    List<String> plugins   // e.g., ["nf-amazon@2.0.0"]
    List<String> modules   // e.g., ["nf-core/samtools@>=1.0.0"]
    List<String> workflows // e.g., ["nf-core/fastq-align@1.0.0"]
}

@CompileStatic
class ToolDefinition {
    String name            // Tool identifier
    String description
    String homepage
    String documentation
    String doi
    List<String> license
    String identifier      // bio.tools identifier
    Map<String, ArgDefinition> args
}

@CompileStatic
class ArgDefinition {
    String flag            // CLI flag, e.g., "-K"
    String type            // boolean, integer, float, string, file, path
    String description
    Object defaultValue
    List<Object> enumValues
    boolean required = false
}
```

**Validation Rules**:
- `version`: Must be valid SemVer (MAJOR.MINOR.PATCH)
- `type` in ArgDefinition: Must be one of: boolean, integer, float, string, file, path
- `enumValues`: If present, configured value must be in this list

---

### 3. ModuleInfo

Module metadata returned from registry API.

```groovy
@CompileStatic
class ModuleInfo {
    String name            // e.g., "nf-core/fastqc"
    String version         // Specific version
    String latestVersion   // Latest available
    String description
    String checksum        // SHA-256 of bundle
    long downloadCount
    Instant publishedAt
    List<String> versions  // All available versions
}
```

---

### 4. InstalledModule

Represents a module in local `modules/` directory.

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
    ModuleManifest manifest

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

### 5. ModuleConfig

Module configuration from `nextflow.config`.

```groovy
@CompileStatic
class ModuleConfig {
    Map<String, String> modules = [:]  // name -> version
    RegistryConfig registry
}

@CompileStatic
class RegistryConfig {
    String url = 'https://registry.nextflow.io'
    List<String> urls = []             // Multiple registries
    Map<String, String> auth = [:]     // registry -> token expression
}
```

**Config Syntax**:
```groovy
modules {
    '@nf-core/fastqc' = '1.0.0'
    '@nf-core/bwa-align' = '1.2.0'
}

registry {
    url = 'https://registry.nextflow.io'
    auth {
        'registry.nextflow.io' = '${NXF_REGISTRY_TOKEN}'
    }
}
```

---

### 6. ToolArgsContext

Runtime context for tool arguments in process scripts.

```groovy
@CompileStatic
class ToolArgsContext {
    private Map<String, ToolArgs> tools = [:]

    ToolArgs getAt(String toolName) {
        return tools[toolName]
    }
}

@CompileStatic
class ToolArgs {
    private Map<String, ArgDefinition> schema
    private Map<String, Object> values

    String getAt(String argName) {
        def def = schema[argName]
        def value = values[argName]
        if (value == null) return ''
        if (def.type == 'boolean') {
            return value ? def.flag : ''
        }
        return "${def.flag} ${value}"
    }

    String toString() {
        schema.keySet()
            .findAll { values.containsKey(it) && values[it] != null }
            .collect { this[it] }
            .findAll { it }
            .join(' ')
    }
}
```

---

### 7. ModuleResolutionResult

Result of module resolution process.

```groovy
@CompileStatic
class ModuleResolutionResult {
    ModuleReference reference
    Path resolvedPath          // Absolute path to main.nf
    ResolutionAction action
    String message             // Warning/info message if any
    ModuleManifest manifest
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
ModuleConfig (1) -----> (*) ModuleReference
     |
     v
RegistryConfig (1) -----> (*) Registry URLs

ModuleReference (1) -----> (0..1) InstalledModule
     |
     v (via registry)
ModuleInfo (1) -----> (1) ModuleManifest

InstalledModule (1) -----> (1) ModuleManifest
                    -----> (*) ToolDefinition
                    -----> (*) ArgDefinition

ToolArgsContext (1) -----> (*) ToolArgs
                    -----> (*) ArgDefinition (schema)
```

---

## Storage Layout

```
project-root/
├── nextflow.config              # modules{}, registry{} blocks
├── main.nf                      # include { X } from '@scope/name'
└── modules/
    └── @scope/
        └── name/
            ├── .checksum        # SHA-256 from registry
            ├── main.nf          # Entry point (required)
            ├── meta.yaml        # Manifest (optional but recommended)
            ├── README.md        # Documentation
            └── [other files]    # Supporting files
```

---

## Validation Summary

| Entity | Field | Validation |
|--------|-------|------------|
| ModuleReference | fullName | Pattern: `^@[a-z0-9][a-z0-9-]*/[a-z][a-z0-9_-]*$` |
| ModuleManifest | version | SemVer: `MAJOR.MINOR.PATCH` |
| ArgDefinition | type | Enum: boolean, integer, float, string, file, path |
| ArgDefinition | enumValues | If set, value must be member |
| InstalledModule | directory | Must contain main.nf |
| ModuleConfig | modules | Keys must be valid module references |
| RegistryConfig | url | Valid HTTPS URL |