# Module System for Nextflow

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2025-11-14
- Tags: modules, dsl, registry, versioning, architecture

## Context and Problem Statement

Nextflow supports local script inclusion via `include` directive but lacks standardized mechanisms for package management, versioning, and distribution of reusable process definitions. This limits code reuse and reproducibility across the ecosystem.

## Decision

Implement a module system with four core capabilities:

1. **Remote module inclusion** via registry
2. **Semantic versioning** with dependency resolution
3. **Unified Nextflow Registry** (rebrand existing plugin registry)
4. **First-class CLI support** (pull, push, search, run)

## Core Capabilities

### 1. Remote Module Inclusion

**DSL Syntax**:
```groovy
// Import from registry (version resolved from meta.yaml)
include { BWA_ALIGN } from module 'bwa/align'

// Existing file-based includes remain supported
include { MY_PROCESS } from './modules/my-process.nf'
```

**Version Resolution**: Module versions and dependencies declared in the **module's own meta.yaml**, not in include statements or separate lock files.

**Resolution**: Modules resolved at workflow parse time (after plugin resolution at startup).

**Caching**: Downloaded modules cached in `$NXF_HOME/modules/cache/` organized by name/version.

### 2. Semantic Versioning

**Version Format**: MAJOR.MINOR.PATCH
- **MAJOR**: Breaking changes to process signatures, inputs, or outputs
- **MINOR**: New processes, backward-compatible enhancements
- **PATCH**: Bug fixes, documentation updates

**Version Declaration**: Module versions and dependency constraints declared in **meta.yaml**:
```yaml
name: bwa/align
version: 1.2.4                    # This module's version

dependencies:
  "samtools/view": "^1.0.0"       # Semantic version constraint
  "samtools/sort": "~2.1.0"       # Tilde constraint
  "picard/mark": "3.0.1"          # Exact version
```

**Version Constraints**:
- `1.2.3`: Exact version
- `^1.2.3`: Compatible with >=1.2.3 <2.0.0 (caret - allow minor/patch)
- `~1.2.3`: Compatible with >=1.2.3 <1.3.0 (tilde - allow patch only)

**Dependency Resolution**: When a module is imported, resolve its dependencies recursively using constraints from each module's meta.yaml. Use minimal version selection algorithm (choose minimum version satisfying all constraints).

### 3. Unified Nextflow Registry

**Architecture Decision**: Extend existing plugin registry at `registry.nextflow.io` to host both plugins and modules.

**Current Plugin API** (reference: https://registry.nextflow.io/openapi/):
```
GET  /api/v1/plugins                                  # List/search plugins
GET  /api/v1/plugins/{pluginId}                       # Get plugin + all releases
GET  /api/v1/plugins/{pluginId}/{version}             # Get specific release
GET  /api/v1/plugins/{pluginId}/{version}/download/{fileName}  # Download artifact
POST /api/v1/plugins/release                          # Create draft release
POST /api/v1/plugins/release/{releaseId}/upload       # Upload artifact
```

**Proposed Module API Extension** (same pattern):
```
GET  /api/v1/modules                                  # List/search modules
GET  /api/v1/modules/{moduleId}                       # Get module + all releases
GET  /api/v1/modules/{moduleId}/{version}             # Get specific release
GET  /api/v1/modules/{moduleId}/{version}/download/{fileName}  # Download source archive
POST /api/v1/modules/release                          # Create draft release
POST /api/v1/modules/release/{releaseId}/upload       # Upload module archive
```

**Internal API** (already exists for modules - currently semantic search):
```
GET /internal/modules          # Natural language search
GET /internal/modules/{name}   # Retrieve module metadata
```

**Registry URL**: `registry.nextflow.io`

**Artifact Types**:
- **Plugins**: JAR files with JSON metadata, resolved at startup
- **Modules**: Source archives (.nf + meta.yaml), resolved at parse time

**Benefits**:
- Reuses existing infrastructure (HTTP service, S3 storage, authentication)
- Consistent API patterns for both artifact types
- Operational simplicity (one service vs. two)
- Internal module API already partially implemented

### 4. First-Class CLI Support

**Commands**:
```bash
nextflow module search <query>         # Search registry
nextflow module pull <name>            # Download to cache
nextflow module push <directory>       # Publish to registry (requires auth)
nextflow module run <name> [args]      # Execute module directly
```

**Key Features**:
- `search`: Find modules by name, description, tags
- `pull`: Download module and transitive dependencies
- `push`: Validate meta.yaml, authenticate, upload to registry
- `run`: Execute module with auto-generated CLI flags from meta.yaml

## Module Structure

**Directory Layout**:
```
my-module/
├── meta.yaml    # Module manifest (metadata, dependencies, I/O specs)
├── main.nf      # Process definitions
├── tests/       # Optional test workflows
└── README.md    # Optional documentation
```

**Module Manifest** (`meta.yaml`):
```yaml
name: bwa/align
version: 1.2.4              # This module's version
description: Align reads using BWA-MEM
author: nf-core community
license: MIT

requires:
  nextflow: ">=24.04.0"
  plugins:
    - nf-amazon@2.0.0

processes:
  - name: BWA_ALIGN
    inputs:
      - name: reads
        type: path
        format: fastq
        ontology: edam:format_1930
    outputs:
      - name: bam
        type: path
        format: bam
        ontology: edam:format_2572

dependencies:                # Module dependencies with version constraints
  "samtools/view": "^1.0.0"
  "samtools/sort": "~2.1.0"
```

## Implementation Strategy

**Phase 1**: Module manifest schema, local module loading, validation tools

**Phase 2**: Extend plugin registry for modules, implement caching, add `pull` and `search` commands

**Phase 3**: Extend DSL parser for `from module` syntax, implement dependency resolution from meta.yaml

**Phase 4**: Implement `push` command with authentication and `run` command

**Phase 5**: Advanced features (search UI, language server integration, ontology validation)

## Technical Details

**Dependency Resolution Flow**:
1. Parse `include` statements → extract module names (no version in include)
2. Check local cache for latest or fetch module meta.yaml from registry
3. Read module's meta.yaml → get version and dependencies
4. Recursively resolve dependencies using version constraints in each meta.yaml
5. Download missing modules (minimal version selection)
6. Parse module scripts → make processes available

**Security**:
- SHA-256 checksum verification for all downloads
- Authentication required for publishing
- Support for private registries

**Integration with Plugin System**:
- Modules can declare plugin dependencies in meta.yaml
- Both plugins and modules query same registry
- Single authentication system
- Separate cache locations: `$NXF_HOME/plugins/` vs `$NXF_HOME/modules/`

## Comparison: Plugins vs. Modules

| Aspect | Plugins | Modules |
|--------|---------|---------|
| Purpose | Extend runtime | Reusable processes |
| Format | JAR files | Source code (.nf) |
| Resolution | Startup | Parse time |
| Metadata | JSON spec | YAML manifest |
| Registry Path | `/v1/artifacts/plugins/` | `/v1/artifacts/modules/` |

## Rationale

**Why unified registry?**
- Reuses battle-tested infrastructure (HTTP API, S3, auth)
- Single discovery experience for ecosystem
- Lower operational overhead
- Type-specific handling maintains separation of concerns

**Why versions in meta.yaml instead of lock files?**
- Each module is self-contained with its own dependencies
- Simpler model: no separate lock file to manage
- Module author controls exact versions of dependencies
- Reproducibility guaranteed by module version itself

**Why parse-time resolution?**
- Modules are source code, not compiled artifacts
- Allows inspection/modification for reproducibility
- Enables dependency analysis before execution

**Why semantic versioning?**
- Clear compatibility guarantees
- Enables automated dependency resolution
- Industry standard (npm, cargo, Go modules)

## Consequences

**Positive**:
- Enables ecosystem-wide code reuse
- Reproducible workflows (versions pinned in module meta.yaml)
- Centralized discovery and distribution
- Minimal operational overhead (single registry)
- No separate lock file needed (versions self-contained in modules)

**Negative**:
- Registry becomes critical infrastructure (requires HA setup)
- Type-specific handling adds registry complexity
- Parse-time resolution adds latency to workflow startup

**Neutral**:
- Modules and plugins conceptually distinct but share infrastructure
- Different resolution timing supported by same API

## Links

- Related: [Plugin Spec ADR](20250922-plugin-spec.md)
- Inspired by: [Go Modules](https://go.dev/ref/mod), [npm](https://docs.npmjs.com), [Cargo](https://doc.rust-lang.org/cargo/)
- Related: [nf-core modules](https://nf-co.re/modules)
