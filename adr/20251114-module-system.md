# Module System for Nextflow

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2025-12-11
- Tags: modules, dsl, registry, versioning, architecture
- Version: 2.1

## Context and Problem Statement

Nextflow supports local script inclusion via `include` directive but lacks standardized mechanisms for package management, versioning, and distribution of reusable process definitions. This limits code reuse and reproducibility across the ecosystem.

Discussion/request goes back to at least 2019, see GitHub issues [#1376](https://github.com/nextflow-io/nextflow/issues/1376), [#1463](https://github.com/nextflow-io/nextflow/issues/1463) and [#4122](https://github.com/nextflow-io/nextflow/issues/4112).

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
// Import from registry (scoped module name, detected by @scope prefix)
include { BWA_ALIGN } from '@nf-core/bwa-align'

// Existing file-based includes remain supported
include { MY_PROCESS } from './modules/my-process.nf'
```

**Module Naming**: NPM-style scoped packages `@scope/name` (e.g., `@nf-core/salmon`, `@myorg/custom`). Unscoped names (eg. local paths) supported for legacy compatibility. No nested paths with the module are allowed - each module must have a `main.nf` as the entry point.

**Version Resolution**: Module versions pinned in `nextflow.config`. If not specified, use the latest available locally in `modules/` directory, or downloaded and cached in the `modules/` directory.

**Resolution Order**:
1. Check `nextflow.config` for pinned version
2. Check local `modules/@scope/name/` for any cached version
3. Query registry for latest version if not found
4. Warn if transitive dependencies are not pinned

**Resolution Timing**: Modules resolved at workflow parse time (after plugin resolution at startup).

**Local Storage**: Downloaded modules stored in `modules/@scope/name@version/` directory in project root (not global cache). Each module must contain a `main.nf` file as the required entry point. It is intended that module source code will be committed to the pipeline git repository.

### 2. Semantic Versioning and Configuration

**Version Format**: MAJOR.MINOR.PATCH
- **MAJOR**: Breaking changes to process signatures, inputs, or outputs
- **MINOR**: New processes, backward-compatible enhancements
- **PATCH**: Bug fixes, documentation updates

**Workflow Configuration** (`nextflow.config`):
```groovy
// Module versions (exact versions only, no ranges)
modules {
    '@nf-core/salmon' = '1.1.0'              // Simple syntax
    '@nf-core/bwa-align' = [                 // Extended syntax (with checksum)
        version: '1.2.0',
        checksum: 'sha256-abc123...'
    ]
}

// Registry configuration (separate block)
registry {
    url = 'https://registry.nextflow.io'     // Default registry

    // allow the use of multiple registry url for resolving module
    // across custom registries, e.g.
    // url = [ 'https://custom.registry.com', 'https://registry.nextflow.io' ]

    auth {
        'registry.nextflow.io' = '${NXF_REGISTRY_TOKEN}'
        'npm.myorg.com' = '${MYORG_TOKEN}'
    }
}
```

**Module Manifest** (`meta.yaml`):
```yaml
name: "@nf-core/bwa-align"
version: 1.2.4                    # This module's version

dependencies:                     # Transitive dependencies (version constraints)
  "@nf-core/samtools-view": "^1.0.0"
  "@nf-core/samtools-sort": "~2.1.0"
```

**Version Constraints** (in module's meta.yaml only):
- `1.2.3`: Exact version
- `>=1.2.3`: Greater or equal
- `>=1.2.3, <2.0.0`: Range (comma-separated) - equivalent to npm's `^1.2.3`
- `>=1.2.3, <1.3.0`: Range (comma-separated) - equivalent to npm's `~1.2.3`

**Version Notation Consistency**:

Modules use the same version constraint syntax already supported by both `nextflowVersion` and plugins:

| Notation | Meaning | nextflowVersion | Plugins | Modules |
| :---- | :---- | :---- | :---- | :---- |
| 1.2.3 | Exact version | ✓ | ✓ | ✓ |
| >=1.2.3 | Greater or equal | ✓ | ✓ | ✓ |
| <=1.2.3 | Less or equal | ✓ | ✓ | ✓ |
| >1.2.3 | Greater than | ✓ | ✓ | ✓ |
| <1.2.3 | Less than | ✓ | ✓ | ✓ |
| >=1.2, <2.0 | Range (comma) | ✓ | ✓ | ✓ |
| !=1.2.3 | Not equal | ✓ | - | - |
| 1.2+ | >=1.2.x <2.0 | ✓ | - | - |
| 1.2.+ | >=1.2.0 <1.3.0 | ✓ | - | - |
| ~1.2.3 | >=1.2.3 <1.3.0 | - | ✓ | - |

Using comparison operators (`>=`, `<`) with comma-separated ranges provides the same expressive power as
npm-style `^` and `~` notation while maintaining consistency with existing Nextflow version constraint syntax.
This avoids introducing new notation that would require additional parser support.

**Dependency Resolution**:
- Workflow's `nextflow.config` specifies exact versions for direct dependencies
- Transitive dependencies resolved from each module's `meta.yaml` using version constraints
- Warn if transitive dependencies are not pinned in workflow config
- Use `nextflow module freeze` to pin all transitive dependencies with checksums

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
GET  /api/v1/modules                                         # List/search modules
GET  /api/v1/modules/{scope}/{name}                          # Get module + all releases
GET  /api/v1/modules/{scope}/{name}/{version}                # Get specific release
GET  /api/v1/modules/{scope}/{name}/{version}/download       # Download source archive
POST /api/v1/modules/release                                 # Create draft release
POST /api/v1/modules/release/{releaseId}/upload              # Upload module archive
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
nextflow module run scope/name              # Run a module directly without a wrapper script
nextflow module search <query>              # Search registry
nextflow module install [scope/name]        # Install all from config, or specific module
nextflow module freeze                      # Pin all versions + checksums to config
nextflow module list                        # Show installed vs configured
nextflow module remove scope/name           # Remove from config + local cache
nextflow module publish scope/name          # Publish to registry (requires api key)
```

#### `nextflow module run scope/name`

Run a module directly without requiring a wrapper workflow script. This command enables standalone execution of any module by automatically mapping command-line arguments to the module's process inputs. If the module is not available locally, it is automatically installed before execution.

**Arguments**:
- `scope/name`: Module identifier to run (required)

**Options**:
- `-version <ver>`: Run a specific version (default: latest or configured version)
- `--<input_name> <value>`: Map value to the corresponding module process input channel
- All standard `nextflow run` options (e.g., `-profile`, `-work-dir`, `-resume`, etc.)

**Behavior**:
1. Checks if module is installed locally; if not, downloads from registry
2. Parses the module's `main.nf` to identify the main process and its input declarations
3. Validates command-line arguments against the process input schema
4. Generates an implicit workflow that wires CLI arguments to process inputs
5. Executes the workflow using standard Nextflow runtime

**Input Mapping**:
- Named arguments (`--reads`, `--reference`) are mapped to corresponding process inputs
- File paths are automatically converted to file channels
- Multiple values can be provided for inputs expecting collections
- Required inputs without defaults must be provided; optional inputs use declared defaults

**Example**:
```bash
# Run BWA alignment module with input files
nextflow module run nf-core/bwa-align \
    --reads 'samples/*_{1,2}.fastq.gz' \
    --reference genome.fa

# Run a specific version with Nextflow options
nextflow module run nf-core/fastqc -version 1.0.0 \
    --input 'data/*.fastq.gz' \
    -profile docker \
    -resume

# Run with work directory and output specification
nextflow module run nf-core/salmon \
    --reads reads.fq \
    --index salmon_index \
    -work-dir /tmp/work \
    --outdir results/
```

---

#### `nextflow module search <query>`

Search the Nextflow registry for available modules matching the specified query. The search operates against module names, descriptions, tags, and author information. Results are displayed with module name, latest version, description, and download statistics.

**Arguments**:
- `<query>`: Search term (required) - matches against module metadata

**Options**:
- `-limit <n>`: Maximum number of results to return (default: 10)
- `-json`: Output results in JSON format for programmatic use

**Example**:
```bash
nextflow module search bwa
nextflow module search "alignment" -limit 50
```

---

#### `nextflow module install [scope/name]`

Download and install modules to the local `modules/` directory. When called without arguments, installs all modules declared in `nextflow.config`. When a specific module is provided, installs that module and adds it to the configuration.

**Arguments**:
- `[scope/name]`: Optional module identifier. If omitted, installs all modules from config

**Options**:
- `-version <ver>`: Install a specific version (default: latest)
- `-force`: Re-download even if already installed locally

**Behavior**:
1. Resolves the module version from `nextflow.config` or queries registry for latest
2. Downloads the module archive from the registry
3. Extracts to `modules/@scope/name@version/` directory
4. Recursively installs transitive dependencies declared in `meta.yaml`
5. Updates `nextflow.config` if installing a new module not already configured

**Example**:
```bash
nextflow module install                      # Install all from config
nextflow module install nf-core/bwa-align    # Install specific module (latest)
nextflow module install nf-core/salmon -version 1.2.0
```

---

#### `nextflow module freeze`

Lock all module versions by writing exact versions and SHA-256 checksums to `nextflow.config`. This ensures fully reproducible builds by capturing the precise state of all dependencies.

**Options**:
- `-verify`: Verify existing checksums without updating

**Behavior**:
1. Scans the `modules/` directory for all installed modules
2. Computes SHA-256 checksums for each module archive
3. Converts simple version syntax to extended syntax with checksums in `nextflow.config`
4. Includes transitive dependencies not explicitly declared

**Output** (in `nextflow.config`):
```groovy
modules {
    '@nf-core/bwa-align' = [
        version: '1.2.4',
        checksum: 'sha256-a1b2c3d4e5f6...'
    ]
}
```

**Example**:
```bash
nextflow module freeze                       # Pin all versions + checksums
nextflow module freeze -verify               # Verify checksums match
```

---

#### `nextflow module list`

Display the status of all modules, comparing what is configured in `nextflow.config` against what is actually installed in the `modules/` directory.

**Options**:
- `-json`: Output in JSON format
- `-outdated`: Only show modules with available updates

**Output columns**:
- Module name (`@scope/name`)
- Configured version (from `nextflow.config`)
- Installed version (from `modules/` directory)
- Latest available version (from registry)
- Status indicator (up-to-date, outdated, missing, not configured)

**Example**:
```bash
nextflow module list
nextflow module list -outdated
```

---

#### `nextflow module remove scope/name`

Remove a module from both the local `modules/` directory and the `nextflow.config` configuration. Also removes orphaned transitive dependencies that are no longer required by other modules.

**Arguments**:
- `scope/name`: Module identifier to remove (required)

**Options**:
- `-keep-config`: Remove local files but keep the entry in `nextflow.config`
- `-keep-files`: Remove from config but keep local files

**Behavior**:
1. Removes the module directory from `modules/@scope/name@version/`
2. Removes the module entry from the `modules {}` block in `nextflow.config`
3. Identifies and optionally removes orphaned transitive dependencies
4. Warns if the module is still referenced in workflow files

**Example**:
```bash
nextflow module remove nf-core/bwa-align
nextflow module remove myorg/custom -keep-files
```

---

#### `nextflow module publish scope/name`

Publish a module to the Nextflow registry, making it available for others to install. Requires authentication via API key and appropriate permissions for the target scope.

**Arguments**:
- `scope/name`: Module identifier to publish (required)

**Options**:
- `-registry <url>`: Target registry URL (default: `registry.nextflow.io`)
- `-tag <tag>`: Additional tags for discoverability
- `-dry-run`: Validate without publishing

**Behavior**:
1. Validates `meta.yaml` schema and required fields (name, version, description)
2. Verifies `main.nf` exists and is valid Nextflow syntax
3. Checks that `README.md` documentation is present
4. Authenticates with registry using configured credentials
5. Creates a release draft and uploads the module archive
6. Publishes the release, making it available for installation

**Requirements**:
- Valid `meta.yaml` with name, version, and description
- `main.nf` entry point file
- `README.md` documentation
- Authentication token configured in `registry.auth` or `NXF_REGISTRY_TOKEN`
- Write permission for the target scope

**Example**:
```bash
nextflow module publish myorg/my-process
nextflow module publish myorg/my-process -dry-run
```

**General Notes**:
- All commands respect the `registry.url` configuration for custom registries
- Unpinned transitive dependencies generate warnings during resolution
- Modules are automatically downloaded on `nextflow run` if missing but configured

## Module Structure

**Directory Layout**:
Everything within the module directory should be uploaded. Module bundle should not exceed 1MB (uncompressed). Typically this is expected to look something like this:
```
my-module/
├── main.nf      # Required: entry point for module
├── meta.yaml    # Optional: Module spec (metadata, dependencies, I/O specs)
├── README.md    # Required: Module description
└── tests/       # Optional test workflows
```

**Module Spec extension** (`meta.yaml`):
```yaml
name: "@nf-core/bwa-align"
version: 1.2.4              # This module's version
description: Align reads using BWA-MEM
author: nf-core community
license: MIT

requires:
  nextflow: ">=24.04.0"
  plugins:
    - nf-amazon@2.0.0

dependencies:
  "@nf-core/samtools-view": "^1.0.0"
  "@nf-core/samtools-sort": "~2.1.0"
```

**Local Storage Structure**:
```
project-root/
├── nextflow.config
├── main.nf
└── modules/                        # Local module cache
    ├── @nf-core/
    │   ├── bwa-align@1.2.4/
    │   │   ├── meta.yaml
    │   │   └── main.nf             # Required entry point
    │   └── samtools-view@1.0.5/
    │       ├── meta.yaml
    │       └── main.nf             # Required entry point
    └── @myorg/
        └── custom-process@2.0.0/
            ├── meta.yaml
            └── main.nf             # Required entry point
```

## Implementation Strategy

**Phase 1**: Module manifest schema, local module loading, validation tools

**Phase 2**: Extend plugin registry for modules, implement caching, add `pull` and `search` commands

**Phase 3**: Extend DSL parser for `from module` syntax, implement dependency resolution from meta.yaml

**Phase 4**: Implement `push` command with authentication and `run` command

**Phase 5**: Advanced features (search UI, language server integration, ontology validation)

## Technical Details

**Dependency Resolution Flow**:
1. Parse `include` statements → extract module names (e.g., `@nf-core/bwa-align`)
2. For each module:
   a. Check `nextflow.config` modules section for pinned version
   b. If not pinned: check local `modules/@scope/name/` for any cached version (use latest local)
   c. If not found locally: query registry for latest version
   d. Warn if module not pinned in config (especially transitive dependencies)
3. Download missing modules to `modules/@scope/name@version/`
4. Read module's `meta.yaml` → resolve transitive dependencies recursively
5. Verify checksums (if present in config)
6. Parse module's `main.nf` file → make processes available

**Security**:
- SHA-256 checksum verification for all downloads
- Authentication required for publishing
- Support for private registries

**Integration with Plugin System**:
- Modules can declare plugin dependencies in meta.yaml
- Both plugins and modules query same registry
- Single authentication system
- Separate cache locations: `$NXF_HOME/plugins/` (global) vs `modules/` (per-project)

## Comparison: Plugins vs. Modules

| Aspect | Plugins | Modules |
|--------|---------|---------|
| Purpose | Extend runtime | Reusable processes |
| Format | JAR files | Source code (.nf) |
| Resolution | Startup | Parse time |
| Metadata | JSON spec | YAML manifest |
| Naming | `nf-amazon` | `@nf-core/salmon` |
| Cache Location | `$NXF_HOME/plugins/` | `modules/@scope/name@version/` |
| Version Config | `plugins {}` in config | `modules {}` in config |
| Registry Path | `/api/v1/plugins/` | `/api/v1/modules/@scope/name` |

## Rationale

**Why unified registry?**
- Reuses battle-tested infrastructure (HTTP API, S3, auth)
- Single discovery experience for ecosystem
- Lower operational overhead
- Type-specific handling maintains separation of concerns

**Why versions in nextflow.config instead of separate lock file?**
- Single source of truth for workflow dependencies
- Simple: exact versions in config, no separate lock file to manage
- Transitive dependencies resolved from module's meta.yaml with version constraints
- Use `nextflow module freeze` to pin all versions + checksums when needed
- Reproducibility via explicit version pinning in config

**Why parse-time resolution?**
- Modules are source code, not compiled artifacts
- Allows inspection/modification for reproducibility
- Enables dependency analysis before execution

**Why NPM-style scoped packages?**
- Organization namespacing prevents name collisions (`@nf-core/salmon` vs `@myorg/salmon`)
- Clear ownership and provenance of modules
- Supports private registries per scope
- Industry-standard pattern (NPM, Terraform, others)
- Enables ecosystem organization by maintainer/organization

**Why semantic versioning?**
- Clear compatibility guarantees
- Enables automated dependency resolution for transitive dependencies
- Industry standard (npm, cargo, Go modules)

## Consequences

**Positive**:
- Enables ecosystem-wide code reuse
- Reproducible workflows (exact versions pinned in nextflow.config)
- Centralized discovery and distribution via unified registry
- Minimal operational overhead (single registry for both plugins and modules)
- NPM-style scoping enables organization namespaces and private registries
- Local `modules/` directory provides project isolation
- Simple config model: no separate lock file unless using `freeze`
- Simple module structure: each module has single `main.nf` entry point

**Negative**:
- Registry becomes critical infrastructure (requires HA setup)
- Type-specific handling adds registry complexity
- Parse-time resolution adds latency to workflow startup
- Local `modules/` directory duplicates storage across projects (unlike global cache)
- Warnings for unpinned transitive dependencies may be noisy initially

**Neutral**:
- Modules and plugins conceptually distinct but share infrastructure
- Different resolution timing supported by same API

## Links

- Related: [Plugin Spec ADR](20250922-plugin-spec.md)
- Inspired by: [Go Modules](https://go.dev/ref/mod), [npm](https://docs.npmjs.com), [Cargo](https://doc.rust-lang.org/cargo/)
- Related: [nf-core modules](https://nf-co.re/modules)
