# Module System for Nextflow

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2025-01-06
- Tags: modules, dsl, registry, versioning, architecture
- Version: 2.2

## Updates

### Version 2.2 (2025-01-06)
- **Structured tool arguments**: Added `args` property to `tools` section for type-safe argument configuration
- **New implicit variables**: `tools.<toolname>.args.<argname>` returns formatted flag+value; `tools.<toolname>.args` returns all args concatenated
- **Deprecation**: All `ext.*` custom directives (e.g., `ext.args`, `ext.args2`, `ext.args3`, `ext.prefix`, `ext.suffix`) deprecated in favor of structured tool arguments

### Version 2.1 (2025-12-11)
- **Unified dependencies**: Consolidated `components`, `dependencies`, and `requires` into single `requires` field
- **New sub-properties**: `requires.modules` and `requires.workflows` for declaring module dependencies
- **Unified version syntax**: `[scope/]name[@constraint]` format across plugins, modules, and workflows
- **Deprecation**: `components` field deprecated (use `requires.modules` instead)

## Context and Problem Statement

Nextflow supports local script inclusion via `include` directive but lacks standardized mechanisms for package management, versioning, and distribution of reusable process definitions. This limits code reuse and reproducibility across the ecosystem.

Discussion/request goes back to at least 2019, see GitHub issues [#1376](https://github.com/nextflow-io/nextflow/issues/1376), [#1463](https://github.com/nextflow-io/nextflow/issues/1463) and [#4122](https://github.com/nextflow-io/nextflow/issues/4112).

## Decision

Implement a module system with four core capabilities:

1. **Remote module inclusion** via registry
2. **Semantic versioning** with dependency resolution
3. **Unified Nextflow Registry** (rebrand existing Nextflow registry)
4. **First-class CLI support** (install, publish, search, list, remove, freeze, run)

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
1. Check `nextflow.config` for declared version
2. Check local `modules/@scope/name/` exists
3. Verify integrity against `.checksum` file
4. If missing, download from registry
5. Warn if local module exist and checksum mismatch
6. Warn if transitive dependencies are not pinned

**Resolution Timing**: Modules resolved at workflow parse time (after plugin resolution at startup).

**Local Storage**: Downloaded modules stored in `modules/@scope/name/` directory in project root (not global cache). Each module must contain a `main.nf` file as the required entry point. It is intended that module source code will be committed to the pipeline git repository.

### 2. Semantic Versioning and Configuration

**Version Format**: MAJOR.MINOR.PATCH
- **MAJOR**: Breaking changes to process signatures, inputs, or outputs
- **MINOR**: New processes, backward-compatible enhancements
- **PATCH**: Bug fixes, documentation updates

**Workflow Configuration** (`nextflow.config`):
```groovy
// Module versions (exact versions only, no ranges)
modules {
    '@nf-core/salmon' = '1.1.0'
    '@nf-core/bwa-align' = '1.2.0'
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
name: nf-core/bwa-align
version: 1.2.4                    # This module's version

requires:
  nextflow: ">=24.04.0"
  modules:                        # Required modules (version constraints)
    - nf-core/samtools/view@>=1.0.0,<2.0.0
    - nf-core/samtools/sort@>=2.1.0,<2.2.0
```

**Version Constraints** (unified `name@constraint` syntax):
- `name`: Any version (latest)
- `name@1.2.3`: Exact version
- `name@>=1.2.3`: Greater or equal
- `name@>=1.2.3,<2.0.0`: Range (comma-separated)

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
- Use `nextflow module freeze` to pin all transitive dependencies to exact versions

### 3. Unified Nextflow Registry

**Architecture Decision**: Extend existing Nextflow registry at `registry.nextflow.io` to host both plugins and modules.

**Current Plugin API** (reference: https://registry.nextflow.io/openapi/):
```
GET  /api/v1/plugins                                  # List/search plugins
GET  /api/v1/plugins/{pluginId}                       # Get plugin + all releases
GET  /api/v1/plugins/{pluginId}/{version}             # Get specific release
GET  /api/v1/plugins/{pluginId}/{version}/download/{fileName}  # Download artifact
POST /api/v1/plugins/release                          # Create draft release
POST /api/v1/plugins/release/{releaseId}/upload       # Upload artifact
```

**Module API** (reference: https://github.com/seqeralabs/plugin-registry/pull/266):
```
GET  /api/modules?query=<text>               # Search modules (semantic search)
GET  /api/modules/{name}                     # Get module + latest release
GET  /api/modules/{name}/releases            # List all releases
GET  /api/modules/{name}/{version}           # Get specific release
GET  /api/modules/{name}/{version}/download  # Download module bundle
POST /api/modules/{name}                     # Publish module version (authenticated)
```

Note: The `{name}` parameter includes the namespace prefix (e.g., "nf-core/fastqc").

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
nextflow module freeze                      # Pin all transitive dependencies to config
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
- `--tools:<toolname>:<argname> <value>`: Configure tool-specific arguments (validated against meta.yaml schema)
- All standard `nextflow run` options (e.g., `-profile`, `-work-dir`, `-resume`, etc.)

**Behavior**:
1. Checks if module is installed locally; if not, downloads from registry
2. Parses the module's `main.nf` to identify the main process and its input declarations
3. Validates command-line arguments against the process input schema
4. Validates tool arguments against the `tools.*.args` schema in `meta.yaml`
5. Generates an implicit workflow that wires CLI arguments to process inputs
6. Executes the workflow using standard Nextflow runtime

**Input Mapping**:
- Named arguments (`--reads`, `--reference`) are mapped to corresponding process inputs
- File paths are automatically converted to file channels
- Multiple values can be provided for inputs expecting collections
- Required inputs without defaults must be provided; optional inputs use declared defaults

**Tool Arguments**:
- Arguments prefixed with `--tools:` configure tool-specific parameters
- Format: `--tools:<toolname>:<argname> <value>` (e.g., `--tools:bwa:K 100000000`)
- Boolean flags can be specified without value (e.g., `--tools:bwa:Y`)
- Arguments are validated against the tool's `args` schema in `meta.yaml`
- Invalid argument names or values that fail type/enum validation produce errors

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

# Run with tool-specific arguments
nextflow module run nf-core/bwa-align \
    --reads 'samples/*_{1,2}.fastq.gz' \
    --reference genome.fa \
    --tools:bwa:K 100000000 \
    --tools:bwa:Y \
    --tools:samtools:output_fmt cram
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
1. If `-version' not specified, resolves the module version from `nextflow.config` or queries registry for latest
2. Downloads the module archive from the registry
3. Extracts to `modules/@scope/name/` directory (replaces existing if version differs)
4. Stores `.checksum` file from registry's X-Checksum response header
5. Recursively installs transitive dependencies declared in `meta.yaml`
6. Updates `nextflow.config` if installing a new module not already configured

**Example**:
```bash
nextflow module install                      # Install all from config
nextflow module install nf-core/bwa-align    # Install specific module (latest)
nextflow module install nf-core/salmon -version 1.2.0
```

---

#### `nextflow module freeze`

Pin all transitive module dependencies to exact versions in `nextflow.config`. This ensures reproducible builds by capturing the precise versions of all dependencies.

**Behavior**:
1. Scans the `modules/` directory for all installed modules
2. Reads version from each module's `meta.yaml`
3. Adds transitive dependencies not explicitly declared to `nextflow.config`

**Output** (in `nextflow.config`):
```groovy
modules {
    '@nf-core/bwa-align' = '1.2.4'
    '@nf-core/samtools/view' = '1.0.5'      // Transitive dependency
    '@nf-core/samtools/sort' = '2.1.0'      // Transitive dependency
}
```

**Example**:
```bash
nextflow module freeze                       # Pin all transitive dependencies
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
1. Removes the module directory from `modules/@scope/name/`
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
name: nf-core/bwa-align
version: 1.2.4              # This module's version
description: Align reads using BWA-MEM
authors:
  - nf-core community
license: MIT

requires:
  nextflow: ">=24.04.0"
  plugins:
    - nf-amazon@2.0.0
  modules:
    - nf-core/samtools/view@>=1.0.0,<2.0.0
    - nf-core/samtools/sort@>=2.1.0,<2.2.0
```

**Local Storage Structure**:
```
project-root/
├── nextflow.config
├── main.nf
└── modules/                        # Local module cache
    ├── @nf-core/
    │   ├── bwa-align/
    │   │   ├── .checksum           # Cached registry checksum
    │   │   ├── meta.yaml
    │   │   └── main.nf             # Required entry point
    │   └── samtools/view/
    │       ├── .checksum
    │       ├── meta.yaml
    │       └── main.nf             # Required entry point
    └── @myorg/
        └── custom-process/
            ├── .checksum
            ├── meta.yaml
            └── main.nf             # Required entry point
```

**Module Integrity Verification**:
- On install: `.checksum` file created from registry's X-Checksum response header
- On run: Local module checksum compared against `.checksum` file
- If match: Proceed without network call
- If mismatch: Report warning (module may have been locally modified)

## Implementation Strategy

**Phase 1**: Module manifest schema, local module loading, validation tools

**Phase 2**: Extend Nextflow registry for modules, implement caching, add `install` and `search` commands

**Phase 3**: Extend DSL parser for `from module` syntax, implement dependency resolution from meta.yaml

**Phase 4**: Implement `publish` command with authentication and `run` command

**Phase 5**: Advanced features (search UI, language server integration, ontology validation)

## Technical Details

**Dependency Resolution Flow**:
1. Parse `include` statements → extract module names (e.g., `@nf-core/bwa-align`)
2. For each module:
   a. Check `nextflow.config` modules section for declared version
   b. Check local `modules/@scope/name/` exists
   c. Verify local module integrity against `.checksum` file
   d. If missing: download from registry; if checksum mismatch: report warning
   e. Warn if module not pinned in config (especially transitive dependencies)
3. On download: store module to `modules/@scope/name/` with `.checksum` file
4. Read module's `meta.yaml` → resolve transitive dependencies recursively
5. Parse module's `main.nf` file → make processes available

**Security**:
- SHA-256 checksum verification on download (stored in `.checksum` file)
- Integrity verification on run (local checksum vs `.checksum` file)
- Authentication required for publishing
- Support for private registries

**Integration with Plugin System**:
- Modules can declare plugin dependencies in meta.yaml
- Both plugins and modules query same registry
- Single authentication system
- Separate cache locations: `$NXF_HOME/plugins/` (global) vs `modules/` (per-project)

## Tool Arguments Configuration

The module system introduces a structured approach to tool argument configuration, replacing the legacy `ext.args` pattern with type-safe, documented argument specifications.

### Current Pattern (Deprecated)

The traditional nf-core pattern uses `ext.args` strings in config files:

```groovy
// Config file
withName: 'BWA_MEM' {
    ext.args = "-K 100000000 -Y -B 3 -R ${meta.read_group}"
    ext.args2 = "--output-fmt cram"
}

// Module script
def args = task.ext.args ?: ''
def args2 = task.ext.args2 ?: ''
bwa mem $args -t $task.cpus $index $reads | samtools sort $args2 -o out.bam -
```

**Limitations:**
- No documentation of available arguments
- No validation or type checking
- Unclear which `ext.argsN` maps to which tool
- No IDE autocompletion support

### New Pattern: Structured Tool Arguments

Modules declare available arguments in `meta.yaml` under each tool's `args` property:

```yaml
tools:
  - bwa:
      description: BWA aligner
      homepage: http://bio-bwa.sourceforge.net/
      args:
        K:
          flag: "-K"
          type: integer
          description: "Process INT input bases in each batch"
        Y:
          flag: "-Y"
          type: boolean
          description: "Use soft clipping for supplementary alignments"

  - samtools:
      description: SAMtools
      homepage: http://www.htslib.org/
      args:
        output_fmt:
          flag: "--output-fmt"
          type: string
          enum: ["sam", "bam", "cram"]
          description: "Output format"
```

### Configuration Usage

Arguments are configured using `tools.<toolname>.args.<argname>`:

```groovy
withName: 'BWA_MEM' {
    tools.bwa.args.K = 100000000
    tools.bwa.args.Y = true
    tools.samtools.args.output_fmt = "cram"
}
```

### Script Usage

In module scripts, access arguments via the `tools` implicit variable:

```groovy
// tools.bwa.args.K → "-K 100000000"
// tools.bwa.args.Y → "-Y"
// tools.bwa.args   → "-K 100000000 -Y"  (all args concatenated)

bwa mem ${tools.bwa.args} -t $task.cpus $index $reads \
    | samtools sort ${tools.samtools.args} -o ${prefix}.bam -
```

### Benefits

| Aspect | `ext.args` (Legacy) | `tools.*.args` (New) |
|--------|---------------------|----------------------|
| Documentation | None | In meta.yaml |
| Type Safety | None | Validated |
| IDE Support | None | Autocompletion |
| Multi-tool | Confusing (`ext.args2`) | Clear (`tools.samtools.args`) |
| Defaults | Manual | Schema-defined |
| Enums | None | Validated |

## Comparison: Plugins vs. Modules

| Aspect | Plugins | Modules |
|--------|---------|---------|
| Purpose | Extend runtime | Reusable processes |
| Format | JAR files | Source code (.nf) |
| Resolution | Startup | Parse time |
| Metadata | JSON spec | YAML manifest |
| Naming | `nf-amazon` | `@nf-core/salmon` |
| Cache Location | `$NXF_HOME/plugins/` | `modules/@scope/name/` |
| Version Config | `plugins {}` in config | `modules {}` in config |
| Registry Path | `/api/v1/plugins/` | `/api/modules/{name}` |

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
- Use `nextflow module freeze` to pin all transitive dependencies when needed
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

## Open Questions

1. **Local vs managed module distinction**: Should local modules use the `@` prefix in include statements, or should a dot file (e.g., `.nf-modules`) be used to distinguish local modules from managed/remote modules?

2. **Tool arguments CLI syntax**: What is the preferred syntax for tool arguments on the command line?
   - Colon-separated: `--tools:<tool>:<arg> <value>`
   - Dot-separated: `--tools.<tool>.<arg> <value>`

3. **Module version configuration**: Should pipeline module versions be specified in `nextflow.config` or in a dedicated pipeline spec file (e.g., `pipeline.yaml`)?

---

## Appendix A: Module Metadata Schema Specification

This appendix defines the JSON schema for module `meta.yaml` files. The schema maintains backward compatibility with existing nf-core module metadata patterns while supporting the new Nextflow module system features.

**Schema File:** [module-spec-schema.json](module-spec-schema.json)
**Published URL:** `https://registry.nextflow.io/schemas/module-spec/v1.0.0`

### Field Reference

#### Core Fields (Existing nf-core Pattern)

These fields are already widely adopted in the nf-core community and remain fully supported:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Module identifier |
| `description` | string | Yes | Brief description of module functionality |
| `keywords` | array[string] | Recommended | Discovery and categorization keywords |
| `authors` | array[string] | Recommended | Original authors (GitHub handles) |
| `maintainers` | array[string] | Recommended | Current maintainers |
| `tools` | array[object] | Conditional | Software tools wrapped by the module |
| `input` | array/object | Recommended | Input channel specifications |
| `output` | object/array | Recommended | Output channel specifications |

#### Extension Fields (Nextflow Module System)

These fields extend the schema to support the new Nextflow module system:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `version` | string | Registry | Semantic version (MAJOR.MINOR.PATCH) |
| `license` | string | Registry | SPDX license identifier for module code |
| `requires` | object | Optional | All requirements: runtime, plugins, and dependencies |
| `requires.nextflow` | string | Optional | Nextflow version constraint |
| `requires.plugins` | array[string] | Optional | Required Nextflow plugins |
| `requires.modules` | array[string] | Optional | Required modules (processes) |
| `requires.workflows` | array[string] | Optional | Required workflows/subworkflows |

### Detailed Field Specifications

#### `name`

The module name must be a fully qualified scoped identifier in `scope/name` format:

```yaml
name: nf-core/fastqc
name: nf-core/bwa-mem
name: myorg/custom-aligner
```

**Naming Rules:**
- Format: `scope/name` (e.g., `nf-core/salmon`, `myorg/custom`)
- Scope: lowercase alphanumeric with hyphens (organization/owner identifier)
- Name: lowercase alphanumeric with underscores/hyphens (module identifier)
- Pattern: `^[a-z0-9][a-z0-9-]*/[a-z][a-z0-9_-]*$`

**Note:** The `@` prefix is used only in Nextflow DSL `include` statements (e.g., `include { FASTQC } from '@nf-core/fastqc'`) to distinguish registry modules from local file paths. The meta.yaml `name` field should not include the `@` prefix.

#### `version`

Semantic version following [SemVer 2.0.0](https://semver.org/):

```yaml
version: "1.0.0"
version: "2.3.1"
version: "1.0.0-beta.1"
```

**Version Semantics:**
- **MAJOR:** Breaking changes to process signatures, inputs, or outputs
- **MINOR:** New processes, backward-compatible enhancements
- **PATCH:** Bug fixes, documentation updates

**Requirement:** Mandatory for registry-published modules (scoped names in `scope/name` format).

#### `requires`

Specifies all requirements for the module: runtime environment, plugins, and dependencies.

```yaml
requires:
  nextflow: ">=24.04.0"
  plugins:
    - nf-amazon@2.0.0
    - nf-wave@>=1.5.0
  modules:
    - nf-core/fastqc@>=1.0.0
    - nf-core/samtools/sort@>=2.1.0,<3.0.0
    - bwa/mem
  workflows:
    - nf-core/fastq-align-bwa@1.0.0
```

**Unified Version Constraint Syntax:**

All requirements (except `nextflow`) use a unified `name@constraint` format:

| Format | Meaning | Example |
|--------|---------|---------|
| `name` | Any version (latest) | `bwa/mem` |
| `name@1.2.3` | Exact version | `nf-core/fastqc@1.0.0` |
| `name@>=1.2.3` | Greater or equal | `nf-core/fastqc@>=1.0.0` |
| `name@>=1.2.3,<2.0.0` | Range constraint | `nf-core/samtools/sort@>=2.1.0,<3.0.0` |

**`requires.nextflow`** - Nextflow version constraint:
```yaml
requires:
  nextflow: ">=24.04.0"           # minimum version
  nextflow: ">=24.04.0,<25.0.0"   # version range
```

**`requires.plugins`** - Required Nextflow plugins:
```yaml
requires:
  plugins:
    - nf-amazon@2.0.0      # exact version
    - nf-wave@>=1.5.0      # minimum version
    - nf-azure             # any version
```

**`requires.modules`** - Required modules (processes):
```yaml
requires:
  modules:
    - nf-core/fastqc@>=1.0.0              # registry module with constraint
    - nf-core/samtools/sort@>=2.1.0       # nested module path
    - bwa/mem                              # local or registry (no constraint)
```

**`requires.workflows`** - Required workflows/subworkflows:
```yaml
requires:
  workflows:
    - nf-core/fastq-align-bwa@1.0.0       # registry workflow
    - my-local-workflow                    # local workflow
```

**Resolution:**
1. The resolver looks up dependencies locally first, then in configured registries
2. Version constraints are resolved transitively
3. Pinned versions are recorded in `nextflow.config` for reproducibility

#### `tools`

Documents the software tools wrapped by the module, including their command-line arguments:

```yaml
tools:
  - bwa:
      description: BWA aligner
      homepage: http://bio-bwa.sourceforge.net/
      license: ["GPL-3.0-or-later"]
      args:
        K:
          flag: "-K"
          type: integer
          description: "Process INT input bases in each batch"
        Y:
          flag: "-Y"
          type: boolean
          description: "Use soft clipping for supplementary alignments"
```

**Tool Properties:**

| Property | Required | Description |
|----------|----------|-------------|
| `description` | Yes | Tool description |
| `homepage` | One of these | Tool homepage URL |
| `documentation` | One of these | Documentation URL |
| `tool_dev_url` | One of these | Development/source URL |
| `doi` | One of these | Publication DOI |
| `arxiv` | No | arXiv identifier |
| `license` | Recommended | SPDX license(s) |
| `identifier` | Recommended | bio.tools identifier |
| `manual` | No | User manual URL |
| `args` | No | Command-line argument specifications |

**Argument Properties (`args.<name>`):**

The `args` object maps argument names to their specifications. Argument names become accessible in scripts via `tools.<toolname>.args.<argname>`.

| Property | Required | Description |
|----------|----------|-------------|
| `flag` | Yes | CLI flag (e.g., `-K`, `--output-fmt`) |
| `type` | Yes | Data type: `boolean`, `integer`, `float`, `string`, `file`, `path` |
| `description` | Yes | Human-readable description |
| `default` | No | Default value |
| `enum` | No | List of allowed values |
| `required` | No | Whether the argument is mandatory (default: false) |

**Argument Type Behavior:**

| Type | Config Example | Output |
|------|----------------|--------|
| `boolean` | `tools.bwa.args.Y = true` | `-Y` |
| `integer` | `tools.bwa.args.K = 100000` | `-K 100000` |
| `string` | `tools.bwa.args.R = "@RG\tID:s1"` | `-R @RG\tID:s1` |
| `string` + `enum` | `tools.samtools.args.output_fmt = "cram"` | `--output-fmt cram` |

#### `input` and `output`

The schema supports both nf-core patterns to ensure backward compatibility:

**Module Pattern (Tuple-based):**
```yaml
input:
  - - meta:
        type: map
        description: Sample metadata
    - reads:
        type: file
        description: Input FastQ files
        ontologies:
          - edam: "http://edamontology.org/format_1930"
  - - index:
        type: directory
        description: Reference index

output:
  bam:
    - - meta:
          type: map
          description: Sample metadata
      - "*.bam":
          type: file
          description: Aligned BAM file
          pattern: "*.bam"
  versions:
    - versions.yml:
        type: file
        description: Software versions
```

**Subworkflow Pattern (Simplified):**
```yaml
input:
  - ch_reads:
      description: |
        Input FastQ files
        Structure: [ val(meta), [ path(reads) ] ]
  - ch_index:
      description: BWA index files
      type: file

output:
  - bam:
      description: Aligned BAM files
  - versions:
      description: Software versions
```

**Channel Element Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `type` | string | Data type: `map`, `file`, `directory`, `string`, `integer`, `float`, `boolean`, `list`, `val` |
| `description` | string | Human-readable description |
| `pattern` | string | File glob pattern or value pattern |
| `optional` | boolean | Whether input is optional (default: false) |
| `default` | any | Default value if not provided |
| `enum` | array | List of allowed values |
| `ontologies` | array | EDAM or other ontology annotations |

### Migration Guide

#### From nf-core Module to Registry Module

**Before (nf-core local):**
```yaml
name: bwa_mem
description: Align reads using BWA-MEM
keywords:
  - alignment
  - bwa
tools:
  - bwa:
      description: BWA software
      homepage: http://bio-bwa.sourceforge.net/
      license: ["GPL-3.0-or-later"]
      identifier: biotools:bwa
authors:
  - "@drpatelh"
maintainers:
  - "@drpatelh"
input:
  # ... existing input spec
output:
  # ... existing output spec
```

**After (Registry-ready):**
```yaml
name: nf-core/bwa-mem              # Added scope prefix
version: "1.0.0"                    # Added version
description: Align reads using BWA-MEM
keywords:
  - alignment
  - bwa
license: MIT                        # Added module license
requires:                           # Added requirements
  nextflow: ">=24.04.0"
  modules:                          # Added module dependencies (if any)
    - nf-core/samtools/sort@>=1.0.0
tools:
  - bwa:
      description: BWA software
      homepage: http://bio-bwa.sourceforge.net/
      license: ["GPL-3.0-or-later"]
      identifier: biotools:bwa
authors:
  - "@drpatelh"
maintainers:
  - "@drpatelh"
input:
  # ... unchanged
output:
  # ... unchanged
```

#### Schema Validation

Use the schema reference in your `meta.yaml`:

```yaml
# yaml-language-server: $schema=https://registry.nextflow.io/schemas/module-spec/v1.0.0

name: nf-core/my-module
version: "1.0.0"
# ...
```

### Compatibility Matrix

| Feature | nf-core Current | Nextflow Module System |
|---------|-----------------|------------------------|
| Simple names | Yes | Yes (local only) |
| Scoped names | No | Yes (registry) |
| Version field | No | Yes (required for registry) |
| `tools` section | Yes | Yes |
| `components` | Yes (subworkflows) | Deprecated → use `requires.modules` |
| `requires` | No | Yes (unified requirements field) |
| I/O specifications | Yes | Yes |
| Ontologies | Yes | Yes |

### Unsupported nf-core Attributes

The following attributes from the nf-core meta schema are **not supported** in the Nextflow module system:

| Attribute | Reason | Future |
|-----------|--------|--------|
| `extra_args` | Not adopted in practice by nf-core modules | Will be redesigned as part of the `tools` schema attribute to document tool-specific arguments and configuration options |
| `components` | Replaced by unified `requires.modules` | Use `requires.modules` for all module dependencies (local and registry) |

### Complete Examples

#### Minimal nf-core Module

```yaml
name: fastqc
description: Run FastQC on sequenced reads
keywords:
  - quality control
  - qc
  - fastq
tools:
  - fastqc:
      description: FastQC quality metrics
      homepage: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
      license: ["GPL-2.0-only"]
      identifier: biotools:fastqc
authors:
  - "@drpatelh"
maintainers:
  - "@drpatelh"
output:
  html:
    - "*.html":
        type: file
        description: FastQC HTML report
  versions:
    - versions.yml:
        type: file
        description: Software versions
```

#### Full Registry Module

```yaml
name: nf-core/bwa-align
version: "1.2.4"
description: Align reads to reference genome using BWA-MEM algorithm
keywords:
  - alignment
  - mapping
  - bwa
  - bam
  - fastq
license: MIT

requires:
  nextflow: ">=24.04.0"
  plugins:
    - nf-wave@1.5.0
  modules:
    - nf-core/samtools/view@>=1.0.0,<2.0.0
    - nf-core/samtools/sort@>=2.1.0,<2.2.0

tools:
  - bwa:
      description: |
        BWA is a software package for mapping DNA sequences
        against a large reference genome.
      homepage: http://bio-bwa.sourceforge.net/
      documentation: https://bio-bwa.sourceforge.net/bwa.shtml
      doi: 10.1093/bioinformatics/btp324
      license: ["GPL-3.0-or-later"]
      identifier: biotools:bwa

authors:
  - "@nf-core"
maintainers:
  - "@drpatelh"
  - "@maxulysse"

input:
  - - meta:
        type: map
        description: Sample metadata map (e.g., [ id:'sample1', single_end:false ])
    - reads:
        type: file
        description: Input FastQ files
        ontologies:
          - edam: "http://edamontology.org/format_1930"
  - - meta2:
        type: map
        description: Reference metadata
    - index:
        type: directory
        description: BWA index directory
        ontologies:
          - edam: "http://edamontology.org/data_3210"

output:
  bam:
    - - meta:
          type: map
          description: Sample metadata
      - "*.bam":
          type: file
          description: Aligned BAM file
          pattern: "*.bam"
          ontologies:
            - edam: "http://edamontology.org/format_2572"
  versions:
    - versions.yml:
        type: file
        description: Software versions
        pattern: "versions.yml"
```

#### Subworkflow with Module Dependencies

```yaml
name: fastq_align_bwa
description: Align reads with BWA and generate statistics
keywords:
  - alignment
  - bwa
  - samtools
  - statistics

requires:
  modules:
    - bwa/mem
    - samtools/sort
    - samtools/index
    - samtools/stats

authors:
  - "@JoseEspinosa"
maintainers:
  - "@JoseEspinosa"

input:
  - ch_reads:
      description: |
        Input FastQ files
        Structure: [ val(meta), [ path(reads) ] ]
  - ch_index:
      description: BWA index files

output:
  - bam:
      description: Sorted BAM files
  - bai:
      description: BAM index files
  - stats:
      description: Alignment statistics
  - versions:
      description: Software versions
```
