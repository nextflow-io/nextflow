# Research: Nextflow Module System Client

**Date**: 2026-01-19
**Feature**: 251117-module-system

## Overview

This document captures technical research and decisions for implementing the Nextflow module system client. All NEEDS CLARIFICATION items from Technical Context have been resolved through codebase exploration.

---

## 1. CLI Command Structure

**Research Question**: How should `nextflow module` CLI commands be implemented?

**Decision**: Follow CmdPlugin pattern with sub-command delegation

**Rationale**:
- CmdPlugin.groovy provides proven pattern for multi-action commands
- Uses JCommander `@Parameters` and `@Parameter` annotations
- Sub-commands (install, search, list, remove, publish, run) handled via positional args
- PluginExecAware interface allows plugin extensibility if needed later

**Reference Implementation**:
```
Location: modules/nextflow/src/main/groovy/nextflow/cli/CmdPlugin.groovy
Pattern:
  - Extends CmdBase
  - @Parameters(commandNames = 'module', commandDescription = '...')
  - @Parameter(names = ['-h', '--help'])
  - args list for sub-command + module name
  - run() method dispatches to install(), search(), etc.
```

**Alternatives Considered**:
- Separate CmdModuleInstall, CmdModuleSearch classes: Rejected - too many entry points, doesn't match existing patterns
- Plugin-based CLI extension: Rejected - module system is core functionality, not optional

---

## 2. DSL Parser Extension for @scope/name

**Research Question**: How to extend `include` statement parsing for registry modules?

**Decision**: Extend ResolveIncludeVisitor to detect `@` prefix and delegate to ModuleResolver

**Rationale**:
- IncludeNode already captures source path as string
- Detection: `source.startsWith('@')` distinguishes registry vs local paths
- Resolution happens at parse time (after plugin resolution) per ADR
- Preserves existing local file include behavior

**Reference Implementation**:
```
Location: modules/nf-lang/src/main/java/nextflow/script/ResolveIncludeVisitor.java
Extension Point: visitInclude() method
Pattern:
  1. Check if source starts with '@'
  2. If yes: call ModuleResolver.resolve(source, configuredVersion)
  3. ModuleResolver returns absolute path to modules/@scope/name/main.nf
  4. Continue with standard include processing
```

**Key Files**:
- `IncludeNode.java` - AST representation
- `IncludeEntryNode.java` - Individual entries
- `ResolveIncludeVisitor.java` - Visitor for resolution

**Alternatives Considered**:
- New ANTLR grammar token for `@`: Rejected - unnecessary parser complexity
- Dot file marker for local modules: Deferred to Open Questions in ADR

---

## 3. Config Parsing for modules{} and registry{} Blocks

**Research Question**: How to add new config DSL blocks?

**Decision**: Create ModulesDsl and RegistryDsl classes following PluginsDsl pattern

**Rationale**:
- PluginsDsl.groovy provides exact template for DSL block handling
- ConfigBuilder already supports dynamic DSL registration
- Groovy's methodMissing enables clean config syntax

**Reference Implementation**:
```
Location: modules/nextflow/src/main/groovy/nextflow/config/parser/v1/PluginsDsl.groovy
Pattern:
  @CompileStatic
  class ModulesDsl {
      private Map<String, String> modules = [:]

      def methodMissing(String name, args) {
          // modules { '@nf-core/fastqc' = '1.0.0' }
          modules[name] = args[0].toString()
      }

      Map<String, String> getModules() { modules }
  }
```

**RegistryDsl Pattern**:
```groovy
class RegistryDsl {
    String url = 'https://registry.nextflow.io'
    List<String> urls = []  // For multiple registries
    Map<String, String> auth = [:]

    void url(String value) { this.url = value }
    void url(List<String> values) { this.urls = values }
    void auth(Closure config) { /* parse auth block */ }
}
```

**Integration Point**: ConfigBuilder.build() instantiates DSL objects

**Alternatives Considered**:
- JSON/YAML config file: Rejected - inconsistent with Nextflow config style
- Dedicated pipeline.yaml: Deferred per ADR Open Questions

---

## 4. Registry HTTP Communication

**Research Question**: How to communicate with module registry API?

**Decision**: Create HttpModuleRepository following HttpPluginRepository pattern

**Rationale**:
- HttpPluginRepository provides robust HTTP client with retry logic
- Uses HxClient from io.seqera.http (already a dependency)
- Handles authentication headers consistently
- Supports connection pooling and timeout configuration

**Reference Implementation**:
```
Location: modules/nf-commons/src/main/nextflow/plugin/HttpPluginRepository.groovy
Pattern:
  class HttpModuleRepository {
      private final URI url
      private final HxClient httpClient
      private final String authToken

      ModuleInfo getModule(String name, String version)
      List<ModuleInfo> search(String query, int limit)
      Path download(String name, String version, Path target)
      void publish(String name, Path bundle)
  }
```

**API Endpoints** (from ADR):
```
GET  /api/modules?query=<text>               # Search
GET  /api/modules/{name}                     # Get module + latest release
GET  /api/modules/{name}/releases            # List all releases
GET  /api/modules/{name}/{version}           # Get specific release
GET  /api/modules/{name}/{version}/download  # Download bundle
POST /api/modules/{name}                     # Publish (authenticated)
```

**Alternatives Considered**:
- Direct HttpClient usage: Rejected - loses retry, pooling benefits
- gRPC protocol: Rejected - registry already uses REST

---

## 5. Authentication Patterns

**Research Question**: How to handle registry authentication?

**Decision**: Support NXF_REGISTRY_TOKEN env var + registry.auth config block

**Rationale**:
- Environment variable provides CI/CD compatibility
- Config block allows per-registry tokens for private registries
- Follows existing plugin auth patterns
- Bearer token in Authorization header (standard HTTP auth)

**Reference Implementation**:
```
Location: modules/nextflow/src/main/groovy/nextflow/cli/CmdAuth.groovy
Pattern:
  1. Check NXF_REGISTRY_TOKEN environment variable
  2. Fall back to registry.auth.'registry.nextflow.io' in config
  3. Add header: Authorization: Bearer <token>
```

**Config Syntax**:
```groovy
registry {
    auth {
        'registry.nextflow.io' = '${NXF_REGISTRY_TOKEN}'
        'private.registry.com' = '${PRIVATE_TOKEN}'
    }
}
```

**Alternatives Considered**:
- Secrets file (~/.nextflow/secrets.json): Possible future enhancement
- OAuth flow: Rejected for CLI - token-based simpler

---

## 6. Checksum Verification

**Research Question**: How to implement module integrity verification?

**Decision**: SHA-256 checksum stored in `.checksum` file, verified on every run

**Rationale**:
- SHA-256 is industry standard, already used for plugin verification
- `.checksum` file stores registry-provided checksum (from X-Checksum header)
- Local checksum computed on-demand and compared
- Mismatch indicates local modification (warn, don't override)

**Implementation Pattern**:
```groovy
class ModuleChecksum {
    static final String ALGORITHM = 'SHA-256'

    static String compute(Path moduleDir) {
        // Hash all files in module directory
        // Exclude .checksum itself
        // Return hex-encoded SHA-256
    }

    static boolean verify(Path moduleDir) {
        def expected = moduleDir.resolve('.checksum').text.trim()
        def actual = compute(moduleDir)
        return expected == actual
    }

    static void save(Path moduleDir, String checksum) {
        moduleDir.resolve('.checksum').text = checksum
    }
}
```

**Checksum Scope**: Covers all files in module directory (main.nf, meta.yaml, README.md, etc.)

**Alternatives Considered**:
- Per-file checksums: Rejected - adds complexity, single checksum sufficient
- MD5: Rejected - SHA-256 more secure

---

## 7. Version Constraint Syntax

**Research Question**: What version constraint syntax to use for module dependencies?

**Decision**: Reuse existing Nextflow plugin version constraint syntax

**Rationale**:
- Already implemented and tested in plugin system
- Users familiar with existing `nextflowVersion` syntax
- Supports ranges, comparisons, exact versions
- No new parser code needed

**Supported Syntax**:
| Notation | Meaning | Example |
|----------|---------|---------|
| `1.2.3` | Exact version | `@nf-core/fastqc@1.0.0` |
| `>=1.2.3` | Greater or equal | `@nf-core/fastqc@>=1.0.0` |
| `<=1.2.3` | Less or equal | `@nf-core/fastqc@<=2.0.0` |
| `>=1.2.0,<2.0.0` | Range | `@nf-core/samtools@>=1.0.0,<2.0.0` |

**Reference**: Version parsing code exists in plugin system; reuse VersionNumber class

**Alternatives Considered**:
- NPM-style `^` and `~`: Rejected - inconsistent with existing Nextflow patterns
- Always latest: Rejected - breaks reproducibility

---

## 8. Tool Arguments Implementation

**Research Question**: How to implement structured tool arguments (`tools.<name>.args`)?

**Decision**: Implement as implicit variable in process scope, validated at parse time

**Rationale**:
- `tools` variable accessible in script block like `task`, `params`
- Validation at parse time catches errors early (per clarification)
- Schema defined in meta.yaml, parsed by ModuleManifest
- Concatenation logic handles flag formatting

**Implementation Pattern**:
```groovy
class ToolArgs {
    private Map<String, ArgDef> schema  // From meta.yaml
    private Map<String, Object> values  // From config

    String getAt(String argName) {
        def def = schema[argName]
        def value = values[argName]
        if (def.type == 'boolean' && value) {
            return def.flag  // e.g., "-Y"
        }
        return "${def.flag} ${value}"  // e.g., "-K 100000"
    }

    String toString() {
        // Concatenate all configured args
        values.collect { name, value ->
            this[name]
        }.join(' ')
    }
}
```

**Config Access**:
```groovy
withName: 'BWA_MEM' {
    tools.bwa.args.K = 100000
    tools.bwa.args.Y = true
}
```

**Script Access**:
```groovy
script:
"""
bwa mem ${tools.bwa.args} -t $task.cpus $index $reads
"""
```

**Alternatives Considered**:
- Runtime validation only: Rejected - late errors waste compute
- String-only values: Rejected - loses type safety benefits

---

## Summary of Key Decisions

| Area | Decision | Key Reference |
|------|----------|---------------|
| CLI | CmdModule extends CmdBase | CmdPlugin.groovy |
| DSL Parser | Extend ResolveIncludeVisitor | ResolveIncludeVisitor.java |
| Config | ModulesDsl + RegistryDsl | PluginsDsl.groovy |
| Registry HTTP | HttpModuleRepository | HttpPluginRepository.groovy |
| Authentication | NXF_REGISTRY_TOKEN + config | CmdAuth.groovy |
| Checksums | SHA-256, .checksum file | Standard Java security |
| Version Syntax | Plugin-compatible constraints | VersionNumber class |
| Tool Args | Implicit variable, parse-time validation | New implementation |

---

## Open Items (Deferred)

These items are noted in the ADR as open questions and do not block implementation:

1. **Local vs managed module distinction**: Whether local modules use `@` prefix or dot file marker
2. **Tool arguments CLI syntax**: Colon vs dot separator (`--tools:bwa:K` vs `--tools.bwa.K`)
3. **Module version location**: nextflow.config vs dedicated pipeline.yaml

Current implementation uses:
- `@` prefix for registry modules only (local paths start with `.` or `/`)
- Colon-separated CLI syntax per ADR assumption
- Versions in nextflow.config per ADR decision