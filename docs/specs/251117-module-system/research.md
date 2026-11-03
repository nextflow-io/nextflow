# Research: Nextflow Module System Client

**Date**: 2026-01-19
**Feature**: 251117-module-system

## Overview

This document captures technical research and decisions for implementing the Nextflow module system client. All NEEDS CLARIFICATION items from Technical Context have been resolved through codebase exploration.

---

## 1. CLI Command Structure

**Research Question**: How should `nextflow module` CLI commands be implemented?

**Decision**: JCommander native subcommands — each subcommand extends `CmdBase` directly; no trait needed

**Rationale**:
- JCommander's subcommand support handles parameter parsing automatically per subcommand
- Each subcommand (install, run, list, remove, search, info, publish) is a separate class extending CmdBase
- `ModuleRun` extends `CmdRun` to reuse pipeline execution logic (PR #6381)
- No custom `ModuleSubCmd` trait needed; cleaner architecture
- `CmdModule` is registered in `Launcher` alongside all other top-level commands

**Implemented Pattern**:
```groovy
@Parameters(commandDescription = "Manage Nextflow modules")
class CmdModule extends CmdBase implements UsageAware {
    static final List<CmdBase> commands = []

    static {
        commands << new ModuleInstall()   // extends CmdBase
        commands << new ModuleRun()       // extends CmdRun
        commands << new ModuleList()      // extends CmdBase
        commands << new ModuleRemove()    // extends CmdBase
        commands << new ModuleSearch()    // extends CmdBase
        commands << new ModuleInfo()      // extends CmdBase
        commands << new ModulePublish()   // extends CmdBase
    }

    void run() {
        final jc = commander()    // JCommander with all subcommands registered
        jc.parse(args as String[])
        final subcommand = jc.getCommands().get(jc.getParsedCommand()).getObjects()[0]
        subcommand.run()
    }
}
```

**Alternatives Considered**:
- CmdFs trait pattern: Considered initially; replaced by JCommander native subcommands — simpler and avoids custom parsing
- Separate top-level Cmd classes (CmdModuleInstall, etc.): Rejected — too many entry points
- Plugin-based CLI extension: Rejected — module system is core functionality, not optional

---

## 2. DSL Parser Extension for @scope/name

**Research Question**: How to extend `include` statement parsing for registry modules?

**Decision**: Extend `ResolveIncludeVisitor` to detect `@` prefix and delegate to a `RemoteModuleResolver` SPI loaded via Java `ServiceLoader`

**Rationale**:
- Keeps `nf-lang` decoupled from runtime module resolution (`nf-lang` has no dependency on `nextflow` module)
- SPI pattern allows plugins or custom implementations to override the default resolver
- Detection: `source.startsWith('@')` distinguishes registry vs local paths — preserves existing include behavior
- Resolution at parse time (after plugin resolution) per ADR

**Implemented Architecture**:
```
include { X } from '@scope/name'
      ↓
ResolveIncludeVisitor.visitInclude()  [nf-lang]
  source.startsWith("@") → RemoteModuleResolverProvider.getInstance().resolve(source, baseDir)
      ↓
RemoteModuleResolverProvider  [nf-lang]
  Java ServiceLoader discovers implementations; picks highest priority
      ↓
DefaultRemoteModuleResolver  [nextflow module]
  Calls ModuleResolver.installModule(reference, version, autoInstall=true)
  Returns Path to modules/@scope/name/main.nf
```

**Key Files**:
- `modules/nf-lang/src/main/java/nextflow/module/spi/RemoteModuleResolver.java` — SPI interface
- `modules/nf-lang/src/main/java/nextflow/module/spi/RemoteModuleResolverProvider.java` — ServiceLoader singleton
- `modules/nf-lang/src/main/java/nextflow/module/spi/FallbackRemoteModuleResolver.java` — error fallback
- `modules/nf-lang/src/main/java/nextflow/script/control/ResolveIncludeVisitor.java` — MODIFIED
- `modules/nextflow/src/main/groovy/nextflow/module/DefaultRemoteModuleResolver.groovy` — default impl

**Alternatives Considered**:
- New ANTLR grammar token for `@`: Rejected — unnecessary parser complexity
- Direct dependency from nf-lang to nextflow module: Rejected — circular dependency risk; SPI decouples cleanly
- Dot file marker for local modules: Deferred in ADR; current impl uses `@` for registry, `.`/`/` for local

---

## 3. Config Parsing for modules{} and registry{} Blocks

**Research Question**: How to add new config DSL blocks?

**Decision**: Create ModulesConfig and RegistryConfig classes implementing ConfigScope interface

**Rationale**:
- ConfigScope is an ExtensionPoint (pf4j) that ConfigBuilder automatically discovers
- Classes implementing ConfigScope and annotated with @ScopeName are automatically parsed
- No need to modify ConfigBuilder or create custom DSL parsers
- Pattern used throughout Nextflow: FusionConfig, CondaConfig, DockerConfig, etc.
- Provides type safety via @CompileStatic and validation via @ConfigOption

**Reference Implementation**:
```
Location: modules/nextflow/src/main/groovy/nextflow/fusion/FusionConfig.groovy
Pattern:
  @ScopeName("modules")
  @Description("Module version declarations")
  @CompileStatic
  class ModulesConfig implements ConfigScope {
      @ConfigOption
      @Description("Module version mappings")
      final Map<String, String> modules = [:]

      ModulesConfig() {}

      ModulesConfig(Map opts) {
          // Parse from config map
      }
  }
```

**ConfigScope Interface**:
```
Location: modules/nf-lang/src/main/java/nextflow/config/spec/ConfigScope.java
public interface ConfigScope extends ExtensionPoint {}
```

**RegistryConfig Pattern**:
```groovy
@ScopeName("registry")
@Description("Module registry configuration")
@CompileStatic
class RegistryConfig implements ConfigScope {
    static final String DEFAULT_REGISTRY_URL = 'https://registry.nextflow.io/api'

    @ConfigOption
    final Collection<String> url   // One or more URLs in priority order

    @ConfigOption
    final String apiKey            // API key; falls back to NXF_REGISTRY_TOKEN env var

    RegistryConfig() {
        url = [DEFAULT_REGISTRY_URL]
        apiKey = null
    }

    RegistryConfig(Map opts) {
        url = opts.url ?: [DEFAULT_REGISTRY_URL]
        apiKey = opts.apiKey as String
    }

    String getUrl() { url ? url[0] : DEFAULT_REGISTRY_URL }
    Collection<String> getAllUrls() { url ?: [DEFAULT_REGISTRY_URL] }
    String getApiKey() { apiKey ?: SysEnv.get('NXF_REGISTRY_TOKEN') }
}
```

**Integration Point**: ConfigBuilder automatically discovers and parses ConfigScope implementations via ExtensionPoint mechanism

**Alternatives Considered**:
- Custom DSL parsers (ModulesDsl/RegistryDsl): Rejected - unnecessary complexity, ConfigScope pattern handles this automatically
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

**Decision**: Support `NXF_REGISTRY_TOKEN` env var + `registry.apiKey` config field

**Rationale**:
- Environment variable provides CI/CD compatibility
- `apiKey` config field allows explicit token configuration
- Authentication is only applied to the primary (first) registry URL
- Bearer token in Authorization header (standard HTTP auth)

**Implementation**:
```
RegistryConfig.getApiKey() returns:
  1. registry.apiKey config value if set
  2. NXF_REGISTRY_TOKEN environment variable as fallback
  3. null if neither is set (unauthenticated requests)
```

**Config Syntax**:
```nextflow
registry {
    apiKey = '${NXF_REGISTRY_TOKEN}'
}
```

**Alternatives Considered**:
- Per-registry token map (`auth {}` block): Was in initial design; simplified to single `apiKey` since only the primary registry uses authentication
- Secrets file (~/.nextflow/secrets.json): Possible future enhancement
- OAuth flow: Rejected for CLI — token-based simpler

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

> **⚠️ REMOVED FROM ADR** — The tool arguments feature (`tools.<name>.args` in meta.yaml and process config) was removed from the module system ADR. It is not implemented and not planned in the current scope. The `meta.yaml` format used in the actual implementation (`ModuleSpec`) does not include tool/argument definitions.

---

## Summary of Key Decisions

| Area | Decision | Key Reference |
|------|----------|---------------|
| CLI | JCommander subcommands; each extends CmdBase (ModuleRun extends CmdRun) | CmdModule.groovy |
| DSL Parser | SPI pattern — ResolveIncludeVisitor delegates to RemoteModuleResolver; DefaultRemoteModuleResolver bridges to ModuleResolver | ResolveIncludeVisitor.java, RemoteModuleResolver.java |
| Config | ModulesConfig + RegistryConfig (ConfigScope) | FusionConfig.groovy, ConfigScope.java |
| Registry HTTP | ModuleRegistryClient using HxClient + npr-api models | HttpPluginRepository.groovy |
| Authentication | `NXF_REGISTRY_TOKEN` env var or `registry.apiKey` config field (primary registry only) | RegistryConfig.groovy |
| Checksums | SHA-256/SHA-512, `.checksum` file, download integrity via X-Checksum header | ModuleChecksum.groovy |
| Version Storage | `nextflow_spec.json` (auto-managed); `modules {}` in nextflow.config (manual alternative) | PipelineSpec.groovy |
| Version Syntax | Plugin-compatible constraints | VersionNumber class |
| Tool Args | ~~Implicit variable, parse-time validation~~ — **Removed from ADR** | N/A |

---

## Open Items (Deferred)

1. **Local vs managed module distinction**: Resolved — `@` prefix for registry modules only; local paths start with `.` or `/`
2. **Tool arguments**: Removed from ADR — not in scope
3. **Module version location**: Resolved — `nextflow_spec.json` (auto-managed by `module install`); `modules {}` block in `nextflow.config` supported as alternative
4. **DSL parser `@scope/name` include**: ✅ Resolved — SPI pattern implemented (T017a-d)