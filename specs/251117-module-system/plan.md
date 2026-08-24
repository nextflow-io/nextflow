# Implementation Plan: Nextflow Module System Client

**Branch**: `251117-module-system` | **Date**: 2026-01-19 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/251117-module-system/spec.md`

## Summary

Implement client-side module system for Nextflow enabling pipeline developers to include remote modules from the Nextflow registry using `@scope/name` syntax, manage versions via `nextflow.config`, configure module parameters via `meta.yaml`, and use CLI commands (install, search, list, remove, publish, run). Implementation extends existing DSL parser, config parser, and follows plugin system patterns for registry communication and authentication.

## Technical Context

**Language/Version**: Groovy 4.0.29 (targeting Java 17 runtime, Java 21 toolchain for development)
**Primary Dependencies**:
- Existing Nextflow DSL parser (nf-lang module, ANTLR)
- Existing config parser (ConfigBuilder, ConfigParser)
- Existing HTTP client (HxClient from io.seqera.http)
- Existing plugin authentication infrastructure
- Existing npr-api (registry data models and schema validation)
**Storage**: Local filesystem (`modules/@scope/name/` per-project, `.checksum` files)
**Testing**: Spock Framework for unit tests, integration tests in `tests/` directory
**Target Platform**: JVM 17+ (same as Nextflow core)
**Project Type**: Multi-module Gradle project extension (core modules + CLI)
**Performance Goals**: Module resolution adds <2 seconds to workflow startup when cached locally (SC-002)
**Constraints**:
- Module bundle size limit: 1MB uncompressed (enforced by registry)
- Backward compatibility: Must not break existing `include` statements
- Offline operation: Must work with locally cached modules
**Scale/Scope**: Ecosystem-wide module distribution; typical project: 5-20 modules

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Modular Architecture | PASS | Module system client belongs in `modules/nextflow` (core CLI) with potential shared utilities in `nf-commons` |
| II. Test-Driven Quality | PASS | Unit tests (Spock), integration tests planned, smoke test support |
| III. Dataflow Programming Model | PASS | Modules are process definitions; include resolution at parse time preserves dataflow semantics |
| IV. Apache 2.0 License | PASS | All new code will include Apache 2.0 headers |
| V. DCO Sign-off | PASS | All commits will use `git commit -s` |
| VI. Semantic Versioning | PASS | Modules use SemVer; plugin-compatible version constraint syntax |
| VII. Groovy Idioms | PASS | Follow existing patterns from CmdPlugin, ConfigBuilder, HttpPluginRepository |

**Gate Status**: PASS - No violations requiring justification

## Project Structure

### Documentation (this feature)

```text
specs/251117-module-system/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
├── contracts/           # Phase 1 output (API contracts)
└── tasks.md             # Phase 2 output (from /speckit.tasks)
```

### Source Code (repository root)

```text
modules/nextflow/src/main/groovy/nextflow/
├── cli/
│   ├── CmdModule.groovy                    # Main module command (uses JCommander)
│   └── module/
│       ├── ModuleInstall.groovy            # Install subcommand (extends CmdBase)
│       ├── ModuleRun.groovy                # Run subcommand (extends CmdRun)
│       ├── ModuleList.groovy               # List subcommand (extends CmdBase)
│       ├── ModuleRemove.groovy             # Remove subcommand (extends CmdBase)
│       ├── ModuleSearch.groovy             # Search subcommand (extends CmdBase)
│       ├── ModuleInfo.groovy               # Info subcommand (extends CmdBase)
│       └── ModulePublish.groovy            # Publish subcommand (extends CmdBase)
├── config/
│   ├── ModulesConfig.groovy                # modules{} config scope
│   └── RegistryConfig.groovy               # registry{} config scope (fields: url, apiKey)
├── module/
│   ├── ModuleReference.groovy              # @scope/name parser
│   ├── ModuleResolver.groovy               # Core resolution logic (version/integrity/install)
│   ├── ModuleStorage.groovy                # Local filesystem operations
│   ├── ModuleRegistryClient.groovy         # HTTP registry client
│   ├── ModuleChecksum.groovy               # SHA-256 integrity verification
│   ├── ModuleSpec.groovy                   # Module manifest (meta.yaml) entity
│   ├── InstalledModule.groovy              # Installed module entity
│   └── DefaultRemoteModuleResolver.groovy  # SPI impl: bridges DSL parser → ModuleResolver
└── pipeline/
    └── PipelineSpec.groovy                 # nextflow_spec.json read/write

modules/nf-lang/src/main/java/nextflow/script/
└── control/ResolveIncludeVisitor.java      # MODIFIED: Delegates @scope/name to SPI resolver

modules/nf-lang/src/main/java/nextflow/module/spi/
├── RemoteModuleResolver.java               # SPI interface (extensible by plugins)
├── RemoteModuleResolverProvider.java       # ServiceLoader wrapper (singleton)
└── FallbackRemoteModuleResolver.java       # Error fallback when no impl found


modules/nextflow/src/test/groovy/nextflow/
├── cli/module/
│   ├── ModuleInstallTest.groovy
│   ├── ModuleRunTest.groovy
│   └── [other subcommand tests]
└── module/
    ├── ModuleResolverTest.groovy
    ├── ModuleStorageTest.groovy
    └── [other module tests]

tests/modules/
├── install-module.nf                       # Integration tests
├── run-module.nf
└── [other integration tests]
```

**Structure Decision**: Implementation extends existing Nextflow core modules following modular architecture. New code in `modules/nextflow` for CLI and core logic. DSL parser extension in `modules/nf-lang` via SPI. No new plugins required.

## Architecture Notes

### Remote Module Inclusion — SPI Pattern

The DSL parser (`ResolveIncludeVisitor`) detects the `@` prefix in `include` statements and delegates resolution to a `RemoteModuleResolver` SPI loaded via Java `ServiceLoader`. This keeps `nf-lang` decoupled from the runtime module resolution logic:

```
include { X } from '@nf-core/fastqc'
      ↓
ResolveIncludeVisitor  (nf-lang)
  source.startsWith("@") → RemoteModuleResolverProvider.getInstance().resolve(...)
      ↓
DefaultRemoteModuleResolver  (nextflow module)
  auto-installs via ModuleResolver if missing → returns Path to main.nf
```

The `RemoteModuleResolver` interface in `nf-lang` can be overridden by plugins with a higher priority value.

## Complexity Tracking

No constitution violations requiring justification.