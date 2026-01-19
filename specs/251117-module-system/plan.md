# Implementation Plan: Nextflow Module System Client

**Branch**: `251117-module-system` | **Date**: 2026-01-19 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/251117-module-system/spec.md`

## Summary

Implement client-side module system for Nextflow enabling pipeline developers to include remote modules from the Nextflow registry using `@scope/name` syntax, manage versions via `nextflow.config`, and use CLI commands (install, search, list, remove, publish, run). Implementation extends existing DSL parser, config parser, and follows plugin system patterns for registry communication and authentication.

## Technical Context

**Language/Version**: Groovy 4.0.29 (targeting Java 17 runtime, Java 21 toolchain for development)
**Primary Dependencies**:
- Existing Nextflow DSL parser (nf-lang module, ANTLR)
- Existing config parser (ConfigBuilder, ConfigParser)
- Existing HTTP client (HxClient from io.seqera.http)
- Existing plugin authentication infrastructure
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
│   └── CmdModule.groovy                    # NEW: Module CLI command
├── config/
│   ├── ConfigBuilder.groovy                # MODIFY: Add modules/registry DSL
│   └── parser/v1/
│       ├── ModulesDsl.groovy               # NEW: modules {} block parser
│       └── RegistryDsl.groovy              # NEW: registry {} block parser
└── module/
    ├── ModuleResolver.groovy               # NEW: Core resolution logic
    ├── ModuleStorage.groovy                # NEW: Local storage management
    ├── ModuleChecksum.groovy               # NEW: Checksum verification
    ├── ModuleManifest.groovy               # NEW: meta.yaml parser
    └── HttpModuleRepository.groovy         # NEW: Registry HTTP client

modules/nf-lang/src/main/java/nextflow/script/
└── ResolveIncludeVisitor.java              # MODIFY: Add @scope/name detection

modules/nextflow/src/test/groovy/nextflow/
├── cli/
│   └── CmdModuleTest.groovy                # NEW: CLI unit tests
├── config/
│   └── ModulesDslTest.groovy               # NEW: Config parsing tests
└── module/
    ├── ModuleResolverTest.groovy           # NEW: Resolution logic tests
    ├── ModuleStorageTest.groovy            # NEW: Storage tests
    └── ModuleChecksumTest.groovy           # NEW: Checksum tests

tests/
└── modules/                                # NEW: Integration tests
    ├── install-module.nf                   # Test module install + include
    ├── version-resolution.nf               # Test version management
    └── checksum-protection.nf              # Test local modification protection
```

**Structure Decision**: Implementation extends existing Nextflow core modules following modular architecture. New code in `modules/nextflow` for CLI and core logic. DSL parser extension in `modules/nf-lang`. No new plugins required.

## Complexity Tracking

No constitution violations requiring justification.