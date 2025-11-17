# Nextflow Development Constitution

<!--
SYNC IMPACT REPORT
==================
Version Change: INITIAL → 1.0.0 (Initial constitution)
Modified Principles: N/A (new constitution)
Added Sections:
  - All principles (I-VII)
  - Development Workflow
  - Quality Standards
  - Governance

Removed Sections: N/A (initial version)

Templates Status:
  ✅ plan-template.md - Reviewed, aligned with modular architecture and testing principles
  ✅ spec-template.md - Reviewed, aligned with user scenario focus and requirements structure
  ✅ tasks-template.md - Reviewed, aligned with test-driven and parallel development principles
  ✅ agent-file-template.md - Reviewed, no agent-specific conflicts
  ✅ checklist-template.md - Reviewed, compatible with quality standards

Follow-up TODOs:
  - None at this time
==================
-->

## Core Principles

### I. Modular Architecture

Nextflow MUST maintain a clear separation between core functionality and extensions through its modular architecture:

- **Core modules** (`modules/`) contain essential functionality: workflow engine (nextflow), shared utilities (nf-commons), language parsing (nf-lang), HTTP filesystem support (nf-httpfs), and lineage tracking (nf-lineage)
- **Plugin system** (`plugins/`) provides cloud provider integrations (AWS, Azure, GCP), execution platforms (Kubernetes), and specialized services (Seqera Platform, Wave container management)
- New features MUST be evaluated for placement: core features belong in `modules/`, specialized/cloud-specific features belong in `plugins/`
- Each module and plugin MUST be independently buildable and testable
- Plugin dependencies MUST be explicitly declared in `build.gradle` with semantic versioning

**Rationale**: This architecture enables independent development of cloud provider features without core engine changes, supports third-party plugin development, and maintains a clean separation of concerns across a large multi-module codebase.

### II. Test-Driven Quality Assurance (NON-NEGOTIABLE)

Testing MUST be comprehensive and multi-layered before any code is merged:

- **Unit tests** MUST use Spock Framework for all Groovy code, be independently executable, and achieve meaningful coverage (measured via JaCoCo)
- **Integration tests** (`tests/` directory) MUST validate end-to-end workflows using actual `.nf` scripts with expected outputs
- **Smoke tests** (`make smoke` or `NXF_SMOKE=1`) MUST be available to skip long-running and cloud-dependent tests during rapid development
- **Cloud validation tests** (`validation/` directory) MUST verify cloud provider integrations end-to-end before release
- **Documentation tests** (`docs/snippets/`) MUST ensure all documentation examples remain functional
- All tests MUST pass before commits, and `make test` MUST be run locally before pushing

**Rationale**: Scientific workflows demand reliability and reproducibility. Multi-layered testing catches issues at appropriate levels: unit tests for logic, integration tests for workflow correctness, and validation tests for cloud provider compatibility.

### III. Dataflow Programming Model

Nextflow's core abstraction MUST adhere to the dataflow programming model:

- Workflows are defined as dataflow graphs where data flows between processes
- Processes MUST be stateless, side-effect-free transformations that communicate via channels
- The DSL MUST prioritize expressiveness for concurrent and parallel pipeline definition
- Changes to the language parser (ANTLR grammars in `nf-lang`) MUST preserve backward compatibility with existing pipelines unless explicitly versioned (DSL1 vs DSL2)
- Concurrency primitives (GPars actors/dataflow) MUST be used correctly to maintain the dataflow semantics

**Rationale**: The dataflow model is Nextflow's fundamental value proposition, enabling automatic parallelization and distribution. Preserving this model ensures existing scientific pipelines continue to work and users can reason about workflow behavior.

### IV. Apache 2.0 License Compliance

All source code MUST include Apache 2.0 license headers:

- Every source file MUST begin with the Apache 2.0 license header
- All contributions MUST comply with Apache 2.0 terms
- Third-party dependencies MUST use compatible licenses
- License compliance MUST be verified during code review

**Rationale**: Legal clarity protects both contributors and users. Consistent licensing enables academic and commercial use, which is critical for scientific software adoption.

### V. Developer Certificate of Origin (DCO) Sign-off

All commits MUST be signed with DCO certification:

- Contributors MUST certify they have the right to submit the code by using `git commit -s` or `git commit --signoff`
- Every commit message MUST include a `Signed-off-by` line
- The DCO bot MUST verify sign-off before any PR can be merged
- Contributors MUST NOT bypass the DCO requirement

**Rationale**: DCO provides legal protection and clear chain of custody for contributions, which is essential for open-source projects with diverse contributors.

### VI. Semantic Versioning and Release Discipline

Version management MUST follow strict semantic versioning with calendar-based releases:

- **Project versions** use calendar-based scheme: `YY.MM.PATCH` where April (`.04.`) and October (`.10.`) are stable releases, all other months use `-edge` suffix (e.g., `25.09.0-edge`)
- **Plugin versions** MUST use semantic versioning (`MAJOR.MINOR.PATCH`)
- Version changes MUST be documented in `changelog.txt` files (both project root and per-plugin)
- Breaking changes MUST increment MAJOR version for plugins and be clearly documented
- Release process MUST follow the documented procedure in `CLAUDE.md` including: updating changelogs, version files, running `make releaseInfo`, using `[release]` tag in commit message

**Rationale**: Predictable versioning enables users to understand compatibility and stability expectations. Calendar-based versioning for the main project makes release timing transparent, while semantic versioning for plugins enables clear communication of breaking changes.

### VII. Groovy Idioms and Code Standards

Code MUST follow Groovy best practices and Nextflow conventions:

- Use Groovy idioms (closures, operator overloading, DSL builders) appropriately
- Follow existing code patterns and conventions from similar modules
- Leverage Groovy's dynamic capabilities judiciously without sacrificing type safety where beneficial
- Use Groovy's `@CompileStatic` where performance is critical or type safety is desired
- AST transformations (in `modules/nextflow`) MUST be well-documented due to their compile-time magic
- Code MUST be formatted consistently (consider CodeNarc configuration in `gradle/codenarc.groovy`)

**Rationale**: Groovy enables powerful DSL capabilities that make Nextflow's language expressive, but requires discipline to maintain readability and debuggability. Consistency across the large codebase improves maintainability.

## Development Workflow

### Build and Development Process

- **Build tool**: Gradle with wrapper (`./gradlew`) is the authoritative build system
- **Quick commands**: Makefile provides convenience targets (`make compile`, `make test`, `make assemble`, `make check`, `make clean`)
- **Development testing**: Use `./launch.sh run script.nf` for testing changes against real workflows without full installation
- **Local installation**: `make install` publishes to Maven local for integration testing
- **Dependency management**: All dependencies MUST be declared in `build.gradle` with explicit versions; use `make deps` to analyze dependency trees

### Git Workflow

- **Branch management**: Work on feature branches, never commit directly to `master`
- **Commit sign-off**: Always use `git commit -s` to add DCO sign-off
- **CI control tags**: Use special commit message tags to control CI behavior:
  - `[ci skip]` - Skip CI tests entirely
  - `[ci fast]` - Run only unit tests, skip integration tests
  - `[e2e stage]` - Run end-to-end tests against Seqera platform staging environment
  - `[e2e prod]` - Run end-to-end tests against production platform
  - `[release]` - Trigger release automation
- **Pull requests**: Must pass all CI checks, require code review, and have DCO verification

### Architecture Decision Records (ADRs)

- Significant structural and technical decisions MUST be documented as ADRs in the `adr/` directory
- ADRs MUST follow the template format: date prefix + descriptive name (e.g., `20251114-module-system.md`)
- ADRs provide historical context for why architectural decisions were made
- When changing fundamental architecture, review existing ADRs and create new ones documenting the rationale

## Quality Standards

### Code Review Requirements

- All changes MUST go through pull request review
- Reviewers MUST verify:
  - Tests are included and passing
  - Code follows Groovy idioms and project conventions
  - License headers are present
  - DCO sign-off is present
  - Changes align with modular architecture principles
  - Breaking changes are appropriately versioned and documented

### Testing Gates

- `make test` MUST pass before committing locally
- All CI tests MUST pass before merging
- Integration tests MUST be run for changes affecting workflow execution
- Cloud validation tests MUST be run before releases touching cloud provider plugins
- Smoke tests enable rapid iteration but MUST NOT replace full test execution

### Performance and Compatibility

- Target platform: Java 17 runtime compatibility (development uses Java 21 toolchain)
- Performance-critical paths SHOULD be profiled and optimized
- Memory usage SHOULD be monitored for large-scale workflows
- Backward compatibility MUST be maintained for existing DSL features unless a new DSL version is introduced

## Governance

### Amendment Process

This constitution supersedes all other development practices. Amendments require:

1. **Proposal**: Submit amendment proposal via GitHub issue or pull request
2. **Discussion**: Community discussion period (minimum 1 week for major changes)
3. **Approval**: Approval from core maintainers
4. **Documentation**: Update this constitution with version bump following semantic versioning:
   - **MAJOR**: Backward incompatible governance changes, principle removal/redefinition
   - **MINOR**: New principle added or materially expanded guidance
   - **PATCH**: Clarifications, wording improvements, typo fixes

### Compliance and Review

- All pull requests and code reviews MUST verify compliance with these principles
- Deviations from principles MUST be explicitly justified in PR description
- Complexity additions MUST be justified against the "simplicity first" principle
- Constitution compliance is enforced through code review and CI automation where possible

### Ratification and Version History

**Version**: 1.0.0 | **Ratified**: 2025-11-17 | **Last Amended**: 2025-11-17

This constitution was derived from the Nextflow project's documented practices in `CLAUDE.md`, `CONTRIBUTING.md`, and the project's existing architectural patterns. It codifies the development principles that have made Nextflow a successful scientific workflow management system.
