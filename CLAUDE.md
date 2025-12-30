# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Nextflow is a scientific workflow management system built primarily in Groovy and Java.
It enables the creation of scalable, portable, and reproducible computational pipelines using a dataflow programming model.
The project follows a modular architecture with a plugin-based system for cloud providers and specialized features.

## Development Commands

### Build System
- **Primary build tool**: Gradle with wrapper (`./gradlew`)
- **Quick commands via Makefile**:
  - `make compile` - Compile sources and export classpath
  - `make assemble` - Build info, compile, and assemble
  - `make test` - Run all tests
  - `make check` - Run verification tasks
  - `make clean` - Clean build artifacts

### Testing
- **Unit tests**: `make test` or `./gradlew test`
  - Uses Spock Framework (Groovy-based testing)
  - JaCoCo integration for code coverage
- **Specific test class**: `./gradlew test --tests "SomeTestClass"`
- **Specific module**: `make test module=nextflow`
- **Smoke tests**: `make smoke` or `NXF_SMOKE=1 ./gradlew test` (skips long-running tests, network-dependent tests, and cloud provider integration tests)
- **Integration tests**: Run from `tests/` directory using test runner scripts
- **Cloud validation tests**: Located in `validation/` directory (requires credentials)

### Development Workflow
- **Development launcher**: `./launch.sh run script.nf` (uses development build)
- **Dependency analysis**: `make deps` or `make deps config=runtime`
- **Install locally**: `make install` (installs to Maven local)

## Architecture

### Core Modules (`modules/`)
- **nextflow**: Main application module with core workflow engine, CLI, AST transformations, executors
- **nf-commons**: Shared utilities, plugin system infrastructure, extension methods
- **nf-httpfs**: HTTP filesystem support and custom providers
- **nf-lang**: Language parsing, ANTLR grammars, AST implementation
- **nf-lineage**: Data lineage tracking and workflow execution history

### Plugin System (`plugins/`)
- **Cloud providers**: nf-amazon (AWS), nf-azure (Azure), nf-google (GCP)
- **Execution platforms**: nf-k8s (Kubernetes)
- **Services**: nf-tower (Seqera Platform), nf-wave (container management)
- **Other**: nf-console (interactive interface), nf-cloudcache (cloud caching)

### Key Technologies
- **Language**: Groovy 4.0.29 (Java-compatible, targeting Java 17)
- **Concurrency**: GPars 1.2.1 (Actor model, parallel/concurrent programming)
- **Build**: Gradle with Java 21 toolchain
- **Parsing**: ANTLR for Nextflow DSL
- **Serialization**: Kryo
- **Database**: LevelDB for local caching
- **Version Control**: JGit integration

### Testing Structure
- **Unit tests**: Each module has `src/test/groovy/` with Spock Framework tests
- **Integration tests**: `tests/` directory with .nf workflows and expected outputs
- **Legacy tests**: `tests-v1/` for DSL v1 compatibility
- **Validation tests**: `validation/` directory for cloud provider end-to-end testing
- **Documentation tests**: `docs/snippets/` for verifying documentation examples

## Development Notes

### Code Standards
- All code must include Apache 2.0 license headers
- Contributions require Developer Certificate of Origin (DCO) sign-off
- Use existing code patterns and conventions from similar modules
- Follow Groovy idioms and leverage the Nextflow DSL patterns

### Common Development Tasks
- **Local development**: Use `make compile` to build and `./launch.sh` to test changes
- **Adding features**: First check modules like `nextflow` for core features or create plugins for specialized functionality
- **Plugin development**: Follow existing plugin patterns in `plugins/` directory
- **Testing changes**: Always run `make test` before committing
- **Cloud testing**: Use validation scripts in `validation/` directory with appropriate credentials

### Build Configuration
- Java toolchain uses version 21 for development, targets Java 17 compatibility
- Uses shadow plugin for creating fat JARs
- Maven publication to S3-based Seqera repositories
- Multi-module project with shared dependencies managed in root build.gradle

### Git conventions

- **DCO sign-off required**: All commits must be signed by adding a `Signed-off-by` line to the commit message or by using the `-s` option (see CONTRIBUTING.md for details).
- **Always use sign-off**: Use `git commit -s` or `git commit --signoff` for commits to avoid DCO bot issues
- **CI control tags**: Use special tags in commit messages to control CI behavior:
  - `[ci skip]` - Skip the execution of CI tests
  - `[ci fast]` - Run only unit tests and skip integration tests
  - `[e2e stage]` - Run end-to-end tests vs Seqera platform stage environment
  - `[e2e prod]` - Same but against production platform
  - `[release]` - Trigger release process

## Important Files
- `VERSION`: Define the current version number
- `nextflow`: Launch wrapper script (updated by build process)
- `.launch.classpath`: Development classpath (generated by `make compile`)
- `build.gradle`: Root build configuration with multi-module setup
- `settings.gradle`: Gradle project structure definition
- `plugins/*/VERSION`: Define the version of the corresponding plugin sub-project.
- `adr/`: Architecture Decision Records (ADRs) documenting significant structural and technical decisions in the project

## Release process

Follow these actions to make a new release:

- Update the `changelog.txt` file in each plugin sub-project (if any change has been done).
- Update the `VERSION` file in in each plugin sub-project.
  Use a semantic version number depending the impact of the change, or do not change
  if no changes have been done to the plugin.
- Update `nextflowVersion` attribute in the `build.gradle` file for plugins requiring specific
  Nextflow versions.
- Commit the version and changelog files changes independently for each plugin. Use as commit
  message the template `Bump plugin-name@version` e.g. `Bump nf-amazon@2.0.0.
- Update `VERSION` file in the project root using a calendar-like versioning scheme. Versions in the 4-th and 10-th month are "stable releases", e.g. `25.10.0`, while versions in all other months are "edge releases", e.g. `25.09.0-edge`.
- Update the project root `changelog.txt` with changes since the past release. Use the git log
  command to determine what changed e.g. `git log v<PREVIOUS VERSION>..`
- Run `make releaseInfo` to update the version number and generate checksums.
- Run this command to stage for commit the release files:
    ```
    git add \
      VERSION \
      changelog.txt \
      nextflow \
      nextflow.md5 \
      nextflow.sha1 \
      nextflow.sha256 \
      modules/nextflow/src/main/resources/META-INF/plugins-info.txt \
      modules/nextflow/src/main/resources/META-INF/build-info.properties
    ```
- Make a commit using the `[release]` and `[e2e prod]` tags in the comment and push it upstream to trigger the release automation with GitHub action:
    ```
    git commit -m "[release] Nextflow version 25.09.0-edge"
    git push origin master
    ```

## Active Technologies
- Groovy 4.0.29, Java 17 target (Java 21 toolchain) + AWS SDK v2 (ECS, EC2, CloudWatch Logs, S3), GPars 1.2.1 (001-ecs-executor)
- AWS S3 via Seqera Fusion filesystem (001-ecs-executor)

## Recent Changes
- 001-ecs-executor: Added Groovy 4.0.29, Java 17 target (Java 21 toolchain) + AWS SDK v2 (ECS, EC2, CloudWatch Logs, S3), GPars 1.2.1
