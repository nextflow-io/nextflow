# Feature Specification: Pipeline Registry CLI and Remote Pipeline Inclusion

**Created**: 2026-11-03
**Status**: Draft
**ADR**: [20261103-pipeline-registry-cli](../../adr/20261103-pipeline-registry-cli.md)

## Overview

This spec covers the Nextflow client side:

- Include a registry pipeline with `include { workflow as X } from 'scope/name'`, and install it automatically when it is missing.
- A `nextflow pipeline` command group: `install`, `list`, `remove`, `search` and `view`.
- Install from the Git commit pinned by a registry release.
- `nextflow run scope/name` runs a registry release, and falls back to Git.

It uses only the registry API that exists today, so it needs no server change. The client needs the `npr-client` read methods from seqeralabs/nextflow-registry#439. Publishing and installing from a registry bundle are specified separately in `docs/specs/261105-pipeline-registry-publish/spec.md`.

## Terms

- **Pipeline include**: an include whose source is `scope/name` and whose entries name `workflow` or `params`.
- **Installed pipeline**: a directory `pipelines/<scope>/<name>/` that holds a `.pipeline-info` file.
- **Registry address**: a registry page URL, or `scope/name` with a registry host as `-hub` or `NXF_HUB` (`-hub registry.nextflow.io`). `nextflow run` resolves it with the registry into a repository URL and a commit.
- **Pipeline root**: the nearest parent directory of a script that holds a `.pipeline-info` file. `scope/name` includes resolve against it.

## User Story 1: Include a registry pipeline (P1)

A meta-pipeline author includes `nf-core/rnaseq` by name and runs the meta-pipeline without installing anything by hand.

**Acceptance scenarios**

1. **Given** `main.nf` has `include { workflow as RNASEQ } from 'nf-core/rnaseq'` and `pipelines/nf-core/rnaseq/` is missing, **when** the user runs `nextflow run main.nf`, **then** Nextflow fetches the latest published release, installs it into `pipelines/nf-core/rnaseq/`, writes `.pipeline-info`, logs `Installing pipeline nf-core/rnaseq@<version>...` and `Pipeline nf-core/rnaseq@<version> installed in <dir>`, and runs.
2. **Given** the pipeline is installed and unchanged, **when** the user runs again, **then** no registry or Git call is made for that pipeline.
3. **Given** a file in the installed pipeline was edited by hand, **when** the user runs, **then** Nextflow warns that the pipeline was modified locally and uses the local files.
4. **Given** a pipeline folder with no `.pipeline-info` (cloned by hand or a Git submodule), **when** the user runs, **then** Nextflow warns that it has no registry origin and uses it as it is.
5. **Given** `include { FASTQC } from 'nf-core/fastqc'` (no `workflow` and no `params` entry), **then** it resolves to `modules/nf-core/fastqc/`, as today.
6. **Given** `include { params as RnaseqParams } from 'nf-core/rnaseq'` alone, **then** it is a pipeline include.
7. **Given** a registry name with no published release, or a revoked release, **then** compilation fails with an error that names the pipeline.

## User Story 2: An included pipeline resolves its own includes (P1)

An installed pipeline that uses registry-style includes behaves as it does when run directly.

**Acceptance scenarios**

1. **Given** `pipelines/nf-core/rnaseq/main.nf` has `include { STAR_ALIGN } from 'nf-core/star/align'`, **then** it resolves to `pipelines/nf-core/rnaseq/modules/nf-core/star/align/`. It does not resolve to the meta-pipeline's `modules/`.
2. **Given** that module is missing, **then** it is auto-installed under `pipelines/nf-core/rnaseq/modules/`.
3. **Given** `pipelines/nf-core/rnaseq/subworkflows/local/align.nf` has a `scope/name` include, **then** it resolves against `pipelines/nf-core/rnaseq/`, the pipeline root, and not against the directory of `align.nf`.
4. **Given** rnaseq includes pipeline `nf-core/other`, **then** it resolves to `pipelines/nf-core/rnaseq/pipelines/nf-core/other/`.
5. **Given** relative includes (`'./modules/...'`), **then** they keep resolving against the including file, with no change.
6. **Given** a script that is not under any pipeline root, **then** module resolution does not change.

## User Story 3: `nextflow pipeline install` (P1)

**Acceptance scenarios**

1. `nextflow pipeline install nf-core/rnaseq` installs the latest published release.
2. `-version 3.14.0` installs that release. It replaces an installed different version only if the installed copy is unmodified.
3. If the installed copy was edited by hand, install refuses and suggests `-force`. With `-force`, it replaces the copy.
4. If the installed version and checksum already match, it prints `Pipeline <name>@<version> is already installed and up to date` and does nothing.
5. The release's `repositoryUrl` is fetched at `commitId` with Git (JGit). `.git` is not kept. If the revision no longer points to `commitId`, install fails. Credentials are not needed, because the registry accepts only public repositories.

## User Story 4: `list`, `remove`, `search`, `view` (P2)

The output follows the matching `nextflow module` command.

**Acceptance scenarios**

1. `pipeline list` prints a table of the installed pipelines under `pipelines/` (top level only), with the columns `Pipeline`, `Version` and `Status` (`ok` or `modified`) and no other text. `-output json` prints `{"pipelines": [ {name, version, status} ]}`.
2. `pipeline remove nf-core/rnaseq` deletes `pipelines/nf-core/rnaseq/`. It refuses a locally modified copy, or one with no `.pipeline-info`, unless `-force` is given.
3. `pipeline search <query>` calls `GET /api/v1/pipelines?search=<query>&page=<n>&size=<limit>`, reading pages until `-limit` results or the last page. It prints `Searching for '<query>'...`, then `Top <N> matching pipeline(s):`, then for each result the name, an indented `Version:` line and an indented `Description:` line. It supports `-limit` and `-output json`.
4. `pipeline view nf-core/rnaseq [-version v]` prints name, version, description, keywords, repository, revision, commit, the registry page URL, and the params taken from the release `schema` when it has one. It supports `-output json`. It prints the usage snippet `nextflow run <scope/name> -hub <registry host> -r <version>`, with the host of the primary registry.

## User Story 5: Run a registry pipeline (P1)

A user runs a pipeline they found in the registry, without knowing its Git repository.

**Acceptance scenarios**

1. **Given** `nf-core/rnaseq` is in the registry, **when** the user runs `nextflow run nf-core/rnaseq -hub registry.nextflow.io`, **then** Nextflow runs the latest published release at its `commitId`, and logs `Using registry release nf-core/rnaseq@<version> (commit <commit>)` before the launch line.
2. **Given** `https://registry.nextflow.io/pipelines/nf-core/rnaseq`, where the URL starts with a configured `registry.url` (without `/api`) or the default public registry, **then** the same release runs. It is resolved only against that registry, and the API key is sent only when it is the primary registry.
   A URL that ends with `@<version>` is an error that points to `-r`.
3. **Given** `NXF_HUB=registry-dev.nextflow.io` (in `registry.url`), **then** a plain `nextflow run jorge-test/nf-demo-fetchngs` uses that registry only, and `-hub github` uses Git. A `-hub` value that is neither a Git provider nor a registry host is an error that lists both.
4. **Given** `-r 3.14.0` and that release exists, **then** Nextflow runs it at its `commitId`.
5. **Given** a registry address and a name not in the registry, a `-r` that is not a release, a revoked release, or a registry that can't be reached, **then** the run fails with an error. There is no fallback to Git.
6. **Given** a plain `scope/name`, a Git URL or `-hub gitlab`, **then** the registry is not called, as today.
7. **Given** a registry name that points to a different Git repository (`jorge-test/nf-demo-fetchngs` → `github.com/jorgee/nf-demo-fetchngs`), **then** the release is cached under the registry name, `nextflow list` shows `jorge-test/nf-demo-fetchngs`, and `workflow.repository` is the Git URL.
8. **Given** a plain name whose cached clone comes from another Git project, **then** the clone is used, and Nextflow logs `Using <name> cached from <Git project>`.
9. **Given** a cached registry pipeline, **when** the user runs `-offline -hub registry-dev.nextflow.io -r <commit>`, or `pull`, `info` or `drop` by the registry name, **then** the cached clone and its own Git repository are used.
10. **Given** a clone cached under `acme/foo` from `github.com/other/foo`, **when** the user runs `https://github.com/acme/foo`, **then** the run fails because the name was downloaded from a different repository.
11. **Given** `nextflow inspect`, **then** the same addresses apply.
12. **Given** `nextflow pull` or `nextflow clone` with a registry page URL or a registry host as `-hub` or `NXF_HUB`, **then** the command fails with an error saying that registry pipelines are only supported by `nextflow run`.

## Functional requirements

- **FR-001**: Path resolution (`ModuleResolver.getIncludeUri` and `RemoteModuleResolver`) receives the include entries and treats an include as a pipeline include when any entry names `workflow` or `params`. The existing `getIncludeUri(URI, String, Path)` stays for the language server and resolves module includes only.
- **FR-002**: The base directory for a `scope/name` include is the nearest of (a) the vendored module directory, using the existing rule, and (b) the pipeline root. If there is neither, the project directory is used. So a workflow module inside an installed pipeline still resolves its own nested `modules/`.
- **FR-003**: A module directory inside an installed pipeline that has no `.module-info` (for example, nf-core modules committed with `modules.json`) is used as it is. It is never reinstalled. The "no registry origin" warning is skipped there, because the pipeline's `.pipeline-info` checksum already covers those files.
- **FR-004**: A pipeline include resolves to `<base>/pipelines/<scope>/<name>/<workflowPath>`, where `workflowPath` comes from the registry (default `main.nf`) and must stay inside the pipeline folder.
- **FR-005**: `.pipeline-info` is a Java properties file with `name`, `version`, `workflowPath`, `repositoryUrl`, `revision`, `commitId`, `registryUrl` and `checksum`. The checksum is a SHA-256 over the installed files without `.pipeline-info`, `.git/`, and subfolders that hold their own `.module-info` or `.pipeline-info`, so modules installed later inside the pipeline don't mark it modified.
- **FR-006**: Auto-install and `install` share one installer. The installer writes to a temporary directory and moves it into place only when it is complete, so that a failed fetch leaves no partial pipeline.
- **FR-007**: The registry client uses `registry.url` and `registry.auth`, as for modules. It calls the `npr-client` methods `getPipeline`, `getPipelineRelease`, `listPipelineReleases` and `listPipelines(search, page, size)` from seqeralabs/nextflow-registry#439. `search` reads further pages until `-limit` results, because the server caps the page size at 100.
- **FR-008**: The language server must not install anything. It reports a missing pipeline like a missing module.
- **FR-009**: `CmdRun` resolves a registry address with `PipelineManager.resolveRelease`, the same lookup as `pipeline install`, before it creates the `AssetManager`. The `AssetManager` gets the release's `repositoryUrl` as the Git address, its `commitId` as the revision, its `workflowPath` as the main script unless `-main-script` is given, and the registry name as the local name.
- **FR-010**: `AssetManager` keeps the local name (cache folder, `nextflow list`, `workflow.projectName`) apart from the remote project (Git provider). For a plain name, the remote project of a multi-revision clone is read from its `origin` remote. For a URL or a local name, a clone of a different repository is an error.
- **FR-011**: The lookup reads `registry.url` and `registry.auth` from `~/.nextflow/config`, the launch folder `nextflow.config` and `-c` files, as the `pipeline` commands do.
- **FR-012**: `PipelineCmdSupport.registryAddress`, used by `run`, `pull` and `clone`, tells whether a name and hub address a registry. A registry URL and a registry host as `-hub` share one host match against the configured `registry.url` values plus the default public registry; a URL must also start with the registry page URL (`PipelineCmdSupport.pageUrl`, the API URL without `/api`) followed by `/pipelines/`. A `-hub` value is checked as a Git provider name first. The lookup client holds only the matched registry, and gets the API key only when that registry is the primary one.

## Out of scope

Publishing and installing from a registry bundle (separate spec); registry addresses in `pull`, `config`, `clone`, `info` and `view`; a short prefix for registry addresses; `pipeline run`, `create`, `spec` and `validate`; a version in the include or in config; loading the included pipeline's config or plugins; `projectDir` for included pipelines; private Git repositories.
