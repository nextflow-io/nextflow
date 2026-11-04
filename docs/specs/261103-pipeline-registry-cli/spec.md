# Feature Specification: Pipeline Registry CLI and Remote Pipeline Inclusion

**Created**: 2026-11-03
**Status**: Draft
**ADR**: [20261103-pipeline-registry-cli](../../adr/20261103-pipeline-registry-cli.md)
**Server side**: seqeralabs/nextflow-registry `adr/20261104-pipeline-publish-and-download.md` (#440)

## Overview

This spec covers the Nextflow client side:

- Include a registry pipeline with `include { workflow as X } from 'scope/name'`, and install it automatically when it is missing.
- A `nextflow pipeline` command group: `install`, `list`, `remove`, `search` and `view` in M0, and `publish` in M1.
- Install from the Git commit of a release (M0) or from a registry bundle (M1).

The work is split into two milestones:

- **M0** uses only the registry API that exists today. It needs no server change.
- **M1** needs the server changes in the registry ADR (seqeralabs/nextflow-registry#440).

## Terms

- **Pipeline include**: an include whose source is `scope/name` and whose entries name `workflow` or `params`.
- **Installed pipeline**: a directory `pipelines/<scope>/<name>/` that holds a `.pipeline-info` file.
- **Pipeline root**: the nearest parent directory of a script that holds a `.pipeline-info` file. `scope/name` includes resolve against it.
- **Release source**: where the files of a release come from. It is `git` (repository + commit) or `artifact` (registry bundle + checksum).

---

## M0: Install and include from Git

### User Story 1: Include a registry pipeline (P1)

A meta-pipeline author includes `nf-core/rnaseq` by name and runs the meta-pipeline without installing anything by hand.

**Acceptance scenarios**

1. **Given** `main.nf` has `include { workflow as RNASEQ } from 'nf-core/rnaseq'` and `pipelines/nf-core/rnaseq/` is missing, **when** the user runs `nextflow run main.nf`, **then** Nextflow fetches the latest published release, installs it into `pipelines/nf-core/rnaseq/`, writes `.pipeline-info`, logs `Installed pipeline nf-core/rnaseq@<version>`, and runs.
2. **Given** the pipeline is installed and unchanged, **when** the user runs again, **then** no registry or Git call is made for that pipeline.
3. **Given** a file in the installed pipeline was edited by hand, **when** the user runs, **then** Nextflow warns that the pipeline was modified locally and uses the local files.
4. **Given** `include { FASTQC } from 'nf-core/fastqc'` (no `workflow` and no `params` entry), **then** it resolves to `modules/nf-core/fastqc/`, as today.
5. **Given** `include { params as RnaseqParams } from 'nf-core/rnaseq'` alone, **then** it is a pipeline include.
6. **Given** a registry name with no published release, or a revoked release, **then** compilation fails with an error that names the pipeline.

### User Story 2: An included pipeline resolves its own includes (P1)

An installed pipeline that uses registry-style includes behaves as it does when run directly.

**Acceptance scenarios**

1. **Given** `pipelines/nf-core/rnaseq/main.nf` has `include { STAR_ALIGN } from 'nf-core/star/align'`, **then** it resolves to `pipelines/nf-core/rnaseq/modules/nf-core/star/align/`. It does not resolve to the meta-pipeline's `modules/`.
2. **Given** that module is missing, **then** it is auto-installed under `pipelines/nf-core/rnaseq/modules/`.
3. **Given** `pipelines/nf-core/rnaseq/subworkflows/local/align.nf` has a `scope/name` include, **then** it resolves against `pipelines/nf-core/rnaseq/`, the pipeline root, and not against the directory of `align.nf`.
4. **Given** rnaseq includes pipeline `nf-core/other`, **then** it resolves to `pipelines/nf-core/rnaseq/pipelines/nf-core/other/`.
5. **Given** relative includes (`'./modules/...'`), **then** they keep resolving against the including file, with no change.
6. **Given** a script that is not under any pipeline root, **then** module resolution does not change.

### User Story 3: `nextflow pipeline install` (P1)

**Acceptance scenarios**

1. `nextflow pipeline install nf-core/rnaseq` installs the latest published release.
2. `-version 3.14.0` installs that release. It replaces an installed different version only if the installed copy is unmodified.
3. If the installed copy was edited by hand, install refuses and suggests `-force`. With `-force`, it replaces the copy.
4. If the installed version and checksum already match, it prints that the pipeline is up to date and does nothing.
5. The release's `repositoryUrl` is fetched at `commitId` through the existing repository providers (GitHub, GitLab, Bitbucket). `.git` is not kept.

### User Story 4: `list`, `remove`, `search`, `view` (P2)

**Acceptance scenarios**

1. `pipeline list` prints each installed pipeline under `pipelines/` (top level only): name, version, source, and status (`ok` or `modified`). It supports `-output json`.
2. `pipeline remove nf-core/rnaseq` deletes `pipelines/nf-core/rnaseq/`. It refuses a locally modified copy unless `-force` is given.
3. `pipeline search <query>` calls `GET /api/v1/pipelines?search=<query>` and prints name, latest version and description. It supports `-limit` and `-output json`, like `module search`.
4. `pipeline view nf-core/rnaseq [-version v]` prints name, version, description, keywords, repository, revision, commit and the registry page URL. It also prints a usage template: the include snippet, plus the params taken from the release `schema` when it has one. It supports `-output json`.

### M0 functional requirements

- **FR-001**: Path resolution (`ModuleResolver.getIncludeUri` and `RemoteModuleResolver`) receives the include entries and treats an include as a pipeline include when any entry names `workflow` or `params`.
- **FR-002**: The base directory for a `scope/name` include is the nearest of (a) the vendored module directory, using the existing rule, and (b) the pipeline root. If there is neither, the project directory is used. So a workflow module inside an installed pipeline still resolves its own nested `modules/`. *(Implementer's call, from the builder check.)*
- **FR-002b**: A module directory inside an installed pipeline that has no `.module-info` (for example, nf-core modules committed with `modules.json`) is used as it is. It is never reinstalled. The "no registry origin" warning is skipped there, because the pipeline's `.pipeline-info` checksum already covers those files. *(Implementer's call.)*
- **FR-003**: A pipeline include resolves to `<base>/pipelines/<scope>/<name>/main.nf`.
- **FR-004**: `.pipeline-info` stores `name`, `version`, the source (`type: git`, `repositoryUrl`, `revision`, `commitId`) and `checksum`, which is a SHA-256 over the installed files without `.pipeline-info`. It uses the same format and checksum code as `.module-info`.
- **FR-005**: Auto-install and `install` share one installer. The installer writes to a temporary directory and moves it into place only when it is complete, so that a failed fetch leaves no partial pipeline.
- **FR-006**: The registry client uses `registry.url` and `registry.auth`, as for modules. M0 calls the new `npr-client` methods `getPipeline`, `getPipelineRelease`, `listPipelineReleases` and `searchPipelines(query, limit)`, which map to `?search=&size=`. This needs an `npr-client` release after 0.29.1 (seqeralabs/nextflow-registry branch `261103-pipeline-client`). It is a client-library change only, with no change to server behaviour.
- **FR-007**: The language server must not install anything. It reports a missing pipeline like a missing module.

### M0 code locations (expected)

- `modules/nf-lang/.../ModuleResolver.java`, `ResolveIncludeVisitor.java`, `module/spi/RemoteModuleResolver.java`: pipeline include detection and the pipeline root.
- `modules/nextflow/.../module/`: the pipeline installer and `.pipeline-info`, reusing `ModuleInfo`, `ModuleChecksum` and `RegistryClientFactory`.
- `modules/nf-cli-v1/.../cli/CmdPipeline.groovy` and `cli/pipeline/CmdPipeline{Install,List,Remove,Search,View}.groovy`, modelled on `CmdModule*`.

---

## M1: Publish and install from bundle

Requires the server changes in the registry ADR (seqeralabs/nextflow-registry#440).

### User Story 5: `nextflow pipeline publish` (P1)

**Acceptance scenarios**

1. **Given** a pipeline folder with `manifest.name = 'myorg/mypipe'` and `manifest.version = '1.2.0'`, **when** the user runs `nextflow pipeline publish`, **then** the client builds a `.tar.gz` bundle and publishes release `1.2.0` of `myorg/mypipe`.
2. `-version` overrides `manifest.version`. A missing name or version, or a name that is not `scope/name`, is an error.
3. The bundle leaves out `.git` and the files listed in `.gitignore`.
4. If the folder is a Git checkout, the request carries `repositoryUrl` (from the `origin` remote), `gitRevision` (the tag at HEAD, or else the branch) and `commitId` (the HEAD SHA).
5. If the Git checkout has uncommitted changes, publish is refused unless `-force` is given.
6. If the folder is not a Git checkout, publish sends only the bundle.
7. A compressed bundle over 50 MB is refused by the client before upload.
8. `-dry-run` builds the bundle and prints name, version, size, checksum and provenance, without uploading.
9. Registry errors (401/403 for scope access, 409 when the version exists, 413 when too large) are shown with the server message.

### User Story 6: Install from bundle (P1)

**Acceptance scenarios**

1. **Given** the release response has `artifact`, **then** install downloads `artifact.url`, checks the `X-Checksum` header against `artifact.checksum`, and unpacks the bundle. No Git call is made.
2. **Given** the checksum does not match, **then** install fails and leaves nothing on disk.
3. **Given** the release has no `artifact`, **then** install uses Git, as in M0.
4. `.pipeline-info` records `type: artifact`, `url` and `checksum`, plus the Git provenance when the release has it.

### M1 functional requirements

- **FR-101**: Bundle creation reuses `ModuleStorage.createBundle`, extended to honour `.gitignore` and to skip `.git`.
- **FR-102**: Publish has two steps, through `npr-client`:
  1. `createPipelineDraft` sends `version`, `checksum` (`sha256:`), `size`, optional `repositoryUrl`, `gitRevision` and `commitId`, and metadata from the manifest (`description`, `homePage`, `keywords` when present).
  2. `uploadPipelineBundle` streams the bundle to `/releases/{releaseId}/upload`.
  If step 2 fails, the draft expires on the server and the client reports the error.
- **FR-103**: Publish requires a registry token (`registry.auth` or `NXF_REGISTRY_TOKEN`).

---

## Out of scope

`pipeline run`, `create`, `spec` and `validate`; a version in the include or in config; loading the included pipeline's config or plugins; `projectDir` for included pipelines; private Git repositories.

## Decisions left to the implementer

- The exact text and JSON shape of the `list`, `search` and `view` output, following the module commands.
- Whether the Git fetch uses a shallow clone or the host's archive download, provided that the result is the exact tree at `commitId`.
