# Feature Specification: Pipeline Registry Publish and Install from Bundle

**Created**: 2026-11-05
**Status**: Draft
**ADR**: [20261105-pipeline-registry-publish](../../adr/20261105-pipeline-registry-publish.md)
**Builds on**: [docs/specs/261103-pipeline-registry-cli](../261103-pipeline-registry-cli/spec.md)
**Server side**: seqeralabs/nextflow-registry #440 (ADR and spec), #441 (implementation)

## Overview

This spec covers the Nextflow client side of:

- `nextflow pipeline publish`: publish a pipeline folder as a bundle, with optional Git provenance.
- Installing a release from its registry bundle when it has one, and from Git otherwise.

## Terms

- **Bundle**: a `.tar.gz` of a pipeline folder, identified by its `sha256:` checksum.
- **Release source**: where the files of an installed release came from: `git` (repository + commit) or `artifact` (registry bundle + checksum).
- **Git provenance**: the repository URL, revision (tag or branch) and commit sent with a publish.

## User Story 1: `nextflow pipeline publish` (P1)

A pipeline author with a registry token publishes a version from a local folder.

**Acceptance scenarios**

1. **Given** a pipeline folder with `manifest.name = 'myorg/mypipe'` and `manifest.version = '1.2.0'`, **when** the user runs `nextflow pipeline publish`, **then** it prints `Publishing pipeline myorg/mypipe@1.2.0...`, creates a draft, uploads the bundle, and prints `Published pipeline myorg/mypipe@1.2.0`.
2. `-version` overrides `manifest.version`. A missing name or version is an error. A name that is not a valid `scope/name`, or a version that is not semver, is an error that says so.
3. **Given** a clean Git checkout whose origin is a public-host URL, **then** the request carries `repositoryUrl` (as `https://host/org/repo`), `gitRevision` (the tag at HEAD, else the branch) and `commitId`.
4. **Given** an origin that cannot be turned into a public `https` URL, or a folder that is not a Git checkout, **then** the pipeline is published as a bundle only, with a one-line note in the first case.
5. **Given** a Git checkout with uncommitted or untracked changes, **then** publish is refused and suggests `-force`. With `-force` it publishes as a bundle only, with a note.
6. **Given** a bundle with no files, **then** publish is refused with a hint that untracked files are not included in a Git checkout.
7. **Given** a bundle over 50 MB, **then** publish is refused before any registry call.
8. `-dry-run` prints the name, version, size, checksum and provenance, and makes no registry call.
9. **Given** the upload fails, **then** it is retried once with the same draft. If it fails again, the message says the version can be published again after the draft expires (1 hour).
10. Registry errors (401/403, 409 when the version exists, 413 when too large) are shown with the server's message.

## User Story 2: Install from bundle (P1)

**Acceptance scenarios**

1. **Given** the release has `artifact`, **when** it is installed (by `install` or by auto-install from an include), **then** the client downloads the bundle, checks it against `artifact.checksum`, and extracts it. No Git call is made.
2. **Given** the checksum does not match, **then** install fails with an error that names the pipeline, and nothing is left under `pipelines/`.
3. **Given** the release has no `artifact`, **then** install uses Git, as before.
4. `.pipeline-info` records `source=artifact`, `artifactUrl` and `artifactChecksum`, plus `repositoryUrl`, `revision` and `commitId` when the release has them. `pipeline list` shows the source.
5. Executable files (for example `bin/` scripts) stay executable after install.

## Functional requirements

- **FR-001**: Bundle contents. In a Git checkout (found with `findGitDir`, also for a pipeline in a repo subfolder): the files tracked in the Git index under the pipeline folder, read from the working tree. Otherwise: every file except `.git/`, paths matching the root `.gitignore`, and folders holding their own `.module-info` or `.pipeline-info`. Never `.pipeline-info`. Symlinks and submodules are skipped, with a warning.
- **FR-002**: Bundle format. `.tar.gz` with POSIX long-name support, mode `0755` for executables and `0644` otherwise (from the Git index in a checkout), a fixed timestamp and no owner names. Checksum `sha256:<hex>`.
- **FR-003**: Extraction keeps `modules/` and `tests/`, applies exec bits, and refuses entries that escape the target folder.
- **FR-004**: Publish has two steps through `npr-client`: `createPipelineRelease` with `version`, `checksum`, `size`, optional Git provenance, and `description`, `homePage` and `keywords` from the manifest; then `uploadPipelineRelease` with the bundle file. The temporary bundle is always deleted.
- **FR-005**: Install from bundle uses `downloadPipelineRelease` and compares the file's SHA-256 with `artifact.checksum` before extracting. A registry error becomes an error that names the pipeline.
- **FR-006**: Publish needs a token (`registry.auth` or `NXF_REGISTRY_TOKEN`). `-registry <url>` overrides `registry.url`.

## Out of scope

Bundles for nf-core releases; private Git repositories; a pipeline spec file; `pipeline run`.

## Open questions

See the ADR: the bundle rule in a Git checkout, publishing from private or unpushed repositories (`-no-git`), and new releases for Git-registered pipelines.
