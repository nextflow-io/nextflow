# Pipeline registry publish and install from bundle

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-05
- Tags: pipelines, registry, cli, publish

## Summary

Add `nextflow pipeline publish`, which publishes a pipeline folder to the Nextflow registry as a `.tar.gz` bundle with optional Git provenance. Let `nextflow pipeline install` and auto-install from an include install a release from its registry bundle when it has one, and fall back to the Git commit otherwise. This builds on [Pipeline registry CLI and remote pipeline inclusion](20261103-pipeline-registry-cli.md).

## Problem Statement

With remote pipeline inclusion, a release is installed from the Git commit it pins. So installing depends on the Git host, and on the repository still existing and being public. There is also no way to publish a pipeline from the CLI: the registry's catalog is filled by the nf-core import and a web form, and the web form can only create a pipeline's first release.

## Goals

- **No hard dependency on Git**: a release can be installed from a bundle stored in the registry.
- **Publish from the CLI**: from a local pipeline folder, with a registry token.
- **Optional provenance**: when the folder is a Git checkout, the release records where the code came from.
- **Keep Git-only releases working**: nf-core and web-form releases have no bundle and still install from Git.

## Non-goals

- Building bundles for nf-core releases (registry side).
- Private Git repositories.
- A spec file for pipelines (like a module's `meta.yml`).

## Considered Options

Where the installed files come from:

1. **Git only**: keep installing from the pinned commit.
2. **Bundle only**: every release must have a bundle.
3. **Both**: a release may carry a bundle; the client uses it when present and falls back to Git.

How the bundle is uploaded:

1. **One JSON request with the bundle base64-encoded inside**, as module publish does.
2. **Two steps, as for plugins**: a JSON draft, then a streamed multipart upload.

## Pros and Cons of the Options

### Source: Git only

- Good, because nothing new is needed.
- Bad, because install keeps depending on the Git host.

### Source: Bundle only

- Good, because it matches modules.
- Bad, because the nf-core import would have to build bundles for every release before anything works.

### Source: Both

- Good, because Git-only releases keep working, and bundles remove the Git dependency where they exist.
- Bad, because the client has two install paths. The cost is small: the source is chosen once, from the release response.

### Upload: one JSON request

- Good, because it is the module code path.
- Bad, because the whole bundle is held in server memory, base64-encoded. On 2026-11-03 the largest latest-release `.tar.gz` across the 106 active nf-core pipelines was 28 MB (scrnaseq 4.2.0). The registry moved plugins away from this for the same reason.

### Upload: two steps

- Good, because memory use stays flat, and the size and checksum are checked before and during the upload.
- Bad, because publishing takes two calls, and a failed upload leaves a draft that expires after 1 hour.

## Solution or decision outcome

**Both sources, and a two-step upload.** A release may have `artifact: { url, checksum }`. The client installs from it when present, and from Git otherwise. Publishing creates a draft with the metadata, the bundle checksum and its size, then streams the bundle. The server design is the registry ADR `adr/20261104-pipeline-publish-and-download.md` (seqeralabs/nextflow-registry#440).

## Rationale & discussion

### Publish

- **Name and version** come from `manifest.name` and `manifest.version` in the pipeline's `nextflow.config`. `-version` overrides the version. The name must match the same `scope/name` pattern as install. The version must be semver, as the registry requires.
- **Bundle contents**:
  - In a Git checkout, the files Git tracks under the pipeline folder, read from the working tree. A pipeline in a subfolder of a repo bundles only its own files.
  - Outside Git, every file except `.git/`, the paths matching the root `.gitignore`, and folders that hold their own `.module-info` (registry-installed modules).
  - Always without `.pipeline-info`. Symlinks and submodules are skipped, with a warning.
- **Bundle format**: `.tar.gz` with POSIX long names, exec bits kept, a fixed timestamp and no owner names, so the same files always give the same `sha256:` checksum.
- **Limits**: an empty bundle is refused. A bundle over 50 MB is refused before any upload.
- **Git provenance** (repository URL, tag at HEAD or else the branch, commit) is sent only when the folder is a Git checkout, the origin can be turned into a public `https://host/org/repo` URL, and the tree has no uncommitted changes. With uncommitted changes, publish is refused unless `-force` is given; then the pipeline is published as a bundle only, because its files don't match the commit.
- **Upload failures**: a failed upload is retried once with the same draft. If it fails again, the version can be published again after the draft expires (1 hour).

### Install from bundle

- When the release has `artifact`, the client downloads the bundle, checks it against `artifact.checksum`, and extracts it, keeping `modules/` and exec bits and refusing paths that escape the folder.
- `.pipeline-info` records `source=artifact`, `artifactUrl` and `artifactChecksum`, plus the Git provenance when the release has it.
- Everything else stays as in remote pipeline inclusion: temporary folder, atomic move, checksum, local edits never overwritten.

### Registry contract

The client relies on:

- `POST /api/v1/pipelines/{name}/releases` (draft, JSON: version, checksum, size, optional repositoryUrl/gitRevision/commitId, description, homePage, keywords), then `POST /api/v1/pipelines/releases/{releaseId}/upload` (multipart, part `payload`).
- `GET /api/v1/pipelines/{name}/{version}/download` with the `X-NF-Module-Checksum` header, as module downloads.
- The optional `artifact` field on releases; nullable `gitRevision` and `repositoryUrl`.
- The `npr-client` methods `createPipelineRelease`, `uploadPipelineRelease` and `downloadPipelineRelease`.

## Open questions

- **Bundle rule in a Git checkout.** Today it is "tracked files", which with `-force` drops new untracked files. The alternative is one rule everywhere: the working tree minus everything Git ignores, with Git used only for provenance.
- **Private or unpushed repositories.** When the origin is a public-host URL, the registry checks that the repository is public and that the commit exists there, so a private repo or an unpushed commit cannot be published. Options: a `-no-git` flag (bundle only), or documenting it with a "push first" hint.
- **New releases for pipelines registered from Git** with the web form (see the registry ADR).

## Links

- Builds on: [Pipeline registry CLI and remote pipeline inclusion](20261103-pipeline-registry-cli.md)
- Spec: [docs/specs/261105-pipeline-registry-publish](../specs/261105-pipeline-registry-publish/spec.md)
- Registry: seqeralabs/nextflow-registry #440 (ADR and spec), #441 (implementation)
