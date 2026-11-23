# Pipeline registry publish and install from bundle

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-05
- Tags: pipelines, registry, cli, publish

## Summary

Add `nextflow pipeline publish`. It validates a pipeline folder and publishes it to the Nextflow registry as a `.tar.gz` bundle. The metadata comes from the manifest. A release with a bundle uses it for install, include, run and pull. A release without one (nf-core imports, web-form releases) keeps using its Git commit. This builds on [Pipeline registry CLI and remote pipeline inclusion](20261103-pipeline-registry-cli.md).

## Problem Statement

With remote pipeline inclusion, a release is installed from the Git commit it pins. So installing depends on the Git host, and on the repository still existing and being public. There is also no way to publish a pipeline from the CLI: the registry's catalog is filled by the nf-core import and a web form, and the web form can only create a pipeline's first release (seqeralabs/nextflow-registry#442 adds later Git releases from the UI).

## Goals

- **No hard dependency on Git**: a release can be installed from a bundle stored in the registry.
- **Publish from the CLI**: from a local pipeline folder, with a registry token. It also works in CI, with no prompts.
- **Publish only what can be run**: the folder is validated before anything is uploaded.
- **Bundle releases work everywhere**: `pipeline install`, includes, `nextflow run` and `nextflow pull`.
- **Keep Git-only releases working**: nf-core and web-form releases have no bundle and still use Git.

## Non-goals

- Building bundles for nf-core releases (registry side).
- Git provenance for CLI releases (see [Git provenance](#git-provenance)).
- A spec file for pipelines (like a module's `meta.yml`).
- Interactive questions during publish. Publish runs in CI too. It prints what it does. `-dry-run` shows the same without uploading.

## Considered Options

Where the installed files come from:

1. **Git only**: keep installing from the pinned commit.
2. **Bundle only**: every release must have a bundle.
3. **Both**: a release may carry a bundle; the client uses it when present and falls back to Git.

How the bundle is uploaded:

1. **One JSON request with the bundle base64-encoded inside**, as module publish does.
2. **Two steps, as for plugins**: a JSON draft, then a streamed multipart upload.

Whether a CLI release records Git provenance (repository URL, revision, commit):

1. **From the checkout**, when it is clean and its origin is a public URL.
2. **Never**: a CLI release is identified by its bundle and checksum only.

How `nextflow run` and `nextflow pull` cache a bundle release:

1. **A local Git repository built from the bundle**, used by `AssetManager` as a `file:` repository.
2. **A separate bundle folder**, `~/.nextflow/assets/.registry/<scope>/<name>/<version>/`, filled like `pipeline install`.

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

### Provenance: from the checkout

- Good, because the registry page links to the source.
- Bad, because the client can't make it reliable:
  - a tag or branch can move
  - a local commit may never have been pushed
  - a clone can have several remotes
- Bad, because the registry checks (public repository, commit exists) catch only some cases. They also block valid ones, such as private or unpushed repositories.

### Provenance: never

- Good, because the release is exactly its bundle, checked by its checksum, and nothing else can change. No clean-tree check, `-force` or `-no-git` is needed.
- Bad, because the registry page has no repository link. The manifest `homePage` can still point to the source.

### Run cache: local Git repository from the bundle

- Good, because most of the run-cache path is reused: offline runs, `list`, `drop` and `info`.
- Good, because each version can be a tag, so `-offline -r <version>` works for bundle releases.
- Bad, because Nextflow creates Git commits in a repository that exists only locally.
- Bad, because one registry name can have both Git releases and bundle releases. Git releases come from the web form, and later releases from seqeralabs/nextflow-registry#442. Bundle releases come from the CLI. `AssetManager` checks that a cached name keeps its provider and its repository, so mixing both needs changes there.
- Bad, because `pull -all` fetches from the clone's `origin`, which would be the local repository, so it never gets new versions from the registry.

### Run cache: separate bundle folder

- Good, because no Git is involved, and the extraction is the same as `pipeline install`.
- Good, because Git and bundle releases of the same name never share a cache entry.
- Bad, because `nextflow run` has to build the script file from that folder, outside `AssetManager`.
- Bad, because `list`, `drop` and `info` don't see these copies.
- Bad, because offline `-r <version>` needs its own lookup.

## Solution or decision outcome

**Both sources, a two-step upload, no provenance, and a local Git repository as the run cache for bundles.** A release may have `artifact: { url, checksum }`. When present, the client uses it for install, includes, run and pull. Otherwise it uses Git. Publishing validates the folder, creates a draft with the metadata, the bundle checksum and its size, then streams the bundle. A CLI release never records Git provenance. For `run` and `pull`, a bundle release is kept in a local Git repository built from the bundle. A prototype must show that it meets the requirements in [Run and pull from bundle](#run-and-pull-from-bundle); if it can't, the separate bundle folder is the fallback. The server design is the registry ADR `adr/20261104-pipeline-publish-and-download.md` (seqeralabs/nextflow-registry#440).

## Rationale & discussion

### Publish

- **Metadata comes from the manifest** in the pipeline's `nextflow.config`. The manifest is in the bundle, so the registry metadata must match it:
  - **name**: `manifest.name`, which must match the `scope/name` pattern. `-name` publishes under another registry name (for example a fork), with a warning that it differs from `manifest.name`.
  - **version**: `manifest.version`, semver as the registry requires. `-version` is accepted only when the manifest has no version; a different value is an error.
  - **main script**: `manifest.mainScript`, default `main.nf`. It is sent as the release's `workflowPath`, which install, includes and run use.
  - **description, homePage, keywords**: from the manifest.
- **Validation**, before anything is uploaded:
  - Errors, which stop the publish:
    - the main script is missing, does not parse with the same parser as `nextflow lint`, or has no entry workflow
    - the name or version is not valid, or conflicts with the options
    - the bundle is empty, or larger than 50 MB
  - Warnings:
    - no `README.md` or `nextflow_schema.json`, which the registry page shows
    - the pipeline does not enable static types, so it can be run but not included, as the [pipeline composition ADR](20260608-pipeline-composition.md) requires for included pipelines
  - Not checked: resolving includes or modules. That needs the network and belongs to `run` and `lint`.
- **Summary instead of questions**: publish prints what it will publish: the files, the size, and where each metadata field comes from. In a Git checkout, it also says that `.git` and the files Git ignores are excluded, and that no Git metadata is sent. `-dry-run` prints the same without uploading.
- **Bundle contents**: everything under the pipeline folder except `.git/`, files ignored by Git, `.pipeline-info`, and folders that hold their own `.module-info` (registry-installed modules). In a Git checkout, "ignored" follows Git's rules: nested `.gitignore` files and `.git/info/exclude`. Outside Git, the root `.gitignore` applies. Untracked files that are not ignored are included. Symlinks and submodules are skipped, with a warning.
- **Bundle format**: `.tar.gz` with POSIX long names, exec bits kept, a fixed timestamp and no owner names, so the same files always give the same `sha256:` checksum.
- **Upload failures**: a failed upload is retried once with the same draft. If it fails again, the version can be published again after the draft expires (1 hour).

### Git provenance

A CLI release records no repository, revision or commit. The client installs the bundle and checks its checksum. Provenance would only add a link, and that link can't be made reliable (see the options above). Git-only releases (nf-core imports, web form) keep their pinned commit, which the registry resolves itself.

### Install from bundle

- When the release has `artifact`, the client downloads the bundle, checks it against `artifact.checksum`, and extracts it, keeping `modules/` and exec bits and refusing paths that escape the folder.
- `.pipeline-info` records `source=artifact`, `artifactUrl` and `artifactChecksum`.
- Everything else stays as in remote pipeline inclusion: temporary folder, atomic move, checksum, local edits never overwritten.

### Run and pull from bundle

- For a release with `artifact`, `nextflow run` and `nextflow pull` download the bundle and check its checksum, as install does. They commit its files to a local Git repository built from the bundle, which `AssetManager` uses as the run cache.
- The cache must:
  - keep Git and bundle releases of the same registry name apart, or handle both in one entry
  - pin each version to the bundle that was checked, not to a name that can move
  - let `-offline -r <version>` find a version that was run or pulled before
  - let `pull -all` get new versions from the registry, not only from the local copy
  - define what `drop` removes
- A release without `artifact` runs from its Git commit, as in the pipeline registry CLI ADR.

### Registry contract

The client relies on:

- `POST /api/v1/pipelines/{name}/releases` (draft, JSON: version, checksum, size, workflowPath, description, homePage, keywords). The CLI does not send the optional repositoryUrl, gitRevision and commitId. Then `POST /api/v1/pipelines/releases/{releaseId}/upload` (multipart, part `payload`).
- `GET /api/v1/pipelines/{name}/{version}/download` with the `X-NF-Module-Checksum` header, as module downloads.
- The optional `artifact` field on releases; nullable `gitRevision` and `repositoryUrl`.
- The `npr-client` methods `createPipelineRelease`, `uploadPipelineRelease` and `downloadPipelineRelease`.

## Open questions

- **Run cache for bundle releases**: how the local Git repository meets the requirements in [Run and pull from bundle](#run-and-pull-from-bundle), in particular Git and bundle releases of the same name, and `pull -all`. To be settled by the prototype.

## Links

- Builds on: [Pipeline registry CLI and remote pipeline inclusion](20261103-pipeline-registry-cli.md)
- Registry: seqeralabs/nextflow-registry #440 (ADR and spec), #441 (implementation)
