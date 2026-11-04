# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-03
- Tags: pipelines, registry, cli, modules

## Summary

Add a `nextflow pipeline` command group and remote pipeline inclusion (`include { workflow as X } from 'nf-core/rnaseq'`). A missing pipeline is installed automatically from the Nextflow registry. M0 installs from the Git commit pinned by each release. M1 adds publishing from the CLI and installing from a bundle stored in the registry.

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path (`'./pipelines/nf-core/rnaseq'`). That ADR deferred remote inclusion and a `nextflow pipeline` CLI to future work.

The registry now has a read-only pipeline catalog (nextflow-registry #418). The catalog is filled from nf-core and from a web form (#419, #421). A release pins a Git repository, revision and commit. It does not store files. The registry has no public publish API: #421 removed it on purpose, so that the API can be designed together with the CLI.

## Goals

- **Remote pipeline inclusion**: `include { params as P ; workflow as X } from 'scope/name'` resolves a registry pipeline. Nextflow installs the pipeline if it is missing, as it does for modules.
- **Same model as modules**: vendored into the project, the version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own module includes resolve inside its own tree.
- **No hard dependency on Git**: a release can be installed from a registry bundle (M1), not only from a Git host.

## Non-goals

- `nextflow pipeline run`, `create`, `spec`, `validate`. Use `nextflow run` to run a pipeline.
- A version in the include statement or in config.
- Loading the included pipeline's config or plugins. The meta-pipeline provides them, as stated in the composition ADR.
- Changing `projectDir` for an included pipeline. See [Open questions](#open-questions).
- Private Git repositories. The registry accepts only public repositories.

## Considered Options

Where the installed files come from:

1. **Git only**: fetch `repositoryUrl` at `commitId`.
2. **Bundle only**: the registry stores a `.tar.gz`, as it does for modules.
3. **Both**: a release may carry a bundle. If it does not, the client falls back to Git.

## Pros and Cons of the Options

### Git only

- Good, because the registry already stores repo + commit for every release, including nf-core. It needs no server work.
- Good, because Nextflow already fetches repositories from GitHub, GitLab and Bitbucket (`AssetManager` and the repository providers).
- Bad, because install depends on the Git host, and on the repository still existing and being public.

### Bundle only

- Good, because it matches modules: one download, a checksum from the registry, no Git host involved.
- Bad, because it needs an upload API, storage, and packaging for every nf-core release before anything works.

### Both

- Good, because M0 can ship with no server change, and M1 removes the Git dependency for pipelines that are published as bundles.
- Bad, because the client has two install paths. The cost is small: the choice is made once, from the release response.

## Solution or decision outcome

**Both, in two phases.** M0 installs from Git and needs no server change. M1 adds a token-only publish endpoint that takes a bundle, and a download endpoint. When the release response has `artifact`, the client installs from it. When it does not, the client uses Git.

## Rationale & discussion

### Pipeline include detection

A `scope/name` include is a **pipeline include** when its entries name `workflow` or `params`. Module files have no entry workflow and no `params` block, so the two kinds can't be confused. The language already requires the alias `workflow as X`, so no new syntax is needed. A pipeline include resolves to `pipelines/<scope>/<name>/main.nf`. Any other `scope/name` include is a module include and resolves to `modules/<scope>/<name>/`, as today. The registry keeps modules and pipelines in separate namespaces, so `nf-core/rnaseq` may exist as both.

Today `ResolveIncludeVisitor` detects a pipeline include only after it has loaded the included file. The check moves earlier, into path resolution, so that the right directory is chosen and installed before loading.

### Storage layout

```
meta-pipeline/
├── main.nf                    # include { workflow as RNASEQ } from 'nf-core/rnaseq'
├── nextflow.config
└── pipelines/
    └── nf-core/
        └── rnaseq/
            ├── .pipeline-info
            ├── main.nf
            ├── nextflow.config
            ├── conf/modules.config
            ├── modules/       # the pipeline's own modules
            └── pipelines/     # pipelines that rnaseq itself includes
```

- The whole repository tree at the pinned commit is installed, without `.git`. This keeps `conf/` available for `includeConfig`, and keeps modules that nf-core commits in the repository.
- Installed pipelines should be committed to the meta-pipeline repository, like vendored modules.
- `.pipeline-info` records the version, the source (`git`: repository URL, revision and commit; or `artifact`: URL and checksum) and a checksum of the installed files.

### Pipeline root

The nearest parent directory that holds a `.pipeline-info` file is the **pipeline root**. A `scope/name` include resolves against the nearest of the vendored module directory (existing rule) and the pipeline root. If there is neither, it resolves against the project directory. Module includes go to `<base>/modules/`, and pipeline includes go to `<base>/pipelines/`. There is no search upward. Modules that a pipeline commits without `.module-info` (nf-core `modules.json` style) are used as they are and never reinstalled.

This is needed for correctness. Today `RemoteModuleResolver.resolveBaseDir` falls back to the project directory for any script that is not under `modules/`. So `include { STAR_ALIGN } from 'nf-core/star/align'` in an installed rnaseq would resolve to the meta-pipeline's `modules/`, which may hold a different version or nothing at all. Relative includes (`'./modules/...'`), which today's nf-core pipelines use, already resolve correctly and do not change.

### Versions and auto-install

The same rules as modules:

| Local state | Behaviour |
|---|---|
| Missing | Install the latest published release and log the version installed |
| Present, checksum matches `.pipeline-info` | Use it |
| Present, checksum differs (edited by hand) | Warn, use the local files, never overwrite |

`nextflow pipeline install <scope/name> -version <v>` pins a version. `-force` replaces locally edited files.

### CLI

| Command | Phase | Registry API |
|---|---|---|
| `pipeline install <scope/name> [-version] [-force]` | M0 | `GET /api/v1/pipelines/{name}[/{version}]` |
| `pipeline list` | M0 | none (reads `pipelines/`) |
| `pipeline remove <scope/name>` | M0 | none |
| `pipeline search <query>` | M0 | `GET /api/v1/pipelines?search=` |
| `pipeline view <scope/name> [-version] [-output json]` | M0 | `GET /api/v1/pipelines/{name}[/{version}]` |
| `pipeline publish [path] [-version] [-force] [-dry-run]` | M1 | two-step upload (new, see below) |

Each command follows its `nextflow module` counterpart, and all of them use `registry.url` and `registry.auth`.

### Publish (M1)

- Name and version come from `manifest.name` and `manifest.version`. `-version` overrides the version. The name must be `scope/name`, and the token must be allowed to write to that scope.
- The bundle is the pipeline folder without `.git` and without the files listed in `.gitignore`. It is built with the module bundle code (`ModuleStorage.createBundle`).
- **Git provenance is optional.** If the folder is a Git checkout, the client sends the repository URL, revision and commit. If the checkout has uncommitted changes, publish is refused unless `-force` is given, because otherwise the bundle would not match the commit it claims.
- The upload has **two steps**, as for plugins. First the client creates a draft release with the metadata, the bundle checksum and its size. Then it streams the bundle. Module publish sends the bundle base64-encoded inside JSON, but pipelines are too large for that. On 2026-11-03 the largest latest-release `.tar.gz` across the 106 active nf-core pipelines was 28 MB (scrnaseq 4.2.0), and the median raw tree was 2.8 MB. The limit is **50 MB** of compressed bundle. It is checked against the declared size before the upload.

### Registry contract (M1)

The server design is the registry ADR `adr/20261104-pipeline-publish-and-download.md` (seqeralabs/nextflow-registry#440). The client relies on these points:

- Publish: `POST /api/v1/pipelines/{name}/releases` (draft, JSON), then `POST /api/v1/pipelines/releases/{releaseId}/upload` (streamed multipart).

- The release response gains an optional `artifact: { url, checksum }`. Old clients ignore it.
- `GET /api/v1/pipelines/{name}/{version}/download` returns the bundle with an `X-Checksum` header, as module downloads do.
- `repositoryUrl` and `commitId` become optional, so that a pipeline can exist without Git.
- The current `contentHash` hashes the revision, commit and launch config. It is not a hash of the files, so it is not used for integrity.

## Open questions

- **`projectDir` in an included pipeline** still points to the meta-pipeline. Pipelines that read `${projectDir}/assets/...`, such as hlatyping and scrnaseq, break when included. The pipeline root could also become the included pipeline's `projectDir`, but that is a change to the language and runtime. It belongs to the composition ADR.
- **Typed pipelines**: inclusion requires `nextflow.enable.types = true`, which most nf-core pipelines don't enable yet.
- **Plugins** of an included pipeline must still be declared by the meta-pipeline (see the composition ADR).

## Links

- Implements: [Pipeline composition](20260608-pipeline-composition.md), section "Pipeline registry and CLI"
- Related: [Module system](20251114-module-system.md), [Workflow modules](20260608-workflow-modules.md)
- Spec: [docs/specs/261103-pipeline-registry-cli](../specs/261103-pipeline-registry-cli/spec.md)
- Registry: seqeralabs/nextflow-registry #418, #419, #421; M0 client methods #439; M1 server ADR #440
