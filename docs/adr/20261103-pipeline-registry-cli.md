# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-03
- Tags: pipelines, registry, cli, modules

## Summary

Add a `nextflow pipeline` command group and remote pipeline inclusion (`include { workflow as X } from 'nf-core/rnaseq'`). A missing pipeline is installed automatically from the Nextflow registry, from the Git commit pinned by the registry release. Publishing from the CLI and installing from a bundle stored in the registry are covered by a separate ADR (see [Future work](#future-work)).

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path (`'./pipelines/nf-core/rnaseq'`). That ADR deferred remote inclusion and a `nextflow pipeline` CLI to future work.

The registry now has a read-only pipeline catalog (nextflow-registry #418). The catalog is filled from nf-core and from a web form (#419, #421). A release pins a Git repository, revision and commit. It does not store files.

## Goals

- **Remote pipeline inclusion**: `include { params as P ; workflow as X } from 'scope/name'` resolves a registry pipeline. Nextflow installs the pipeline if it is missing, as it does for modules.
- **Same model as modules**: vendored into the project, the version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own module includes resolve inside its own tree.
- **No server change**: use only the registry API that exists today. The client needs read methods in `npr-client` (nextflow-registry #439).

## Non-goals

- `nextflow pipeline publish` and installing from a registry bundle (see [Future work](#future-work)).
- `nextflow pipeline run`, `create`, `spec`, `validate`. Use `nextflow run` to run a pipeline.
- A version in the include statement or in config.
- Loading the included pipeline's config or plugins. The meta-pipeline provides them, as stated in the composition ADR.
- Changing `projectDir` for an included pipeline. See [Open questions](#open-questions).
- Private Git repositories. The registry accepts only public repositories.

## Considered Options

Where the installed files come from:

1. **Git**: fetch the release's `repositoryUrl` at its `commitId`.
2. **A bundle stored in the registry**, as for modules.

## Pros and Cons of the Options

### Git

- Good, because the registry already stores repo + commit for every release, including nf-core. It needs no server work.
- Good, because JGit is already a Nextflow dependency.
- Bad, because install depends on the Git host, and on the repository still existing and being public.

### Bundle stored in the registry

- Good, because it matches modules: one download, a checksum from the registry, no Git host involved.
- Bad, because it needs an upload API, storage, and packaging for every release before anything works.

## Solution or decision outcome

**Git now.** Install from the Git commit pinned by the release. This works with the registry as it is. The install code chooses the source in one place, so a bundle source can be added later without changing anything else (see [Future work](#future-work)).

## Rationale & discussion

### Pipeline include detection

A `scope/name` include is a **pipeline include** when its entries name `workflow` or `params`. Module files have no entry workflow and no `params` block, so the two kinds can't be confused. The language already requires the alias `workflow as X`, so no new syntax is needed. A pipeline include resolves to `pipelines/<scope>/<name>/main.nf`. Any other `scope/name` include is a module include and resolves to `modules/<scope>/<name>/`, as today. The registry keeps modules and pipelines in separate namespaces, so `nf-core/rnaseq` may exist as both.

Before this change, `ResolveIncludeVisitor` detected a pipeline include only after it had loaded the included file. The check moves earlier, into path resolution, so that the right directory is chosen and installed before loading.

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
- `.pipeline-info` records the version, the repository URL, revision and commit, and a checksum of the installed files.

### Pipeline root

The nearest parent directory that holds a `.pipeline-info` file is the **pipeline root**. A `scope/name` include resolves against the nearest of the vendored module directory (existing rule) and the pipeline root. If there is neither, it resolves against the project directory. Module includes go to `<base>/modules/`, and pipeline includes go to `<base>/pipelines/`. There is no search upward. Modules that a pipeline commits without `.module-info` (nf-core `modules.json` style) are used as they are and never reinstalled.

This is needed for correctness. Before this change, `RemoteModuleResolver.resolveBaseDir` fell back to the project directory for any script that is not under `modules/`. So `include { STAR_ALIGN } from 'nf-core/star/align'` in an installed rnaseq would resolve to the meta-pipeline's `modules/`, which may hold a different version or nothing at all. Relative includes (`'./modules/...'`), which today's nf-core pipelines use, already resolve correctly and do not change.

### Versions and auto-install

The same rules as modules:

| Local state | Behaviour |
|---|---|
| Missing | Install the latest published release and log the version installed |
| Present, checksum matches `.pipeline-info` | Use it |
| Present, checksum differs (edited by hand) | Warn, use the local files, never overwrite |
| Present, no `.pipeline-info` (cloned or submodule) | Warn, use it as it is |

`nextflow pipeline install <scope/name> -version <v>` pins a version. `-force` replaces locally edited files.

### CLI

| Command | Registry API |
|---|---|
| `pipeline install <scope/name> [-version] [-force]` | `GET /api/v1/pipelines/{name}[/{version}]` |
| `pipeline list [-output json]` | none (reads `pipelines/`) |
| `pipeline remove <scope/name> [-force]` | none |
| `pipeline search <query> [-limit] [-output json]` | `GET /api/v1/pipelines?search=&page=&size=` |
| `pipeline view <scope/name> [-version] [-output json]` | `GET /api/v1/pipelines/{name}[/{version}]` |

Each command follows its `nextflow module` counterpart in options and output. All of them use `registry.url` and `registry.auth`. `view` does not show a usage snippet yet; that waits until we decide between running a registry pipeline with `nextflow run` or with `nextflow pipeline run`.

## Future work

- **Publish and install from a bundle**: `nextflow pipeline publish`, and installing a release from a bundle stored in the registry so that install does not depend on a Git host. This is decided in a separate ADR, `20261105-pipeline-registry-publish.md`, together with the server design in nextflow-registry #440.

## Open questions

- **`projectDir` in an included pipeline** still points to the meta-pipeline. Pipelines that read `${projectDir}/assets/...`, such as hlatyping and scrnaseq, break when included. The pipeline root could also become the included pipeline's `projectDir`, but that is a change to the language and runtime. It belongs to the composition ADR.
- **Typed pipelines**: inclusion requires `nextflow.enable.types = true`, which most nf-core pipelines don't enable yet.
- **Plugins** of an included pipeline must still be declared by the meta-pipeline (see the composition ADR).

## Links

- Implements: [Pipeline composition](20260608-pipeline-composition.md), section "Pipeline registry and CLI"
- Related: [Module system](20251114-module-system.md), [Workflow modules](20260608-workflow-modules.md)
- Spec: [docs/specs/261103-pipeline-registry-cli](../specs/261103-pipeline-registry-cli/spec.md)
- Registry: seqeralabs/nextflow-registry #418, #419, #421; client read methods #439
