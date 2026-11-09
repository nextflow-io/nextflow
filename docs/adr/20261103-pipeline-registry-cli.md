# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-03
- Tags: pipelines, registry, cli, modules

## Summary

Add a `nextflow pipeline` command group, remote pipeline inclusion (`include { workflow as X } from 'nf-core/rnaseq'`), and registry lookup in `nextflow run nf-core/rnaseq`. A missing pipeline is installed automatically from the Nextflow registry, from the Git commit pinned by the registry release. Publishing from the CLI and installing from a bundle stored in the registry are covered by a separate ADR (see [Future work](#future-work)).

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path (`'./pipelines/nf-core/rnaseq'`). That ADR deferred remote inclusion and a `nextflow pipeline` CLI to future work.

The registry now has a read-only pipeline catalog (nextflow-registry #418). The catalog is filled from nf-core and from a web form (#419, #421). A release pins a Git repository, revision and commit. It does not store files.

## Goals

- **Remote pipeline inclusion**: `include { params as P ; workflow as X } from 'scope/name'` resolves a registry pipeline. Nextflow installs the pipeline if it is missing, as it does for modules.
- **Same model as modules**: vendored into the project, the version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Run a registry pipeline**: `nextflow run scope/name [-r version]` runs a registry release, and falls back to Git as today when the name or version is not in the registry.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own module includes resolve inside its own tree.
- **No server change**: use only the registry API that exists today. The client needs read methods in `npr-client` (nextflow-registry #439).

## Non-goals

- `nextflow pipeline publish` and installing from a registry bundle (see [Future work](#future-work)).
- `nextflow pipeline run`, `create`, `spec`, `validate`. `nextflow run` runs a pipeline.
- Registry lookup in `pull`, `config`, `clone`, `info`, `view` and `drop` (see [Future work](#future-work)).
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

Each command follows its `nextflow module` counterpart in options and output. All of them use `registry.url` and `registry.auth`. `view` shows the usage snippet `nextflow run <scope/name> -r <version>`.

The `pipeline` commands and the existing project commands have similar names, but they manage different things:

| Command | What it manages | Used for |
|---|---|---|
| `pull`, `drop`, `info`, `list`, `clone` | the global Git cache `~/.nextflow/assets` | `nextflow run` |
| `pipeline install`, `remove`, `list` | pipelines vendored into the project's `pipelines/`, with `.pipeline-info` | `include` |
| `view` | prints the cached `main.nf` | reading the code |
| `pipeline view` | registry metadata (versions, description, params) | choosing a release |

### Running a registry pipeline

There is no `nextflow pipeline run`. `nextflow run` looks the name up in the registry first and falls back to Git.

**Which names are looked up.** Only a name of the form `scope/name`, when:
- no local file or folder has that name (the local check already runs first)
- it is not a URL and does not end with `.nf`
- `-hub` is not given
- `-offline` is not given

A single word (`hello`), a URL (`https://github.com/nf-core/rnaseq`) and `-hub gitlab` never call the registry. A bare `github.com/nf-core/rnaseq` is not supported today and is not added:
- `https://` or `-hub` already skip the registry.
- It would be ambiguous with registry names. In module names the scope is the first segment and the name can have more than one segment (`nf-core/star/align`), so `github.com/nf-core/rnaseq` could also read as scope `github.com` and name `nf-core/rnaseq`.

**Version.** Without `-r`, Nextflow runs the latest published release. This changes behaviour: today it runs the default branch. For nf-core the result is the same, because their default branch is the latest release. With `-r`, the value is matched against the registry versions.

**Resolution.**

| Registry answer | What `nextflow run` does |
|---|---|
| Release found | Run it from the Git cache, at the release's `repositoryUrl` and `commitId` |
| Name not found (`404`) | Run from Git, as today. Logged at debug level only |
| `-r` is not a registry version | Warn, and run that Git revision of the release's `repositoryUrl` |
| Release revoked | Warn, and run from Git, as today |
| Registry down, timeout or `5xx` | Warn, and run from Git, as today |

The release runs at its pinned `commitId`, not at its tag, because a tag can be moved. The release's `workflowPath` is used as the main script unless `-main-script` is given.

**Warn, not fail.** nf-core pipelines are imported into the registry, so failing would break existing runs such as `nextflow run nf-core/rnaseq -r dev` or `-r <commit>`. Considered:
- **Always fail** for a revoked release or an unknown `-r`. Rejected: it breaks existing runs.
- **Fail by default, with an environment variable to opt out.** Rejected: it breaks existing runs until users learn the variable.
- **Warn by default, with an environment variable to opt in to failing.** Deferred: it can be added later without changing anything for existing users.

**Cache.** A registry release is fetched with `AssetManager`, with the release's `repositoryUrl` as the project and its `commitId` as the revision. It is cached in `~/.nextflow/assets/.repos/<project>/clones/<commit>/` like any Git revision. So `pull`, `list`, `info` and `drop` keep working, and `scm` credentials apply. A new registry cache was rejected: it would repeat what `AssetManager` does for Git, and the existing commands would ignore it.
- The cache uses the Git name. `jorge-test/nf-demo-fetchngs`, which points to `github.com/jorgee/nf-demo-fetchngs`, is listed as `jorgee/nf-demo-fetchngs`. For nf-core both names are the same.
- `-offline` skips the registry. So a pipeline whose registry name differs from its Git name must be run offline by its Git name.
- A bundle release with no Git (M1) can't use `AssetManager`. Its cache is decided in the M1 ADR.

**Registry config.** `registry.url` and `registry.auth` are read as for the `pipeline` commands: `~/.nextflow/config`, `nextflow.config` in the launch folder, and `-c` files. The pipeline's own config can't be used, because it is not downloaded yet.

**Output.** Before the launch line, Nextflow logs `Using registry release nf-core/rnaseq@3.14.0 (commit 1a2b3c4d5e)`. The launch line does not change: it shows the Git repository and commit, as for any Git run.

## Future work

- **Registry lookup in other commands**: `pull`, `config`, `clone`, `info`, `view` and `drop` use the Git name. With the lookup in `pull`, `-offline` would also work with a registry name that differs from the Git name.
- **Strict mode**: an opt-in environment variable to fail, instead of warning, on a revoked release or a `-r` that is not a registry version.
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
