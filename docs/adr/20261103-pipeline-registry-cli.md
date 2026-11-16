# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-03
- Tags: pipelines, registry, cli, modules

## Summary

Add a `nextflow pipeline` command group, remote pipeline inclusion (`include { workflow as X } from 'nf-core/rnaseq'`), and running a registry pipeline with `nextflow run` and an explicit registry address. A missing pipeline is installed automatically from the configured registries, from the Git commit pinned by the registry release. Publishing from the CLI and installing from a bundle stored in the registry are left to a separate ADR (see [Future work](#future-work)).

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path (`'./pipelines/nf-core/rnaseq'`). That ADR deferred remote inclusion and a `nextflow pipeline` CLI to future work. There is also no way to run a registry pipeline: `nextflow run` only knows Git.

The registry now has a read-only pipeline catalog, filled from nf-core and from a web form. A release pins a Git repository, revision and commit. It does not store files.

## Goals

- **Remote pipeline inclusion**: `include { params as P ; workflow as X } from 'scope/name'` resolves a registry pipeline. Nextflow installs the pipeline if it is missing, as it does for modules.
- **Same model as modules**: vendored into the project, the version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Run a registry pipeline**: `nextflow run` runs a registry release when the registry is addressed explicitly, with a registry page URL or with a registry host as `-hub` (`-hub registry.nextflow.io`). A plain `scope/name` keeps going to Git, as today.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own module includes resolve inside its own tree.
- **No server change**: use only the registry API that exists today, through new read methods in `npr-client`.

## Non-goals

- `nextflow pipeline publish` and installing from a registry bundle (see [Future work](#future-work)).
- `nextflow pipeline run`, `create`, `spec`, `validate`. `nextflow run` runs a pipeline.
- Registry addresses in `config`, `info`, `view` and `drop`. They work by registry name once the pipeline is cached. `clone` rejects registry addresses.
- A version in the include statement or in config.
- Loading the included pipeline's config or plugins, and changing its `projectDir`. The composition ADR covers these.
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

**Git now.** Install from the Git commit pinned by the release. This works with the registry as it is. Vendored installs use a small JGit fetcher (public repositories, no credentials), and `run` and `pull` use `AssetManager` and its cache. The install code chooses the source in one place, so a bundle source can be added later (see [Future work](#future-work)).

A release without a commit ID is fetched at its revision, which can move.

## Rationale & discussion

### Pipeline include detection

A `scope/name` include is a **pipeline include** when its entries name `workflow` or `params`. Module files have no entry workflow and no `params` block, so the two kinds can't be confused. The language already requires the alias `workflow as X`, so no new syntax is needed. A pipeline include resolves to `pipelines/<scope>/<name>/<workflowPath>`, where `workflowPath` comes from the release (default `main.nf`) and must stay inside the pipeline directory. Any other `scope/name` include is a module include and resolves to `modules/<scope>/<name>/`, as today. The registry keeps modules and pipelines in separate namespaces, so `nf-core/rnaseq` may exist as both.

Path resolution detects a pipeline include by its entry names, before loading, so that the right directory is chosen and installed. The check after loading stays, for the typing rules.

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
- `.pipeline-info` (Java properties) records the name, version, main script path, registry URL, Git repository, revision and commit, and a checksum of the installed files.

### Pipeline root

The nearest parent directory, up to the project directory, that holds a `.pipeline-info` file is the **pipeline root**. A `scope/name` include resolves against the nearest of the vendored module directory (existing rule) and the pipeline root. If there is neither, it resolves against the project directory. Module includes go to `<base>/modules/`, and pipeline includes go to `<base>/pipelines/`. An include that is not found under its base is not looked up in an outer base.

Modules that a pipeline commits without `.module-info` (nf-core `modules.json` style) are used as they are and never reinstalled. The pipeline checksum covers them, and skips subdirectories with their own `.module-info` or `.pipeline-info`, so modules installed later inside a pipeline don't mark it as modified.

This is needed for correctness. Before this change, `RemoteModuleResolver.resolveBaseDir` fell back to the project directory for any script that is not under `modules/`. So `include { STAR_ALIGN } from 'nf-core/star/align'` in an installed rnaseq would resolve to the meta-pipeline's `modules/`, which may hold a different version or nothing at all. Relative includes (`'./modules/...'`), which today's nf-core pipelines use, already resolve correctly and do not change.

### Versions and auto-install

The same rules as modules:

| Local state | Behaviour |
|---|---|
| Missing | Install the latest published release and log the version installed |
| Present, checksum matches `.pipeline-info` | Use it |
| Present, checksum differs (edited by hand) | Warn, use the local files, never overwrite |
| Present, no `.pipeline-info` (cloned or submodule) | Warn, use it as it is |

`nextflow pipeline install <scope/name> -version <v>` pins a version. An unmodified copy of another version is replaced. `-force` is needed to replace or remove a copy that is edited or has no `.pipeline-info`.

The checksum of an installed pipeline is computed once per process, and a failing include (including its auto-install) is tried once per compile.

### CLI

`pipeline install`, `list`, `remove`, `search` and `view` follow their `nextflow module` counterparts in options and output, use `registry.url` and `registry.auth`, and call only existing read endpoints. `view` shows the usage snippet `nextflow run <scope/name> -hub <primary registry host> -r <version>`.

The `pipeline` commands manage vendored copies under the project's `pipelines/`, for `include`. The existing project commands (`pull`, `drop`, `list`, ...) keep managing the global run cache in `~/.nextflow/assets`.

### Running a registry pipeline

There is no `nextflow pipeline run`. `nextflow run` runs a registry release only when the registry is addressed explicitly. A plain `scope/name` keeps its meaning: a Git project on the default hub.

**Addresses.**

| Address | Example | Source |
|---|---|---|
| Registry page URL | `nextflow run https://registry.nextflow.io/pipelines/nf-core/rnaseq` | that registry |
| Registry host as hub | `nextflow run nf-core/rnaseq -hub registry.nextflow.io` | that registry |
| Registry host as default hub | `NXF_HUB=registry.nextflow.io nextflow run nf-core/rnaseq` | that registry |
| Plain name, Git URL or Git hub | `nextflow run nf-core/rnaseq`, `https://github.com/nf-core/rnaseq`, `-hub gitlab` | Git, as today |

- A registry URL is the registry page address: a configured `registry.url` without `/api`, followed by `/pipelines/<scope>/<name>`, ignoring the scheme. It can be copied from the browser, and `/pipelines/` keeps it apart from modules (`nf-core/rnaseq` may be both) and from other registry paths. The default public registry is always recognized, because the user typed its address.
- A `-hub` (or `NXF_HUB`) value is read in this order: a Git provider name (Git, as today), then the host (or `host:port`) of a configured or the public registry, then an error that lists both. Provider names come first, so existing setups behave the same. A hostname can't be confused with Git, because `-hub` only takes provider names.
- A URL on the host of a known Git provider is a Git URL, even when its path has a `pipelines` segment.
- An address names exactly one registry, and is resolved only against it. The API key is sent only to the primary registry in `registry.url`.
- An address carries no version: `@<version>` is an error, also with `-offline`. The version is given with `-r`.
- `pull` supports registry addresses like `run`, so that a node with network can prepare offline runs. `clone` rejects them. The other project commands don't get registry support, since they may be deprecated once the pipeline CLI exists.
- A bare `github.com/nf-core/rnaseq` (without `https://`) is not added. It is not supported today.

**Version.** Without `-r`, the latest published release runs. With `-r`, the value must be a registry version. The release runs at its pinned `commitId`, not at its tag, because a tag can be moved. The release's `workflowPath` is the main script unless `-main-script` is given.

**No fallback.** With an explicit registry address, every failure stops the run with an error: name not found, `-r` not a release, release revoked, registry not reachable. The syntax is new, so no existing run depends on a fallback.

**Offline.** With `-offline`, the registry is not called. The run uses the copy cached under the registry name, and `-r` is used as a Git revision (a commit, or a tag equal to the version). This works for nf-core, whose versions are their Git tags. It fails for versions that are not tag names (a `v` prefix, or a release registered from a branch or commit), and the error suggests the commit. A local version-to-commit map (in the bare repository's Git config, a file next to the clone, local refs, or one central file) was not added. The Git config variant is the candidate if one is needed.

**Registry config.** `registry.url` and `registry.auth` are read as for the `pipeline` commands: `~/.nextflow/config`, `nextflow.config` in the launch folder, and `-c` files. The pipeline's own config can't be used, because it is not downloaded yet. The registry config is read only for a possible registry address, never for a plain name, a Git URL or a Git provider hub, so Git runs read their config once, as before. A launch config may be readable only once (Seqera Platform).

**Considered addressing options.**
- **Registry first, Git fallback for a plain `scope/name`.** Rejected: the same name can resolve to different code depending on the registry answer, the network and the version, so a user can run an unintended pipeline without noticing.
- **Explicit registry address (chosen).** A registry URL and a registry host as `-hub`, as above.
- **A reserved hub name, `-hub registry`.** Rejected: it may clash with an SCM provider already called `registry`, and it would use all configured registries in order, while a URL names one. A host keeps one rule for both: an address names exactly one registry.
- **A short prefix**, such as a scheme or a leading character. Deferred: `NXF_HUB=<registry host>` already gives a short form, and a prefix can be added later on top of the same lookup.

### Run cache

A registry release is fetched with `AssetManager` from its `repositoryUrl` at its `commitId`, and cached under its **registry name**: `~/.nextflow/assets/.repos/<scope>/<name>/`.

- `AssetManager` keeps two names. The project name is the local name: the cache folder, `nextflow list`, and `workflow.projectName`. The remote project is the Git one, used only by the Git provider (API calls and clone URL), and so by `workflow.repository`.
- For a registry address, a clone of a different repository under the registry name is an error ("already downloaded from a different repository").
- For `run -offline` with a registry address, the remote project is read from the clone's `origin` remote (multi-revision clones only), and Nextflow logs `Using <name> cached from <Git project>` when it differs from the registry name.
- Plain Git names and URLs are not checked against the clone, as today. Comparing them is not reliable: some providers (Bitbucket Server) use a different path in the clone URL than in the project URL.

**Registry names mirror Git names.** A registry pipeline is expected to have the same `<scope>/<name>` as its Git repository, as nf-core pipelines do. Then the registry run and a plain Git run of the same name share one cache entry and run the same code. When a registry name differs from its Git repository, the registry run still works, but a plain Git run of the registry name uses the cached clone, which comes from another repository, without a check. This is a registry best practice, not enforced by Nextflow.

**Considered cache options.**
- **Cache under the Git name.** Simplest: no `AssetManager` change. Rejected: the registry name is lost (`nextflow list` shows the Git name), and offline runs by registry name can't find the copy.
- **`AssetManager` with a local name (chosen).** A small change to `AssetManager`. The existing commands keep working.
- **A separate registry folder** (`~/.nextflow/assets/.registry/<scope>/<name>/<version>/`) managed by `PipelineManager`, with `.pipeline-info`. It also fits bundle releases, which have no Git. Not chosen for now: `list`, `pull` and `drop` would not see it, and `CmdRun` would build the script file itself. Bundle releases need it or a similar cache; their ADR decides.

## Future work

- **A short prefix** for registry addresses.
- **Publish and install from a bundle**: `nextflow pipeline publish`, and installing a release from a bundle stored in the registry so that install does not depend on a Git host. This is left to a separate ADR.

## Links

- Implements: [Pipeline composition](20260608-pipeline-composition.md), section "Pipeline registry and CLI"
- Related: [Module system](20251114-module-system.md), [Workflow modules](20260608-workflow-modules.md)
- Spec: [docs/specs/261103-pipeline-registry-cli](../specs/261103-pipeline-registry-cli/spec.md)
