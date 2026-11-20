# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque
- Status: accepted
- Date: 2026-11-03
- Tags: pipelines, registry, cli, modules

## Summary

Add a `nextflow pipeline` command group, remote pipeline inclusion (`include { workflow as X } from 'nf-core/rnaseq'`), and running a registry pipeline with `nextflow run` and an explicit registry address. Pipelines are installed from the Git commit pinned by the registry release.

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path (`'./pipelines/nf-core/rnaseq'`). That ADR deferred remote inclusion and a `nextflow pipeline` CLI to future work. There is also no way to run a registry pipeline: `nextflow run` only knows Git.

The registry now has a read-only pipeline catalog, filled from nf-core and from a web form. A release pins a Git repository, revision and commit. It does not store files.

## Goals

- **Remote pipeline inclusion**: `include { params as P ; workflow as X } from 'scope/name'` resolves a registry pipeline, installed automatically when missing.
- **Same model as modules**: vendored into the project, the version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Run a registry pipeline** with `nextflow run`, without changing what a plain `scope/name` means.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own module includes resolve inside its own tree.
- **No server change**: use only the registry API that exists today.

## Non-goals

- `nextflow pipeline publish` and installing from a registry bundle (see [Future work](#future-work)).
- `nextflow pipeline run`, `create`, `spec`, `validate`. `nextflow run` runs a pipeline.
- Registry addresses in `clone`, `config`, `info`, `view` and `drop`. These commands may be deprecated once the pipeline CLI exists. They work by registry name once the pipeline is cached, and `clone` rejects registry addresses.
- A bare `github.com/nf-core/rnaseq` address (without `https://`). It is not supported today.
- A version in the include statement or in config.
- Loading the included pipeline's config or plugins, and changing its `projectDir`. The composition ADR covers these.
- Private Git repositories. The registry accepts only public repositories.

## Considered Options

This ADR makes three decisions:

1. **Install source**: Git, or a bundle stored in the registry.
2. **Run addressing**: how `nextflow run` knows that a name refers to the registry. Registry first with a Git fallback, an explicit registry address, a reserved hub name, or a short prefix.
3. **Run cache**: where `run` keeps a registry pipeline. Under its Git name, under its registry name through `AssetManager`, or in a separate registry folder.

## Pros and Cons of the Options

### Install source

**Git**: fetch the release's `repositoryUrl` at its `commitId`.

- Good, because the registry already stores repo + commit for every release, including nf-core. It needs no server work.
- Good, because JGit is already a Nextflow dependency.
- Bad, because install depends on the Git host, and on the repository still existing and being public.

**Bundle stored in the registry**, as for modules.

- Good, because it matches modules: one download, a checksum from the registry, no Git host involved.
- Bad, because it needs an upload API, storage, and packaging for every release before anything works.

### Run addressing

**Registry first, Git fallback for a plain `scope/name`.**

- Bad, because the same name can resolve to different code depending on the registry answer, the network and the version, so a user can run an unintended pipeline without noticing.

**Explicit registry address**: a registry page URL, or a registry host as `-hub` / `NXF_HUB`.

- Good, because a plain `scope/name` keeps its meaning, so existing runs don't change.
- Good, because a page URL can be copied from the browser.
- Good, because a URL and a host follow one rule: an address names exactly one registry.

**A reserved hub name, `-hub registry`.**

- Bad, because it may clash with an SCM provider already called `registry`.
- Bad, because it would use all configured registries in order, while a URL names one.

**A short prefix**, such as a scheme or a leading character.

- Deferred, because `NXF_HUB=<registry host>` already gives a short form, and a prefix can be added later on top of the same lookup.

### Run cache

**Cache under the Git name.**

- Good, because it is the simplest: no `AssetManager` change.
- Bad, because the registry name is lost (`nextflow list` shows the Git name), and offline runs by registry name can't find the copy.

**Cache under the registry name, with a local name in `AssetManager`.**

- Good, because it is a small change, and `list`, `pull` and `drop` keep working.
- Bad, because when a registry name differs from its Git repository, a plain Git run of the registry name uses the cached clone without a check.

**A separate registry folder** (`~/.nextflow/assets/.registry/<scope>/<name>/<version>/`) with `.pipeline-info`.

- Good, because it also fits bundle releases, which have no Git.
- Bad, because `list`, `pull` and `drop` would not see it, and `run` would have to build the script path itself.

## Solution or decision outcome

- **Install source: Git now.** Install from the Git commit pinned by the release. The install code chooses the source in one place, so a bundle source can be added later.
- **Run addressing: an explicit registry address**, with no fallback to Git. A plain `scope/name` still runs from Git.
- **Run cache: the registry name**, through a local name in `AssetManager`. Bundle releases need a separate registry folder or a similar cache; their ADR decides.

## Rationale & discussion

### Pipeline includes

A `scope/name` include is a **pipeline include** when its entries name `workflow` or `params`. Module files have no entry workflow and no `params` block, so the two kinds can't be confused, and the language already requires the alias `workflow as X`, so no new syntax is needed. Path resolution makes this check before loading, so that the right directory is chosen and installed. The registry keeps modules and pipelines in separate namespaces, so `nf-core/rnaseq` may exist as both.

A pipeline include resolves to `pipelines/<scope>/<name>/`, with the main script given by the release. The whole repository tree at the pinned commit is installed, without `.git`. This keeps `conf/` available for `includeConfig`, and keeps the modules that nf-core commits in the repository. Installed pipelines should be committed to the meta-pipeline repository, like vendored modules.

```
meta-pipeline/
├── main.nf                    # include { workflow as RNASEQ } from 'nf-core/rnaseq'
└── pipelines/
    └── nf-core/
        └── rnaseq/            # pipeline root
            ├── .pipeline-info
            ├── main.nf
            ├── conf/
            ├── modules/       # the pipeline's own modules
            └── pipelines/     # pipelines that rnaseq itself includes
```

Versions, auto-install, local edits and `-force` follow the module rules. `.pipeline-info` records the release and a checksum of the installed files.

### Pipeline root

The nearest parent directory, up to the project directory, that holds a `.pipeline-info` file is the **pipeline root**. A `scope/name` include resolves against the nearest of the vendored module directory and the pipeline root, else the project directory, and is not looked up in an outer base.

This is needed for correctness. Before this change, a `scope/name` include in any script not under `modules/` resolved against the project directory. So `include { STAR_ALIGN } from 'nf-core/star/align'` in an installed rnaseq would resolve to the meta-pipeline's `modules/`, which may hold a different version or nothing at all. Relative includes (`'./modules/...'`), which today's nf-core pipelines use, already resolve correctly and do not change.

Modules that a pipeline commits without `.module-info` (nf-core `modules.json` style) are used as they are and never reinstalled, since the pipeline checksum covers them.

### CLI

`pipeline install`, `list`, `remove`, `search` and `view` follow their `nextflow module` counterparts, use `registry.url` and `registry.apiKey`, and call only existing read endpoints. They manage the vendored copies under the project's `pipelines/`, for `include`. The existing project commands (`pull`, `drop`, `list`, ...) keep managing the global run cache in `~/.nextflow/assets`.

### Running a registry pipeline

`nextflow run` and `nextflow pull` resolve a registry address to the release's repository and pinned commit. `pull` supports it so that a node with network can prepare offline runs.

- **Registry URL**: a configured `registry.url` without `/api`, followed by `/pipelines/<scope>/<name>`. `/pipelines/` keeps it apart from modules and other registry paths. The default public registry is always recognized, because the user typed its address. A URL on the host of a known Git provider stays a Git URL.
- **`-hub` / `NXF_HUB`**: a Git provider name is checked first, so existing setups behave the same, then a registry host. A hostname can't be confused with Git, because `-hub` only takes provider names.
- **One registry**: an address is resolved only against the registry it names, and the API key is sent only to the primary registry.
- **Version**: an address carries no version (`@<version>` is an error); the version is given with `-r`. The release runs at its pinned `commitId`, not at its tag, because a tag can be moved. A release without a commit ID runs at its revision.
- **No fallback**: with a registry address, every registry failure stops the run. The syntax is new, so no existing run depends on a fallback.

**Offline.** With `-offline`, the registry is not called. The run uses the copy cached under the registry name, and `-r` is used as a Git revision (a commit, or a tag equal to the version). This works for nf-core, whose versions are their Git tags, and fails for versions that are not tag names (a `v` prefix, or a release registered from a branch or commit). A local version-to-commit map (in the bare repository's Git config, a file next to the clone, local refs, or one central file) was not added. The Git config variant is the candidate if one is needed.

**Registry config.** `registry.url` and `registry.apiKey` are read as for the `pipeline` commands. The pipeline's own config can't be used, because it is not downloaded yet. The registry config is read only for a possible registry address, never for a plain name, a Git URL or a Git provider hub, because a launch config may be readable only once (Seqera Platform).

### Run cache

A registry release is cached under its **registry name**. `AssetManager` keeps two names: the local name (cache folder, `nextflow list`, `workflow.projectName`) and the remote Git project (Git provider, `workflow.repository`).

- A registry run is checked against the clone: a clone of a different repository under the registry name is an error, and offline the Git project is read from the clone.
- Plain Git names and URLs are not checked, as today. Comparing them is not reliable: some providers (Bitbucket Server) use a different path in the clone URL than in the project URL.

**Registry names mirror Git names.** A registry pipeline is expected to have the same `<scope>/<name>` as its Git repository, as nf-core pipelines do. Then the registry run and a plain Git run of the same name share one cache entry and run the same code. This is a registry best practice, not enforced by Nextflow.

## Future work

- **A short prefix** for registry addresses.
- **Publish and install from a bundle**: `nextflow pipeline publish`, and installing a release from a bundle stored in the registry so that install does not depend on a Git host. This is left to a separate ADR.

## Links

- Implements: [Pipeline composition](20260608-pipeline-composition.md), section "Pipeline registry and CLI"
- Related: [Module system](20251114-module-system.md), [Workflow modules](20260608-workflow-modules.md)
