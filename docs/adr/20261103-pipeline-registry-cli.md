# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque, Ben Sherman
- Status: accepted
- Date: 2026-11-25
- Tags: pipelines, registry, cli, modules

## Summary

Support registry pipelines in two ways, distinguished by use rather than by kind. To **include** a pipeline, the `nextflow module` commands install it into the project's `modules/` like a module, and `include { workflow as X } from 'nf-core/rnaseq'` resolves it like any registry include. To **run** a pipeline, the existing project commands (`run`, `pull`, `clone`, `info`, `view`, `drop`) accept an explicit registry address. Pipelines are installed from the Git commit pinned by the registry release.

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path. That ADR deferred remote inclusion and a pipeline CLI to future work. There is also no way to run a registry pipeline: `nextflow run` only knows Git.

The registry has a read-only pipeline catalog, filled from nf-core and from a web form. A release pins a Git repository, revision and commit. It does not store files.

## Goals

- **Remote pipeline inclusion**: `include { workflow as X } from 'scope/name'` resolves a registry pipeline, installed automatically when missing.
- **Same model as modules**: vendored into the project, version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Run a registry pipeline** with `nextflow run`, without changing what a plain `scope/name` means.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own registry includes resolve inside its own tree.
- **Small CLI surface**: no new command group for users to learn.

## Non-goals

- Publishing pipelines, and installing from a registry bundle (see [Future work](#future-work)).
- Running a pipeline with `nextflow module run`. `nextflow run` runs pipelines.
- The registry as the default source for plain names in the project commands.
- A version in the include statement or in config.
- Loading the included pipeline's config or plugins, and changing its `projectDir`. The composition ADR covers these.
- Private Git repositories. The registry accepts only public repositories.

## Principle: distinguish by use, not by kind

A registry name can be used in two ways, and each already has a home:

- **Include** is project-level. Included code is vendored into `modules/` and committed with the project. The `nextflow module` commands manage it.
- **Run** is user-level. Run code is cached in `$NXF_ASSETS`. The project commands (`run`, `pull`, `list`, `info`, `view`, `drop`, `clone`) manage it.

Whether the name is a module or a pipeline matters only when it is downloaded: a module comes from a registry bundle, a pipeline from Git. Everywhere else, the kind is an implementation detail.

## Considered Options

1. CLI surface
2. Install folder
3. Kind detection
4. Info file
5. Include base directory
6. Install source
7. Run addressing
8. Run cache

## Pros and Cons of the Options

### CLI surface

**A `nextflow pipeline` command group** (`install`, `list`, `remove`, `search`, `view`), mirroring `nextflow module`.

- Good, because each command handles one kind.
- Bad, because it duplicates five commands, and users must know the kind of a name before they can install it.
- Bad, because it overlaps with the project commands, which already manage pipelines in the run cache. Two commands would list "pipelines" in different places.

**The `module` commands cover both kinds; the project commands run.**

- Good, because the command follows from what the user wants to do, which they always know, and not from the kind, which they may not.
- Good, because no command is added. `module install`, `list`, `remove`, `search` and `view` gain pipelines; `run`, `pull`, `clone`, `info`, `view` and `drop` gain registry addresses.
- Bad, because the `module` name is less accurate. `module run` and `-update-deps` stay module-only.

### Install folder

**`pipelines/<scope>/<name>`.**

- Good, because the kind is visible in the tree.
- Bad, because an include must know the kind to pick the folder, which forces kind detection before loading (see below).
- Bad, because a second vendoring folder needs its own resolution and nesting rules.

**`modules/<scope>/<name>`.**

- Good, because `include ... from 'scope/name'` always resolves to `modules/<scope>/<name>/main.nf`, whatever the kind.
- Good, because the nesting rules of [workflow modules](20260608-workflow-modules.md) apply unchanged.

### Kind detection

**From the include entries**: an include whose entries are `workflow` or `params` is a pipeline include.

- Bad, because it rules out including other entries of a pipeline, such as a record type.
- Bad, because the same name could resolve to different folders depending on what is included.

**Only when downloading**: try the module lookup first, and the pipeline lookup on a 404.

- Good, because includes and every other command ignore the kind.
- Good, because a module install costs no extra call. A pipeline install costs one failed module lookup.
- Bad, because it requires modules and pipelines to share one registry name space (see [Prerequisite](#prerequisite-one-registry-name-space)).

### Info file

**`.pipeline-info`** next to the pipeline.

- Bad, because `module list`, `remove` and the resolver would look for two files.

**`.module-info` with `kind=Pipeline`.**

- Good, because one file marks every installed folder, and a file without `kind` stays a module, so existing installs are unaffected.

### Include base directory

Before this change, a `scope/name` include in any script not under `modules/` resolved against the project directory. An installed pipeline keeps its subworkflows outside `modules/` (`subworkflows/`, `workflows/`), so `include { STAR_ALIGN } from 'nf-core/star/align'` in an installed rnaseq would resolve to the calling project's `modules/`, which may hold a different version or nothing.

**A pipeline-specific root**, found from `.pipeline-info`.

- Bad, because it is a rule for one kind.

**The nearest folder of the including file, up to the project directory, that holds `.module-info`**, takes the place of the project directory; otherwise the previous rule applies.

- Good, because one rule covers installed pipelines and any installed module, and a module's own `main.nf` resolves as before.
- Good, because hand-vendored folders without `.module-info` keep the previous behavior.

### Install source

**Git**: fetch the release's repository at its pinned commit.

- Good, because the registry already stores repository and commit for every release, including nf-core. It needs no server work.
- Bad, because install depends on the Git host, and on the repository still existing and being public.

**Bundle stored in the registry**, as for modules.

- Good, because it matches modules: one download, no Git host involved.
- Bad, because it needs an upload API, storage, and packaging for every release before anything works.

### Run addressing

**Registry first, Git fallback for a plain `scope/name`.**

- Bad, because the same name can resolve to different code depending on the registry answer, the network and the version, so a user can run an unintended pipeline without noticing.

**Explicit registry address**: a registry page URL, or a registry host as `-hub` / `NXF_HUB`.

- Good, because a plain `scope/name` keeps its meaning, so existing runs don't change.
- Good, because a page URL can be copied from the browser.
- Good, because a URL and a host follow one rule: an address names exactly one registry.

**A reserved hub name, `-hub registry`.**

- Bad, because it may clash with an SCM provider already called `registry`, and it would use all configured registries while a URL names one.

**A short prefix**, such as a scheme or a leading character.

- Deferred, because `NXF_HUB=<registry host>` already gives a short form, and a prefix can be added later on top of the same lookup.

### Run cache

**Cache under the Git name.**

- Bad, because the registry name is lost (`nextflow list` shows the Git name), and offline runs by registry name can't find the copy.

**Cache under the registry name, with a local name in `AssetManager`.**

- Good, because it is a small change, and `list`, `pull`, `info`, `view` and `drop` keep working.
- Bad, because when a registry name differs from its Git repository, a plain Git run of the registry name uses the cached clone without a check. Registry names are expected to mirror Git names, as nf-core pipelines do, so that both runs share one cache entry and run the same code. This is a registry best practice, not enforced by Nextflow.

**A separate registry folder** in the assets directory.

- Bad, because the project commands would not see it, and `run` would have to build the script path itself.

## Solution or decision outcome

- **CLI**: no `pipeline` command group. The `module` commands install, list, remove, search and view both kinds; `run`, `pull`, `clone`, `info`, `view` and `drop` accept registry addresses.
- **Install folder**: `modules/<scope>/<name>`, holding the whole Git tree at the pinned commit without `.git`, so `conf/` stays available for `includeConfig`.
- **Kind**: found only when downloading, module first, pipeline on a 404.
- **Info file**: `.module-info`, with `kind=Pipeline` for pipelines. A pipeline has no `meta.yml`, so its version is read from `.module-info`, and its checksum covers the whole tree, excluding `.git` and nested installed folders.
- **Base directory**: the nearest folder holding `.module-info`.
- **Install source**: Git now; a registry bundle is future work.
- **Run addressing**: an explicit registry address, with no fallback to Git.
- **Run cache**: the registry name. `AssetManager` keeps two names: the local name (cache folder, `nextflow list`, `workflow.projectName`) and the remote Git project (`workflow.repository`). A clone of a different repository under the registry name is an error.

## Rationale & discussion

### Prerequisite: one registry name space

Finding the kind on download requires that a name is either a module or a pipeline, never both. The registry kept them apart, so the name space is unified in nextflow-registry, as a separate change. The only known nf-core collision is `isoseq`, which is being raised with nf-core.

### Installing and including a pipeline

An include always resolves `main.nf`, so a pipeline without a top-level `main.nf` cannot be installed: it would leave a folder that fails on every later use. If `main.nf` exists but the release declares another main script, `module install` warns. Versions, auto-install, local edits and `-force` follow the module rules. `module run` on a pipeline aborts before downloading anything and points to `nextflow run`. `module remove -keep-files` is rejected for a pipeline, because without `.module-info` the folder can no longer be included.

Modules that a pipeline commits without `.module-info` (nf-core `modules.json` style) are used as they are and never reinstalled, since the pipeline checksum covers them. Relative includes, which nf-core pipelines use today, already resolve correctly and do not change.

### Running a registry pipeline

`run`, `pull` and `clone` resolve a registry address to the release's repository and pinned commit. `pull` supports it so that a node with network can prepare offline runs; `clone` keeps `.git` at the pinned commit, with the full history unless `-deep` is given, as for a Git project.

- **Registry URL**: a configured `registry.url` without `/api`, followed by `/pipelines/<scope>/<name>`. The default public registry is always recognized, because the user typed its address.
- **`-hub` / `NXF_HUB`**: a Git provider name is checked first, so existing setups behave the same, then a registry host.
- **One registry**: an address is resolved only against the registry it names, and the API key is sent only to the primary registry.
- **Version**: an address carries no version; the version is given with `-r`. The release runs at its pinned commit, not at its tag, because a tag can be moved.
- **No fallback**: with a registry address, every registry failure stops the run. The syntax is new, so no existing run depends on a fallback.
- **Registry config**: read from the user config, the launch directory config and `-c` files, never from the pipeline config, which is not downloaded yet. It is read only for a possible registry address, because a launch config may be readable only once (Seqera Platform).

**Offline and cache-only commands.** With `-offline`, the registry is not called: the run uses the copy cached under the registry name, and `-r` is a Git revision (a commit, or a tag equal to the version). This works for nf-core, whose versions are their Git tags. `info`, `view` and `drop` only read the cache, so they never call the registry and their `-r` is always a Git revision. A local version-to-commit map was not added.

## Future work

- **A short prefix** for registry addresses.
- **Publish and install from a bundle**, so that install does not depend on a Git host. This is left to a separate ADR.
- **An intent-based name** for the `module` commands, now that they also manage pipelines.

## Links

- Implements: [Pipeline composition](20260608-pipeline-composition.md), section "Pipeline registry and CLI"
- Related: [Module system](20251114-module-system.md), [Workflow modules](20260608-workflow-modules.md)
