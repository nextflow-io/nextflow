# Pipeline registry CLI and remote pipeline inclusion

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-11-03
- Tags: pipelines, registry, cli, modules

## Summary

Add a `nextflow pipeline` command group, remote pipeline inclusion (`include { workflow as X } from 'nf-core/rnaseq'`), and running a registry pipeline with `nextflow run` and an explicit registry address. A missing pipeline is installed automatically from the Nextflow registry, from the Git commit pinned by the registry release. Publishing from the CLI and installing from a bundle stored in the registry are covered by a separate ADR (see [Future work](#future-work)).

## Problem Statement

[Pipeline composition](20260608-pipeline-composition.md) lets a pipeline include another pipeline and call it like a named workflow. Today the included pipeline must be vendored by hand, by cloning or with Git submodules, and included by relative path (`'./pipelines/nf-core/rnaseq'`). That ADR deferred remote inclusion and a `nextflow pipeline` CLI to future work.

The registry now has a read-only pipeline catalog (nextflow-registry #418). The catalog is filled from nf-core and from a web form (#419, #421). A release pins a Git repository, revision and commit. It does not store files.

## Goals

- **Remote pipeline inclusion**: `include { params as P ; workflow as X } from 'scope/name'` resolves a registry pipeline. Nextflow installs the pipeline if it is missing, as it does for modules.
- **Same model as modules**: vendored into the project, the version read from the installed files, integrity checked against a local checksum, local edits never overwritten.
- **Run a registry pipeline**: `nextflow run` runs a registry release when the registry is addressed explicitly, with a registry page URL or with a registry host as `-hub` (`-hub registry.nextflow.io`). A plain `scope/name` keeps going to Git, as today.
- **Reproducibility**: an included pipeline behaves as it does when run directly. Its own module includes resolve inside its own tree.
- **No server change**: use only the registry API that exists today. The client needs read methods in `npr-client` (nextflow-registry #439).

## Non-goals

- `nextflow pipeline publish` and installing from a registry bundle (see [Future work](#future-work)).
- `nextflow pipeline run`, `create`, `spec`, `validate`. `nextflow run` runs a pipeline.
- Resolving registry addresses in `pull`, `config`, `clone`, `info`, `view` and `drop` (see [Future work](#future-work)). They already work on a registry pipeline once it is cached, by its registry name.
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

Each command follows its `nextflow module` counterpart in options and output. All of them use `registry.url` and `registry.auth`. `view` shows the usage snippet `nextflow run <scope/name> -hub <registry host> -r <version>`.

The `pipeline` commands and the existing project commands have similar names, but they manage different things:

| Command | What it manages | Used for |
|---|---|---|
| `pull`, `drop`, `info`, `list`, `clone` | the global Git cache `~/.nextflow/assets` | `nextflow run` |
| `pipeline install`, `remove`, `list` | pipelines vendored into the project's `pipelines/`, with `.pipeline-info` | `include` |
| `view` | prints the cached `main.nf` | reading the code |
| `pipeline view` | registry metadata (versions, description, params) | choosing a release |

### Running a registry pipeline

There is no `nextflow pipeline run`. `nextflow run` runs a registry release only when the registry is addressed explicitly. A plain `scope/name` keeps its meaning: a Git project on the default hub.

**Addresses.**

| Address | Example | Source |
|---|---|---|
| Registry page URL | `nextflow run https://registry.nextflow.io/pipelines/nf-core/rnaseq` | that registry |
| Registry host as hub | `nextflow run nf-core/rnaseq -hub registry.nextflow.io` | that registry |
| Registry host as default hub | `NXF_HUB=registry-dev.nextflow.io nextflow run jorge-test/nf-demo-fetchngs` | that registry |
| Plain name | `nextflow run nf-core/rnaseq` | Git, default hub (`NXF_HUB` or `github`), as today |
| Git URL or Git hub | `https://github.com/nf-core/rnaseq`, `-hub gitlab` | Git, as today |

- A URL is a registry address when it starts with one of the configured `registry.url` values (without `/api`) followed by `/pipelines/`. This is the same rule as for Git URLs, where the host picks the provider. The default public registry (`https://registry.nextflow.io`) is always recognized, even when `registry.url` lists only other registries, because the user typed its address.
- An address names exactly one registry, whether it is a URL or a host. It is resolved only against that registry, not against the other configured ones. The API key is sent only when that registry is the primary one in `registry.url`, since the key belongs to the primary registry.
- A URL carries no version. A URL that ends with `@<version>` is an error: the version is given with `-r`.
- A `-hub` (or `NXF_HUB`) value is read in this order:
  1. a Git provider name (`github`, `gitlab`, a name from the SCM file): Git, as today. This is checked first, so existing setups behave the same.
  2. the host of a configured `registry.url`, or of the default public registry: that registry.
  3. anything else: an error that lists the known providers and registry hosts.
- A hostname can't be confused with Git, because `-hub` only takes provider names. Through `NXF_HUB`, a user or a site can make a registry the default for plain names.
- A bare `github.com/nf-core/rnaseq` (without `https://`) is not added. It is not supported today, and it would be ambiguous with registry names, whose name part can have more than one segment (`nf-core/star/align`).

**Version.** Without `-r`, the latest published release runs. With `-r`, the value must be a registry version. The release runs at its pinned `commitId`, not at its tag, because a tag can be moved. The release's `workflowPath` is the main script unless `-main-script` is given.

**No fallback.** With an explicit registry address, every failure stops the run with an error: name not found, `-r` not a release, release revoked, registry not reachable. The syntax is new, so no existing run depends on a fallback.

**Offline.** With `-offline`, the registry is not called. The run uses the copy cached under the registry name, and `-r` is used as a Git revision (a commit, or a tag equal to the version).

**Registry config.** `registry.url` and `registry.auth` are read as for the `pipeline` commands: `~/.nextflow/config`, `nextflow.config` in the launch folder, and `-c` files. The pipeline's own config can't be used, because it is not downloaded yet.

**Output.** Before the launch line, Nextflow logs `Using registry release nf-core/rnaseq@3.14.0 (commit 1a2b3c4d5e)`.

**Considered addressing options.**
- **Registry first, Git fallback for a plain `scope/name`.** Rejected: the same name can resolve to different code depending on the registry answer, the network and the version, so a user can run an unintended pipeline without noticing.
- **Explicit registry address (chosen).** A registry URL and a registry host as `-hub`, as above.
- **A reserved hub name, `-hub registry`.** Rejected: it may clash with an SCM provider already called `registry`, and it would use all configured registries in order, while a URL names one. A host keeps one rule for both: an address names exactly one registry.
- **A short prefix**, such as a scheme or a leading character. Deferred: `NXF_HUB=<registry host>` already gives a short form, and a prefix can be added later on top of the same lookup.

### Run cache

A registry release is fetched with `AssetManager` from its `repositoryUrl` at its `commitId`, and cached under its **registry name**: `~/.nextflow/assets/.repos/<scope>/<name>/`.

- `AssetManager` keeps two names. The project name is the local name: the cache folder, `nextflow list`, and `workflow.projectName`. The remote project is the Git one, used only by the Git provider (API calls and clone URL), and so by `workflow.repository`.
- When a cached multi-revision clone is used by its plain name, as in `pull`, `info`, `drop` or `run -offline`, the remote project is read from the clone's `origin` remote. So these commands work on a registry pipeline by its registry name, with no change in the commands.
- For a URL or a registry address, a clone of a different repository under the same name is an error ("already downloaded from a different repository"). This extends today's check, which only compares the host.
- For a plain name, the cached clone decides, as today: a plain `foo/bar` already uses a clone from another hub when one is cached under that name. When the clone's Git project differs from the name, Nextflow logs `Using <name> cached from <Git project>`.

**Considered cache options.**
- **Cache under the Git name.** Simplest: no `AssetManager` change. Rejected: the registry name is lost (`nextflow list` shows the Git name), and offline runs by registry name can't find the copy.
- **`AssetManager` with a local name (chosen).** About 15 lines in `AssetManager`. The existing commands keep working.
- **A separate registry folder** (`~/.nextflow/assets/.registry/<scope>/<name>/<version>/`) managed by `PipelineManager`, with `.pipeline-info`. It also fits M1 bundle releases, which have no Git. Not chosen for M0: `list`, `pull` and `drop` would not see it, and `CmdRun` would build the script file itself. M1 bundle releases need it or a similar cache; the M1 ADR decides.

## Future work

- **Registry addresses in other commands**: `pull`, `clone`, `config` and `info` could accept a registry URL or a registry host as `-hub`, so that a release can be fetched before an offline run.
- **A short prefix** for registry addresses.
- **Publish and install from a bundle**: `nextflow pipeline publish`, and installing a release from a bundle stored in the registry so that install does not depend on a Git host. This is decided in a separate ADR, `20261105-pipeline-registry-publish.md`, together with the server design in nextflow-registry #440.

## Open questions

- **`projectDir` in an included pipeline** still points to the meta-pipeline. Pipelines that read `${projectDir}/assets/...`, such as hlatyping and scrnaseq, break when included. The pipeline root could also become the included pipeline's `projectDir`, but that is a change to the language and runtime. It belongs to the composition ADR.
- **Typed pipelines**: inclusion requires `nextflow.enable.types = true`, which most nf-core pipelines don't enable yet.
- **Plugins** of an included pipeline must still be declared by the meta-pipeline (see the composition ADR).

### To be decided: running a registry pipeline

These points come from the prototype of explicit registry addresses and must be settled before merging.

- **Shape of a registry URL.** Recommended: keep the `/pipelines/` path and derive the page address from the API URL.
  - The prototype requires `/pipelines/`: `https://registry.nextflow.io/pipelines/nf-core/rnaseq` is the address of the registry page, so it can be copied from the browser. The registry keeps modules and pipelines in separate namespaces (`nf-core/rnaseq` may be both), and the path says which one is meant, as `/modules/` does for modules.
  - Without it (`https://registry.nextflow.io/nf-core/rnaseq`) the URL is shorter, but it is not a page address and it may clash with other registry paths.
  - The page address is the API URL without `/api`. This holds for every npr deployment, which serves pages and API from one host, and `pipeline view` builds its links the same way. A registry with another layout can still be used with its host as `-hub`, or a page URL option can be added to the `registry` scope later.
- **`-r` offline.** Offline, the registry can't be asked, so the prototype uses `-r` as a Git revision: a commit, or a tag equal to the version (`jorge-test` uses `v1.0.0` for `1.0.0`). Options:
  - A. A file next to the clone, `.repos/<scope>/<name>/.registry-releases`, with `version=commit` for each release run online. It sits beside `bare/` and `clones/`, outside any Git tree. The legacy layout would need multi-revision for registry runs.
  - B. The same mapping in the bare repository's Git config (`[registry "1.0.0"] commit = …`), which is never versioned and is removed with the clone. Recommended if this design is kept.
  - C. A local Git ref such as `refs/registry/1.0.0`. It needs extra code to resolve `-r`, and a local tag would look like a remote one in `nextflow info`.
  - D. One central file for all registry pipelines. `drop` would leave stale entries.
  - E. Keep `-r` as a Git revision offline, with a clear error when it does not resolve.
- **A registry host as hub in other commands.** Only `run` and `inspect` understand a registry host as `-hub` or `NXF_HUB`. `pull`, `clone`, `config` or `info` of a plain name that is not cached fail with "Unknown repository provider: registry.nextflow.io". Options:
  - A. A clear error in `AssetManager.checkHubProvider`: a registry host is only supported by `nextflow run`. Recommended for now.
  - B. Support registry addresses in `pull` now, so that a release can be fetched before an offline run.
  - C. Ignore a registry host in `NXF_HUB` in the other commands and use GitHub. Rejected: a silent switch.

## Links

- Implements: [Pipeline composition](20260608-pipeline-composition.md), section "Pipeline registry and CLI"
- Related: [Module system](20251114-module-system.md), [Workflow modules](20260608-workflow-modules.md)
- Spec: [docs/specs/261103-pipeline-registry-cli](../specs/261103-pipeline-registry-cli/spec.md)
- Registry: seqeralabs/nextflow-registry #418, #419, #421; client read methods #439
