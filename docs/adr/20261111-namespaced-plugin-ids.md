# Namespaced plugin IDs

- Authors: Jorge Ejarque
- Status: proposed
- Deciders: Ben Sherman, Phil Ewels, Paolo Di Tommaso
- Date: 2026-11-12
- Tags: plugins, registry

Technical Story: https://github.com/nextflow-io/nextflow/issues/8077

## Summary

Support plugin IDs of the form `<namespace>/<name>` (e.g. `myorg/nf-foo`), next to the existing unscoped IDs (e.g. `nf-amazon`). Reuse the namespace model that modules already use in the Nextflow registry.

This decision needs changes in three repositories: Nextflow, the plugin registry and the Gradle plugin (`nextflow-plugin-gradle`). Each repository has a copy of this ADR. This copy gives the full detail for Nextflow and a short summary for the other two.

## Problem Statement

Plugin IDs are global and flat. Two organisations cannot publish a plugin with the same name, and there is no ownership of a name beyond "first to claim it". Private registries make this worse: each one needs its own names, and an organisation cannot reserve its names in the public registry.

Modules already solved this with scoped names (`nf-core/fastqc`, see [module system ADR](20251114-module-system.md)). Plugins should follow the same model.

Today a `/` in a plugin ID fails in several places:

- Nextflow uses the ID in file paths (lock file, plugins store, local plugin root, offline repository).
- The registry rejects a `/` in the ID: `PluginId` validation fails on claim, publish, metadata and download.
- The Gradle plugin takes the ID from `project.name`, which cannot contain `/`.

## Goals or Decision Drivers

- An organisation can own a namespace and publish plugins under it, in the public registry and in private registries.
- The same plugin name can exist in different namespaces, and both can be used in the same run.
- Existing unscoped plugins keep working without changes, for plugin authors and for users.
- Old and new versions of Nextflow, the registry and the Gradle plugin work together for unscoped plugins. A namespaced plugin used with an old component fails with a clear error.
- One naming model for plugins and modules: same ownership, same API shape.

## Non-goals

- Authentication for private registries. Nextflow sends no token when it fetches plugin metadata or downloads a plugin. This is needed for private registries but is a separate decision.
- Moving existing plugins into namespaces (e.g. `nf-amazon` to `nextflow/nf-amazon`).
- Making the namespace required for new plugins.
- Lookup across several registries. It already works with `registry.url` and does not change.

## Considered Options

### 1. Keep flat IDs, use a naming convention

Plugin authors add a prefix to the name, e.g. `myorg-nf-foo`.

- **Pro:** No change in Nextflow, the registry or the Gradle plugin.
- **Con:** No ownership. Anyone can claim `myorg-*`.
- **Con:** Different from modules.

### 2. Namespace only in the registry

The registry stores plugins under a namespace, but the plugin ID stays flat (`nf-foo`). The registry maps the flat ID to the namespace.

- **Pro:** No change in Nextflow file paths.
- **Con:** A pipeline cannot say which `nf-foo` it wants. Two plugins with the same name cannot be used in the same run, because pf4j keys plugins by ID.
- **Con:** The ID in the MANIFEST and the ID in the registry differ, which is confusing.

### 3. Namespaced plugin IDs, same model as modules

The namespace is part of the plugin ID everywhere: config, MANIFEST, registry and local store.

- **Pro:** Same model as modules. Clear ownership and no clashes.
- **Pro:** pf4j needs no change. It only checks that the ID is not empty, and it parses `Plugin-Dependencies: org/foo@>=1.0` correctly.
- **Con:** Changes in Nextflow, the registry and the Gradle plugin.
- **Con:** Older Nextflow versions cannot use namespaced plugins.

## Solution

Option 3.

### Plugin ID

- Format: `[<namespace>/]<name>`, with at most one `/`.
- Unscoped IDs are not changed. `nf-foo` and `myorg/nf-foo` are different plugins.
- The namespace is optional.
- Naming rules for a namespaced ID are owned by the registry: the namespace follows the registry namespace rule (`^[a-z0-9]([a-z0-9-]{1,62}[a-z0-9])?$`), and the name follows the current plugin name rule. Unscoped IDs keep their current rules in each component, so no existing plugin is refused.
- The ID is used as is in config (`plugins { id 'myorg/nf-foo@1.0.0' }`), in `include from 'plugin/myorg/nf-foo'`, in `nextflow plugin myorg/nf-foo:<cmd>` and in `Plugin-Dependencies`.
- The `nf-core` namespace is protected in the registry, so nobody can publish `nf-core/...` plugins. The `nextflow` namespace is not protected. Seqera claims it as a normal namespace ( same as modules or pipelines ).

### Nextflow

Nextflow does not check the naming rules. It only checks the structure of the ID, because it turns the ID into file paths, and some sources (offline mode, test repositories, local plugins) never call the registry.

- **ID check:** `PluginRef.parse()` refuses an ID with more than one `/`, an empty part, or a part equal to `.` or `..`.
- **Plugins store:** `$NXF_HOME/plugins/<namespace>/<name>-<version>`. Unscoped plugins stay at `$NXF_HOME/plugins/<name>-<version>`. This is the same layout as modules and pipelines. The current code `pluginsStore.resolve("$id-$version")` already builds this path. Nextflow must create the parent directory before the atomic move.
- **Lock file:** use a flat name: `nextflow-plugin-<namespace>~<name>-<version>.lock`.
- **Local plugin root** (`LocalPluginManager.createLinkFromPath`): the symlink name must be flat and unique, e.g. `<namespace>~<name>-<version>`. Today it uses only the directory name, so `org1/nf-foo` and `org2/nf-foo` would overwrite each other.
- **Offline mode** (`LocalUpdateRepository`): scan one level deeper, into namespace directories. The ID still comes from the MANIFEST.
- **Download** (`PluginUpdater.download0`): no change. The unzipped directory is named after the zip file (`<name>-<version>.zip`), and the check compares only the last path element.
- **Embedded mode** (`EmbeddedPluginManager`): it loads every directory of the store root as a plugin. It must treat a directory without a MANIFEST as a namespace directory, the same as offline mode.
- **Registry refuses a namespaced ID:** an old registry answers 400 to a metadata request with a `/` in an ID. A new registry can also refuse one (bad name, old Nextflow). Nextflow shows a clear message with the registry response. Today metadata errors are only logged as warnings, and the run fails later with "Cannot find version". Nextflow adds the metadata error to that final error.

### Registry (summary)

The registry reuses its namespace model for modules and pipelines: ID validation with an optional `<namespace>/` prefix, a `namespace_id` column on plugins, publish permission through namespace membership, `%2F`-encoded IDs in API routes, download file name `<name>-<version>.zip`, a `plugins/` prefix in OCI storage, namespaced plugins left out of the legacy `plugins.json`, and namespace revoke, delete and "My plugins" support. It refuses a namespaced ID when the request's `nextflowVersion` is older than a configured minimum, with a clear message. See the ADR in the registry repository.

### Gradle plugin (summary)

A new optional `nextflowPlugin { namespace = 'myorg' }` setting. The ID becomes `myorg/<project.name>` in the MANIFEST, the release request and the local install path (`plugins/myorg/<name>-<version>`). The zip stays `<project.name>-<version>.zip`. See the ADR in the `nextflow-plugin-gradle` repository.

### Compatibility

| Nextflow | Registry | Unscoped plugin | Namespaced plugin |
|---|---|---|---|
| old | old | works | not possible |
| new | old | works | fails: clear message from Nextflow |
| old | new | works | fails: the registry message (`nextflowVersion` check) is shown as a warning, then "Cannot find version" |
| new | new | works | works |

| Gradle plugin | Registry | Unscoped plugin | Namespaced plugin |
|---|---|---|---|
| old | any | works | not possible (no `namespace` setting) |
| new | old | works | fails: the registry refuses the ID |
| new | new | works | works |

Nextflow versions that read the legacy `plugins.json` index never see namespaced plugins.

### Rollout

1. Registry and Gradle plugin, in any order.
2. Nextflow. It can be done in parallel and tested with a local `plugins.json` index, because it changes nothing for unscoped plugins.
3. Docs: ID format, store layout, minimum Nextflow version.

## Rationale

- Plugins and modules live in the same registry. One naming and ownership model is simpler for users and for the registry.
- Putting the namespace in the ID (instead of only in the registry) lets a pipeline choose the exact plugin, and lets two plugins with the same name run together.
- pf4j does not need changes. All the Nextflow work is in places where the ID becomes a file name.

## Consequences

- Pipelines that use namespaced plugins need a minimum Nextflow version. They should set `manifest.nextflowVersion`.
- With an old registry, only the namespaced plugin fails. Nextflow asks again one plugin per request, so unscoped plugins of the same run still resolve.
- Old Nextflow versions show the registry message only as a warning. This cannot be fixed in versions already released.
- Unscoped plugins, their store paths and older Nextflow versions are not affected.

