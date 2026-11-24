# Plugin namespace in the MANIFEST

- Authors: Jorge Ejarque
- Status: proposed
- Deciders: Ben Sherman, Phil Ewels, Paolo Di Tommaso
- Date: 2026-11-24
- Tags: plugins, registry

Technical Story: https://github.com/seqeralabs/nextflow-registry/pull/451#issuecomment-7130306191 (proposal), https://github.com/seqeralabs/nextflow-registry/pull/451 (link and freeze ADR)

## Summary

A plugin keeps one record and one list of releases when it moves into a namespace. The zip MANIFEST keeps the plain name in `Plugin-Id`. A new attribute, `Plugin-Namespace`, holds the namespace:

```
Plugin-Id: nf-foo
Plugin-Namespace: myorg
```

- Older Nextflow ignores `Plugin-Namespace` and loads the zip as `nf-foo`, as today.
- Newer Nextflow loads the same zip as `myorg/nf-foo`.
- The registry resolves `nf-foo` and `myorg/nf-foo` to the same record and the same releases, including the releases published before the move.
- The registry tells Nextflow the namespaced ID of a moved plugin. Newer Nextflow changes every plain ID of that plugin to the namespaced ID. A run then loads the plugin only once.
- After an enforcement date, an existing plugin must move before it can publish a new version, and new plugins must be namespaced.

This builds on the [namespaced plugin IDs ADR](20261111-namespaced-plugin-ids.md) and is an alternative to the link and freeze ADR (#451), which is option 1 below. It changes some rules of the namespaced plugin IDs ADR (see [Changes to the namespaced plugin IDs ADR](#changes-to-the-namespaced-plugin-ids-adr)).

## Problem Statement

With namespaced plugin IDs, an author who wants a namespace must publish the plugin again under a new ID. The new ID has no releases, no download stats and no link to the old plugin. Users are not told about the new ID, and every pipeline must change the ID by hand.

The registry alone cannot make `nf-foo` serve the releases of `myorg/nf-foo`:

- pf4j, the plugin library in Nextflow, registers a plugin under the `Plugin-Id` in its zip MANIFEST. It cannot map one plugin ID to another.
- Nextflow downloads a plugin by the ID it asked for, then calls `startPlugin(id)` with that ID (`PluginUpdater.load0`).
- In an experiment, a zip with `Plugin-Id: myorg-nf-foo` was served as `nf-foo`, and the reverse. Both failed with `PluginNotFoundException: Plugin '<requested id>' not found`.

So a zip can only be loaded under the ID in its MANIFEST, unless Nextflow changes that ID when it reads the MANIFEST.

## Goals or Decision Drivers

- **One plugin, one record.** A moved plugin has one list of releases, one set of stats and one page. The author publishes each release once.
- **Fast adoption of namespaces.** Moving is cheap for authors and required after a set date. After that date, no new plugins without a namespace are added.
- **Backward compatibility.** Pipelines that use `nf-foo`, including pinned versions such as `nf-foo@1.2.0`, and Nextflow versions without namespace support keep working and keep getting new releases, with no change on their side.
- **Namespaces stay part of the ID** on newer Nextflow, so `org1/nf-foo` and `org2/nf-foo` can exist and run together.

## Non-goals

- Moving a plugin to a different name (`nf-foo` to `myorg/nf-bar`).
- Moving a plugin from one namespace to another.
- Authentication for private registries.
- An audit trail of who moved a plugin.

## Considered Options

### 1. Link, then freeze

`nf-foo` and `myorg/nf-foo` stay two plugins, linked in the registry. Both accept releases until the author stops releases for `nf-foo` (the freeze). The zip `Plugin-Id` is namespaced for `myorg/nf-foo`.

- **Pro:** no change to how Nextflow loads plugins, beyond a warning and a duplicate check.
- **Pro:** the author decides how long to support older Nextflow versions.
- **Con:** two records and two release lists. While linked, the author builds and publishes each release two times.
- **Con:** after the freeze, pipelines that use `nf-foo` and older Nextflow versions get no new releases.
- **Con:** a run could ask for both IDs. Nextflow must detect this case.
- **Con:** moving is optional and costly, so most plugins would stay unscoped and both models would stay forever.

Not chosen: it keeps two models at the same time, and the author pays for it on every release.

### 2. Link, plus an alias in Nextflow

As option 1, and newer Nextflow treats the linked pair as one plugin. It merges their release lists and maps one ID to the other when it loads a zip.

- **Pro:** users keep `nf-foo` in their config and get the releases of both IDs on newer Nextflow.
- **Con:** all the costs of option 1 (two records, two publishes), plus the Nextflow alias work.
- **Con:** older Nextflow still gets only what the author publishes under `nf-foo`.

Not chosen: the chosen option needs the same Nextflow work, without the second record.

### 3. Duplicate storage

As option 1, and the registry stores a copy of each `myorg/nf-foo` release with `Plugin-Id: nf-foo`, served as an `nf-foo` release.

- **Pro:** no Nextflow change. The author publishes once.
- **Con:** the registry changes the author's zip. The copy has a different checksum from the file the author built.
- **Con:** each release needs two OCI artifacts. The registry must download, unpack and repack the zip. Today it only streams the upload.
- **Con:** a run could ask for both IDs. Nextflow must detect this case.

Not chosen: the registry would rewrite and store every release twice.

### 4. Rename

Change the ID of all rows from `nf-foo` to `myorg/nf-foo`, with the namespace in `Plugin-Id`.

- **Con:** old zips say `Plugin-Id: nf-foo` and cannot be loaded as `myorg/nf-foo`.
- **Con:** every pipeline with `nf-foo`, and every plugin with `Plugin-Dependencies: nf-foo@...`, breaks at once.
- **Con:** Nextflow versions without namespace support lose the plugin.

Not chosen: it breaks pipelines and older Nextflow versions.

### 5. Unique names, namespace for ownership only

A plugin name is unique in the whole registry. The namespace only says who owns the plugin. `myorg/nf-foo` in a config is another way to write `nf-foo`.

- **Pro:** old pipelines and old Nextflow versions keep working. A plugin cannot be loaded twice, because it has one ID.
- **Con:** two organisations cannot use the same name, an organisation cannot reserve its names, and private registries can reuse names from the public one. These are the problems that namespaces solve.

Not chosen: it makes moving easy, but it removes the reason for namespaces.

### 6. `Plugin-Namespace` in the MANIFEST, one record with two IDs (chosen)

`Plugin-Id` stays the plain name and `Plugin-Namespace` holds the namespace. A move sets the namespace on the existing record. Newer Nextflow builds the ID from both attributes and maps plain IDs of moved plugins to the namespaced one.

- **Pro:** one record, one release list, one stored zip per release, no link or freeze states.
- **Pro:** older Nextflow loads the same zip as `nf-foo` and keeps getting new releases.
- **Pro:** the namespace stays part of the ID on newer Nextflow.
- **Pro:** pf4j needs no change. Nextflow subclasses the pf4j descriptor finder, whose methods are `protected`.
- **Con:** the ID in the MANIFEST and the ID on newer Nextflow differ.
- **Con:** the largest Nextflow change: the ID mapping, the rewrite of dependencies, and the saved mapping for offline runs. This is the alias work of option 2.
- **Con:** Nextflow must get the namespaced ID of every moved plugin from the registry before it loads anything. Without it, a plugin can be loaded twice.
- **Con:** some rules of the namespaced plugin IDs work (#443) change.

## Solution

Option 6.

### Plugin MANIFEST

| Plugin | `Plugin-Id` | `Plugin-Namespace` | `Plugin-Dependencies` |
|---|---|---|---|
| Unscoped, not moved | `nf-foo` | none | plain IDs |
| Moved, release built before the move | `nf-foo` | none | plain IDs |
| Moved, release built after the move | `nf-foo` | `myorg` (may be missing, see Publish) | plain IDs while the author supports older Nextflow |
| New namespaced plugin (never had a plain name) | `nf-foo` | `myorg` | namespaced IDs, or the plain ID of a moved plugin |

- `Plugin-Id` is always the plain name and never contains `/`. pf4j reads only this attribute, so older Nextflow is not affected by the new one.
- A moved plugin keeps plain dependencies while it supports older Nextflow: older Nextflow cannot resolve a namespaced dependency.
- A new namespaced plugin depends on namespaced IDs. A plain ID works only for a plugin that has a plain name (unscoped or moved).

### Registry

**Move.**

- **Who:** the owner of `nf-foo` who is a member of the target namespace, or an admin.
- **What:** sets `namespace_id` on the existing `npr_plugin` row. The row ID stays `nf-foo`, so releases, claims and usage rows do not change.
- **Rules:**
  - Same name only. Only an unscoped plugin moves.
  - The publish rules of the namespace apply, so a move into a protected namespace (`nf-core`) is refused.
  - A move is refused if a row `myorg/nf-foo` already exists. After a move, creating `myorg/nf-foo` is refused.
  - After a move, a new claim on `nf-foo` is refused. The owners are the namespace members.
  - The author cannot undo a move. An admin can, only before the first release with `Plugin-Namespace`, because newer Nextflow loads those releases as `myorg/nf-foo`.
- **Where:** a UI action on the plugin page. No public API.

**Both IDs resolve to one record.** Every read and write path first resolves the requested ID to the row, then uses the row ID and `namespace_id`. This covers metadata, download, `/dependencies`, publish (draft and upload), validation, permissions, OCI tags, import, pages, mail and usage.

- `nf-foo` finds the row by its ID.
- `myorg/nf-foo` finds the row with ID `myorg/nf-foo`, or else the row with ID `nf-foo` and `namespace_id` of `myorg`.
- Permission comes from the row's `namespace_id`, not from the requested ID. An older Gradle plugin that sends `nf-foo` for a moved plugin gets namespace-membership permission.
- Both IDs return the same releases and download the same zip.

**Hint for Nextflow.** `/dependencies` returns one entry per requested ID, with that ID in `id`, so older Nextflow finds it under the ID it asked for. Each entry of a moved plugin also has `canonicalId: "myorg/nf-foo"`.

- When one request has both IDs (for example `myorg/nf-foo` from the config and `nf-foo` from a dependency), the response has both entries.
- The plain name is the part of `canonicalId` after the `/`, because a move keeps the name. So one field maps the plain name to the namespaced ID, and back.
- New namespaced plugins and unscoped plugins have no `canonicalId`.
- `/dependencies` already returns the dependencies of the requested plugins (up to 5 levels), so Nextflow gets every hint in one call.
- The field is added in a new TypeSpec API version with `@added`. Nextflow decodes the response with Gson, which ignores unknown fields.

**`nextflowVersion` check.** The minimum Nextflow version for namespaced IDs (namespaced plugin IDs ADR) is the first Nextflow release that reads `Plugin-Namespace`. The check already applies only to the requested ID, so a moved plugin requested by its plain name passes. Dependencies found during the walk that have no plain name are left out of the response for older clients.

**Publish.**

- `Plugin-Id` must equal the plugin name and must not contain `/`.
- `Plugin-Namespace`, when present, must equal the namespace of the row. An unscoped plugin that did not move must not have it.
- A moved plugin may publish without `Plugin-Namespace`. This allows patch releases built from older branches (for example for core plugins). Newer Nextflow loads them as `myorg/nf-foo` through the hint.
- The registry warns when a release of a moved plugin depends on a plugin with no plain name, because older Nextflow cannot load that release.

**Enforcement.** After the enforcement date:

- A publish to an existing unscoped plugin that has not moved is refused, with a message that says how to move it.
- A claim of a new unscoped plugin is refused.

**Namespace revoked or deleted.**

- **Revoked:** `myorg/nf-foo` returns 410, as in the namespaced plugin IDs ADR. `nf-foo` keeps serving its releases, so pipelines and older Nextflow that use the plain name are not affected. Publishing is refused while the namespace is revoked.
- **Deleted:** refused while the namespace holds moved plugins. An admin first undoes the moves or deletes the plugins.

**Pages and lists.** One page and one list entry per plugin, reached from both IDs. The page shows the namespaced ID once the plugin has moved. Stats need no merge, because there is one record.

**Legacy `plugins.json`.** It keeps listing moved plugins under the plain name, with the same release rules as today. Older Nextflow loads them as before.

### Gradle plugin

- `nextflowPlugin { namespace = 'myorg' }` writes `Plugin-Namespace: myorg` and sends `myorg/<project.name>` in the release request.
- `Plugin-Id` stays `<project.name>`. The zip name stays `<project.name>-<version>.zip`.
- The local install path is `plugins/myorg~<project.name>-<version>`, as in the namespaced plugin IDs work.
- A refused publish shows the registry message, as today.

### Nextflow

**Plugin ID.** A `NamespaceDescriptorFinder` (subclass of pf4j `ManifestPluginDescriptorFinder`) sets the ID when it reads the MANIFEST:

1. `Plugin-Id` as is, when it contains `/` (a zip built with a namespaced `Plugin-Id`);
2. `Plugin-Namespace` + `/` + `Plugin-Id`, when `Plugin-Namespace` is present;
3. the namespace from the store folder name, for a release built without `Plugin-Namespace` (`myorg~nf-foo-1.0.0` gives `myorg/nf-foo`);
4. otherwise, the plain `Plugin-Id`.

If the request was namespaced and the zip declares a different `Plugin-Namespace`, Nextflow refuses the plugin. A different namespace is accepted only for a plain request, because that is how a moved plugin answers.

`LocalPluginManager`, `EmbeddedPluginManager` and the offline repository (`LocalUpdateRepository`) use this finder.

**ID mapping.** Nextflow keeps a map from plain name to namespaced ID. Sources:

- the `canonicalId` hints of `/dependencies`;
- a zip served for a plain request whose `Plugin-Namespace` gives a namespaced ID;
- in embedded mode, the loaded plugins that declare `Plugin-Namespace`;
- the store links (see below), online and offline.

Nextflow never maps `nf-foo` to `myorg/nf-foo` from names alone, because `nf-foo` can be a different plugin. Store links are trusted because only Nextflow writes them, and only from a hint or a served zip.

Order in a run:

1. Read the store links.
2. Fetch the metadata of the config, `-plugins` and default plugins. The hints arrive here.
3. Map the IDs, then merge the config and default plugins (the config ID wins), then check `NXF_PLUGINS_ALLOWED` (both sides mapped).
4. Load the plugins. In the dependency loop, map each `Plugin-Dependencies` entry in memory before pf4j resolves it. The MANIFEST file is not changed.

The metadata of a plugin is kept under both IDs, so a lookup by either ID finds it, also with a `plugins.json` index that has only the plain ID.

The map also applies to `include ... from 'plugin/nf-foo'`, on-demand starts (`startIfMissing`), `isStarted`, `getPlugin`, `getExtensions`, `nextflow plugin install`, `nextflow plugin nf-foo:<cmd>`, updates, and the core plugin versions pinned in `DefaultPlugins` (which answer for both IDs).

Nextflow skips the registry call for pinned plugins already in the store. When a run has a plain ID and a namespaced ID with the same name and no link between them, it makes the call anyway, to learn whether they are one plugin.

**Notice.** When the config or `-plugins` uses the plain name of a moved plugin, Nextflow logs a warning once per plugin: "Plugin `nf-foo` has moved to `myorg/nf-foo`. Update your config." The run goes on.

**Plugin store.** A moved plugin is stored once. A symlink gives it the other name:

- Online, with the hint, Nextflow stores the plugin in `myorg~nf-foo-<ver>` and links `nf-foo-<ver>` to it.
- If an older Nextflow already stored `nf-foo-<ver>`, Nextflow links `myorg~nf-foo-<ver>` to that folder and does not download again.
- A link always points to a real, complete folder, never to another link.
- Newer Nextflow takes one download lock for both names (the lock of the namespaced ID). Older Nextflow takes the lock of the plain ID. If both run at once, older Nextflow replaces the link with a real folder (`REPLACE_EXISTING`). This is harmless: the store then keeps two copies.
- If the link cannot be created (read-only store, no symlink support), Nextflow logs it and goes on. The plugin still loads.

**Offline.**

- Nextflow reads the map from the store links, one link per version. It skips broken links.
- A store filled only by an older Nextflow has no links, so it gives no map.
- The offline docs say: fill the store online with the Nextflow version used offline, and copy it with links kept (`cp -a`, `rsync -a`).
- `nextflow plugin install` asks the registry for every requested plugin, also for plugins already in the store, so it can add the missing links.

**Core plugins.** The hardcoded IDs change to `nextflow-io/...` in the same Nextflow release as the namespace support. Core zips start writing `Plugin-Namespace` only from that release.

**Dev mode.** `DevPluginManager` does not use the new finder. A plugin author testing a namespaced plugin locally refers to it by its plain name.

### Compatibility

| Nextflow | Asks for `nf-foo` (moved) | Asks for `myorg/nf-foo` |
|---|---|---|
| reads `plugins.json` only | releases as today, loaded as `nf-foo` | not possible |
| `/dependencies`, no namespace support | all releases, loaded as `nf-foo` | not possible |
| namespace support | loaded as `myorg/nf-foo`, with the notice | loaded as `myorg/nf-foo` |

What older Nextflow does not get:

- new namespaced plugins, which have no plain name;
- releases of a moved plugin that depend on a plugin with no plain name (the registry warns at publish).

Older Nextflow sharing a plugin store with newer Nextflow, offline, with a version range:

- It reads every folder in the store and keys it by `Plugin-Id`. A folder `org1~nf-foo-2.0.0` of a new namespaced plugin shows as a release of `nf-foo`. Older Nextflow then fails with "is not a regular file".
- A broken link (for example after the user deletes the folder it points to) makes it fail with "Cannot find the manifest path".
- Pinned versions are not affected.

### Worked example: core plugins (`nf-amazon` to `nextflow-io/nf-amazon`)

Core plugins are started by Nextflow itself. Each Nextflow release has their IDs and versions in its code.

Where Nextflow names core plugins (checked at Nextflow commit `f2e306cbe`):

| Where | IDs | Started when |
|---|---|---|
| `FileHelper.PLUGINS_MAP` | `nf-amazon`, `nf-google`, `nf-azure`, `nf-tower` | a path uses `s3://`, `gs://`, `az://` or `seqera://` |
| `PluginsFacade` (default plugins) | `nf-tower`, `nf-wave`, `nf-seqera`, `nf-cloudcache`, `nf-amazon`, `nf-google`, `nf-azure`, `nf-k8s`, `nf-agent`, `nf-agent-pi` | config such as `tower.enabled`, `wave.enabled` or `fusion.enabled`, a cloud executor, or a cloud work dir |
| `RepositoryFactory` | `nf-codecommit` | a CodeCommit repository URL |
| `CmdAuth`, `CmdLaunch` | `nf-tower` | `nextflow auth`, `nextflow launch` |
| `META-INF/plugins-info.txt` (`DefaultPlugins`) | all core plugins | pins the version of each core plugin for that Nextflow release |

Steps:

1. **Registry:** move `nf-amazon` to `nextflow-io`.
2. **Nextflow:** in the release with namespace support, change the IDs above to `nextflow-io/nf-amazon`, and build the core plugins with `Plugin-Namespace: nextflow-io`.

Effect:

- Older Nextflow keeps starting `nf-amazon` at its pinned version, and gets new releases through the plain name.
- Newer Nextflow starts `nextflow-io/nf-amazon`. A config that names `nf-amazon` maps to it, with the notice.
- Plugins that depend on `nf-amazon` (for example `xpack-amzn`, `requirePlugins = ['nf-amazon@>=3.13.0']`) keep working: Nextflow maps the dependency and keeps the pinned core version. They need no new release.
- Patch releases for older Nextflow lines are published to the same record without `Plugin-Namespace`.

### Worked example: xpack (`xpack-amzn` to `seqera/xpack-amzn`)

`BaseXpackPlugin` uses the plugin ID as the license product code (`product = wrapper.descriptor.pluginId`). On newer Nextflow the ID is `seqera/xpack-amzn`, so the license server would get the wrong code, and users would see "Seqera/xpack-amzn license" in errors.

Required change in xpack, before the first build with `Plugin-Namespace`: use the name without the namespace as the product code. `core/xpack` is bundled into each plugin, so every `seqera/xpack-*` build must include it.

```groovy
product = wrapper.descriptor.pluginId.tokenize('/').last()
```

### Rollout

1. Registry and Gradle plugin releases, with enforcement off.
2. Nextflow release with namespace support. It contains the namespaced plugin IDs work and this ADR together, and the core plugin ID change.
3. Move the core plugins to `nextflow-io` and the Seqera plugins to `seqera`. Apply the xpack product-code change before the first build with `Plugin-Namespace`.
4. Announce the enforcement date to plugin authors, with a how-to for claiming a namespace and moving a plugin.
5. Turn enforcement on.

### Changes to the namespaced plugin IDs ADR

- `Plugin-Id` is always the plain name. The namespace goes in `Plugin-Namespace`.
- The Gradle plugin writes `Plugin-Namespace` instead of a namespaced `Plugin-Id`. Its install path is `plugins/myorg~<name>-<version>`.
- The registry checks `Plugin-Id` against the name and `Plugin-Namespace` against the namespace.
- Nextflow takes the ID from `Plugin-Namespace`, or from the store folder name. So `LocalUpdateRepository`, `EmbeddedPluginManager` and `LocalPluginManager` change: they use the new finder.
- The minimum Nextflow version for namespaced IDs is the first release that reads `Plugin-Namespace`.
- The store and lock file layout (`<namespace>~<name>-<version>`) and the ID check in `PluginRef.parse` do not change.

## Prototype results

A prototype in Nextflow (on top of nextflow-io/nextflow#8082) tested the Nextflow part. The registry hint was simulated with a map.

- pf4j loads a zip as `<Plugin-Namespace>/<Plugin-Id>` through the `protected` descriptor hooks. Released Nextflow 26.10.4 loads a zip with `Plugin-Namespace` with no warning.
- With the hint, a moved plugin loads once, whatever ID the config, `include`, a dependency or an on-demand start uses.
- Without the hint:
  - a plain ID requested after the namespaced one fails with "There is an already loaded plugin", or starts the plugin twice for a release built without `Plugin-Namespace`;
  - a plain dependency on a moved plugin that is not loaded yet fails with `DependenciesNotFoundException`.
- Offline, the store links give the map. Released Nextflow 26.10.4 works offline on a store with links, for pinned versions and version ranges.

Other store layouts considered for the offline map: a marker file inside the plugin folder (safe, but older and newer Nextflow keep two copies) and a separate store for newer Nextflow. Links keep one copy and need no change in older Nextflow.

## Rationale

- One record and one zip per release remove the cost that options 1 to 3 put on authors and on the registry.
- `Plugin-Id` stays plain, so older Nextflow loads every release without change. This is what lets moving become required.
- The namespace still becomes part of the ID on newer Nextflow, which keeps the reasons for namespaces (option 5 does not).
- The price is in Nextflow: the ID mapping, which needs the registry hint, and the saved map for offline runs.

## Consequences

- Newer Nextflow needs the `canonicalId` hint before it loads a plugin. A registry or index without it can load a moved plugin twice.
- Offline runs need a store filled online by a Nextflow version with namespace support.
- Older Nextflow sharing a store with newer Nextflow can fail offline with a version range (see Compatibility).
- The ID shown by Nextflow (`myorg/nf-foo`) differs from `Plugin-Id` in the MANIFEST (`nf-foo`).
- Authors of moved plugins keep plain dependencies while they support older Nextflow.
- After the enforcement date, no new unscoped plugins are added.

## Open decisions

- The enforcement date, and how long after the Nextflow release it falls.
- How an admin allows a single release that enforcement would refuse.
- Whether an admin can move a plugin whose author does not respond.
- Docs: ID format, the `Plugin-Namespace` attribute, the migration guide, and offline store preparation.
