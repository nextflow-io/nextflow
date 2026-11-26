# Namespaced plugin IDs

- Authors: Jorge Ejarque
- Status: proposed
- Deciders: Ben Sherman, Phil Ewels, Paolo Di Tommaso
- Date: 2026-11-12
- Tags: plugins, registry

Technical Story: https://github.com/nextflow-io/nextflow/issues/8077, https://github.com/seqeralabs/nextflow-registry/pull/451 (link and freeze, option B2), https://github.com/seqeralabs/nextflow-registry/pull/451#issuecomment-7130306191 (`Plugin-Namespace` proposal)

## Summary

Support plugin IDs of the form `<namespace>/<name>` (e.g. `myorg/nf-foo`), next to the existing unscoped IDs (e.g. `nf-amazon`). Reuse the namespace model that modules already use in the Nextflow registry.

- The MANIFEST (`META-INF/MANIFEST.MF` inside the plugin zip) keeps the plain name in `Plugin-Id`. A new attribute, `Plugin-Namespace`, holds the namespace.
- An existing unscoped plugin can move into a namespace. It keeps one registry record and one list of releases, and the registry accepts both IDs (`nf-foo` and `myorg/nf-foo`).
- Older Nextflow ignores `Plugin-Namespace` and loads a moved plugin as `nf-foo`, as today. Newer Nextflow loads the same zip as `myorg/nf-foo`.
- The registry tells Nextflow the namespaced ID of a moved plugin. Newer Nextflow changes every plain ID of that plugin to the namespaced ID, so a run loads the plugin only once.
- Once the registry turns on enforcement (a setting that makes namespaces required), an existing plugin must move before it can publish a new version, and new plugins must be namespaced.

This decision needs changes in three repositories: Nextflow, the plugin registry and the Gradle plugin (`nextflow-plugin-gradle`). Each repository has a copy of this ADR. This copy gives the full detail for Nextflow and a short summary for the other two.

## Problem Statement

Plugin IDs are global and flat. Two organisations cannot publish a plugin with the same name. Nobody owns a name: the first author to claim it gets it. Private registries make this worse. Each one needs its own names, and an organisation cannot reserve its names in the public registry.

Modules already solved this with scoped names (`nf-core/fastqc`, see [module system ADR](20251114-module-system.md)). Plugins should follow the same model.

Today a `/` in a plugin ID fails in several places:

- Nextflow uses the ID in file paths (lock file, plugins store, local plugin root, offline repository).
- The registry rejects a `/` in the ID: `PluginId` validation fails on claim, publish, metadata and download.
- The Gradle plugin takes the ID from `project.name`, which cannot contain `/`.

Second problem: existing plugins. With namespaced IDs alone, an author who wants a namespace must publish the plugin again under a new ID. The new ID has no releases, no download stats and no link to the old plugin, and every pipeline must change the ID by hand. The registry alone cannot make `nf-foo` serve the releases of `myorg/nf-foo`:

- pf4j, the plugin library in Nextflow, registers a plugin under the `Plugin-Id` in its zip MANIFEST. It cannot map one plugin ID to another.
- Nextflow downloads a plugin by the ID it asked for, then calls `startPlugin(id)` with that ID (`PluginUpdater.load0`).
- In an experiment, a zip with `Plugin-Id: myorg-nf-foo` was served as `nf-foo`, and the reverse. Both failed with `PluginNotFoundException: Plugin '<requested id>' not found`.

So a zip can only be loaded under the ID in its MANIFEST, unless Nextflow changes that ID when it reads the MANIFEST.

## Goals or Decision Drivers

- An organisation can own a namespace and publish plugins under it, in the public registry and in private registries.
- The same plugin name can exist in different namespaces, and both can be used in the same run.
- An existing plugin can move into a namespace and keep one record: one list of releases, one set of stats and one page. The author publishes each release once.
- Pipelines that use a plain ID keep working, including pinned versions such as `nf-foo@1.2.0`. Nextflow versions without namespace support also keep working. Both keep getting new releases of moved plugins, with no change on their side.
- A new namespaced plugin used with an old component fails with a clear error.
- Fast adoption: moving is cheap for authors, and the registry can make namespaces required.
- One naming model for plugins and modules: same ownership, same API shape.

## Non-goals

- Authentication for private registries. Nextflow sends no token when it fetches plugin metadata or downloads a plugin. This is needed for private registries but is a separate decision.
- Moving a plugin to a different name (`nf-foo` to `myorg/nf-bar`), or from one namespace to another.
- An audit trail of who moved a plugin.
- Changing how Nextflow combines several registries. It already works with `registry.url`. This ADR only says how a moved plugin fits in (see "Several registries").

## Considered Options

There are two questions:

- **A.** Should a plugin ID include a namespace? (A1 to A3)
- **B.** If yes, how does the zip carry the namespace, and how does an existing plugin move into a namespace? (B1 to B5)

### A. Plugin ID model

#### A1. Keep flat IDs, use a naming convention

Plugin authors add a prefix to the name, e.g. `myorg-nf-foo`.

- **Pro:** No change in Nextflow, the registry or the Gradle plugin.
- **Con:** No ownership. Anyone can claim `myorg-*`.
- **Con:** Different from modules.

#### A2. Namespace only in the registry

The registry stores plugins under a namespace, but the plugin ID stays flat (`nf-foo`). The registry maps the flat ID to the namespace.

- **Pro:** No change in Nextflow file paths.
- **Con:** A pipeline cannot say which `nf-foo` it wants. Two plugins with the same name cannot be used in the same run, because pf4j keys plugins by ID.

#### A3. Namespaced plugin IDs, same model as modules (chosen)

The namespace is part of the plugin ID in config, `Plugin-Dependencies`, the registry and the local store.

- **Pro:** Same model as modules. Clear ownership and no clashes.
- **Pro:** pf4j needs no change. It parses `Plugin-Dependencies: org/foo@>=1.0` correctly.
- **Con:** Changes in Nextflow, the registry and the Gradle plugin.
- **Con:** Older Nextflow versions cannot use new namespaced plugins.

### B. Namespace in the zip, and moving existing plugins

#### B1. Namespace inside `Plugin-Id`

The MANIFEST has `Plugin-Id: myorg/nf-foo`.

- **Pro:** The ID is the same in the MANIFEST, in Nextflow and in the registry.
- **Con:** A zip loads under one ID only. A moved plugin cannot serve `nf-foo` to pipelines and older Nextflow, so moving needs a second record (B2) or a second copy of each zip (B3).

#### B2. Link, then freeze

`nf-foo` and `myorg/nf-foo` stay two plugins, linked in the registry. Both accept releases until the author stops releases for `nf-foo` (the freeze).

- **Pro:** No change to how Nextflow loads plugins, beyond a warning and a duplicate check.
- **Pro:** The author decides how long to support older Nextflow versions.
- **Con:** Two records and two release lists. While linked, the author builds and publishes each release two times.
- **Con:** After the freeze, pipelines that use `nf-foo` and older Nextflow versions get no new releases.
- **Con:** A run could ask for both IDs. Nextflow must detect this case.
- **Con:** Moving is optional and costly, so most plugins would stay unscoped and both models would stay forever.

#### B3. Duplicate storage

As B2, and the registry stores a copy of each `myorg/nf-foo` release with `Plugin-Id: nf-foo`, served as an `nf-foo` release.

- **Pro:** No Nextflow change. The author publishes once.
- **Con:** The registry changes the author's zip. The copy has a different checksum from the file the author built.
- **Con:** Each release needs two OCI artifacts (the container-registry objects where the registry stores zips). The registry must download, unpack and repack the zip. Today it only streams the upload.
- **Con:** A run could ask for both IDs. Nextflow must detect this case.

#### B4. Unique names, namespace for ownership only

A plugin name is unique in the whole registry. The namespace only says who owns the plugin. `myorg/nf-foo` in a config is another way to write `nf-foo`.

- **Pro:** Old pipelines and old Nextflow versions keep working. A plugin cannot be loaded twice, because it has one ID.
- **Con:** Two organisations cannot use the same name, an organisation cannot reserve its names, and private registries can reuse names from the public one. These are the problems that namespaces solve (A3).

#### B5. `Plugin-Namespace` in the MANIFEST, one record with two IDs (chosen)

`Plugin-Id` stays the plain name and `Plugin-Namespace` holds the namespace. A move sets the namespace on the existing record. Newer Nextflow builds the ID from both attributes and maps plain IDs of moved plugins to the namespaced one.

- **Pro:** One record, one release list, one stored zip per release, no link or freeze states.
- **Pro:** Older Nextflow loads the same zip as `nf-foo` and keeps getting new releases.
- **Pro:** The namespace stays part of the ID on newer Nextflow.
- **Pro:** pf4j needs no change. Nextflow subclasses the pf4j descriptor finder, whose methods are `protected`.
- **Con:** The ID in the MANIFEST and the ID on newer Nextflow differ.
- **Con:** The largest Nextflow change: the ID mapping, the rewrite of dependencies, and the saved mapping for offline runs.
- **Con:** Nextflow must get the namespaced ID of every moved plugin from the registry before it loads anything. Without it, a plugin can be loaded twice.

## Solution

Options A3 and B5.

### Plugin ID

- Format: `[<namespace>/]<name>`, with at most one `/`.
- Unscoped IDs are not changed. `nf-foo` and `myorg/nf-foo` are different plugins, unless `nf-foo` moved into `myorg`.
- The namespace is optional until the registry turns enforcement on.
- The registry owns the naming rules for a namespaced ID. The namespace follows the registry namespace rule (`^[a-z0-9]([a-z0-9-]{1,62}[a-z0-9])?$`). The name follows the current plugin name rule. Unscoped IDs keep their current rules, so no existing plugin is refused.
- The ID is used as is in config (`plugins { id 'myorg/nf-foo@1.0.0' }`), in `include from 'plugin/myorg/nf-foo'`, in `nextflow plugin myorg/nf-foo:<cmd>` and in `Plugin-Dependencies`.
- `nf-core` is a protected namespace: nobody can publish `nf-core/...` plugins.
- `nextflow-io` and `seqera` are not protected. They are already claimed in production. Core and open-source plugins go in `nextflow-io`. Seqera plugins, such as `xpack`, go in `seqera`.

### Plugin MANIFEST

```
Plugin-Id: nf-foo
Plugin-Namespace: myorg
```

| Plugin | `Plugin-Id` | `Plugin-Namespace` | `Plugin-Dependencies` |
|---|---|---|---|
| Unscoped, not moved | `nf-foo` | none | plain IDs |
| Moved, release built before the move | `nf-foo` | none | plain IDs |
| Moved, release built after the move | `nf-foo` | `myorg`, or none for patch releases built from older branches | plain IDs while the author supports older Nextflow |
| New namespaced plugin (never had a plain name) | `nf-foo` | `myorg` | namespaced IDs, or the plain ID of a moved plugin |

- `Plugin-Id` is always the plain name and never contains `/`. pf4j reads only this attribute, so older Nextflow is not affected by the new one.
- A moved plugin keeps plain dependencies while it supports older Nextflow: older Nextflow cannot resolve a namespaced dependency.
- A new namespaced plugin depends on namespaced IDs. A plain ID works only for a plugin that has a plain name (unscoped or moved).

### Registry hint: `canonicalId`

`canonicalId` is a new optional field in the `GET /api/v1/plugins/dependencies` response. It gives the namespaced ID of a moved plugin. The registry fills it and Nextflow needs it.

- The response has one entry per requested ID. The `id` field is the requested ID, so older Nextflow finds the entry under the ID it asked for.
- Each entry of a moved plugin has `canonicalId: "myorg/nf-foo"`, whichever ID was requested.
- When one request has both IDs (for example `myorg/nf-foo` from the config and `nf-foo` from a dependency), the response has both entries.
- The plain name is the part of `canonicalId` after the `/`, because a move keeps the name. So one field maps the plain name to the namespaced ID, and back.
- New namespaced plugins and unscoped plugins have no `canonicalId`.
- `/dependencies` already returns the dependencies of the requested plugins (up to 5 levels), so Nextflow gets every hint in one call.
- `canonicalId` is a plain optional field in the existing API version. Nextflow decodes the response with Gson, which ignores unknown fields.
- Nextflow accepts a `canonicalId` only if it is a valid namespaced ID. For a plain request it must have the same name. For a namespaced request it must be the same ID.

### Nextflow

**ID check.** Nextflow does not check the naming rules. The registry does. Nextflow checks only the structure of a namespaced ID. It needs this check because it turns the ID into file paths, and some sources never call the registry (offline mode, test repositories, local plugins). For a namespaced ID, `PluginRef.parse()` refuses more than one `/`, an empty part, a part equal to `.` or `..`, and a `~`. Unscoped IDs are not checked, as today. This means an unscoped ID with a `~`, such as `a~b`, gets the same store folder as `a/b`. The public registry refuses `~` in plugin names, so this can only happen with local or test repositories.

**File names.** Every file name built from the ID uses a flat name, with `~` in place of `/` (`PluginRef.dirName`): `myorg/nf-foo` becomes `myorg~nf-foo`.

- Plugins store: `$NXF_HOME/plugins/<namespace>~<name>-<version>`. Unscoped plugins stay at `$NXF_HOME/plugins/<name>-<version>`.
- Lock file: `nextflow-plugin-<namespace>~<name>-<version>.lock`.
- The local plugin root symlink needs no change: the store stays one level deep.
- Why a flat store and not `<namespace>/<name>-<version>`: a nested store, the layout of modules and pipelines, needs extra code in offline mode, embedded mode and the symlink name, plus creating the namespace folder. The flat store reuses the lock file scheme and keeps that code small.
- Download (`PluginUpdater.download0`): the unzipped directory is named after the zip file (`<name>-<version>.zip`), so for a namespaced plugin it does not match the store folder name. Nextflow renames the folder when it puts it in the store.

**Plugin ID from the MANIFEST.** `NamespaceDescriptorFinder`, a subclass of pf4j `ManifestPluginDescriptorFinder`, sets the ID when it reads the MANIFEST:

1. `Plugin-Id` as is, when it contains `/` (a zip built with a namespaced `Plugin-Id`);
2. `<Plugin-Namespace>/<Plugin-Id>`, when `Plugin-Namespace` is present;
3. `<namespace>/<Plugin-Id>`, when the store folder name is `<namespace>~...` (a release built without `Plugin-Namespace`);
4. otherwise, the plain `Plugin-Id`.

`LocalPluginManager`, `EmbeddedPluginManager`, `DevPluginManager` and the offline repository (`LocalUpdateRepository`) use this finder. If the request was namespaced and the zip declares a different namespace, Nextflow refuses the plugin. For a plain request, a declared namespace is normal: it is a moved plugin.

**Plugin store.** A moved plugin is stored once. A symlink gives it the other name:

- Online, with the hint, Nextflow stores the plugin in `myorg~nf-foo-<ver>` and links `nf-foo-<ver>` to it.
- If an older Nextflow already stored `nf-foo-<ver>`, Nextflow links `myorg~nf-foo-<ver>` to that folder and does not download again.
- A link always points to a real, complete folder, never to another link.
- Newer Nextflow takes one download lock for both names: the lock of the namespaced ID. Older Nextflow takes the lock of the plain ID.
- If both run at once, older Nextflow can replace the link with a real folder. This is harmless: the store then keeps two copies.
- If the link cannot be created (read-only store, no symlink support), Nextflow logs it and goes on. The plugin still loads.

**ID mapping.** Nextflow keeps a map from plain name to namespaced ID. Sources:

- the `canonicalId` hints of `/dependencies`;
- a zip served for a plain request whose `Plugin-Namespace` gives a namespaced ID;
- in embedded mode, the loaded plugins that declare `Plugin-Namespace`;
- the store links, online and offline.

Nextflow never maps `nf-foo` to `myorg/nf-foo` from names alone, because `nf-foo` can be a different plugin. Store links are trusted because only Nextflow writes them, and only from a hint or a served zip.

Order in a run:

1. Read the store links.
2. Fetch the metadata of the config, `-plugins` and default plugins. The hints arrive here.
3. Map the IDs.
4. Merge the config and default plugins. The config ID wins.
5. Check `NXF_PLUGINS_ALLOWED` against the mapped IDs.
6. Load the plugins. In the dependency loop, map each `Plugin-Dependencies` entry in memory before pf4j resolves it. The MANIFEST file is not changed.

The metadata of a plugin is kept under both IDs, so a lookup by either ID finds it. This also works when the `plugins.json` index has only the plain ID.

The map applies everywhere Nextflow looks up a plugin by ID: `include ... from 'plugin/nf-foo'`, on-demand starts, CLI commands, updates, and the pinned core plugin versions.

Nextflow skips the registry call for pinned plugins already in the store. When a run has a plain ID and a namespaced ID with the same name and no link between them, it makes the call anyway, to learn whether they are one plugin.

**Notice.** When the config or `-plugins` uses the plain name of a moved plugin, Nextflow logs a warning once per plugin: "Plugin `nf-foo` has moved to `myorg/nf-foo` -- update your config". The run goes on.

**Offline.**

- Nextflow reads the map from the store links, one link per version. It skips broken links.
- A store filled only by an older Nextflow has no links, so it gives no map.
- The offline docs say: fill the store online with the Nextflow version used offline, and copy it with links kept (`cp -a`, `rsync -a`).
- `nextflow plugin install` asks the registry for every requested plugin, also for plugins already in the store, so it can add the missing links.

**Registry refuses a namespaced ID.** An old registry answers 400 to a metadata request with a `/` in an ID. A new registry can also refuse one (bad name, old Nextflow). Nextflow shows the registry response as a warning, with a note that the registry may not support namespaced plugins. The run then fails with "Cannot find version", as today. The note is added only to a 400 answer, not to errors such as 429 or 503.

**Several registries.** `registry.url` can list several registries. Nextflow asks each one for the requested IDs and combines the releases of each ID. When two registries serve the same version, the first one in the list wins. A registry that refuses an ID is skipped, and the next one is asked. With a moved plugin:

- One registry moved `nf-foo` to `myorg`, another has plain `nf-foo` and no move (a mirror, or an older registry). A request for `nf-foo` gives one plugin, `myorg/nf-foo`, with the releases of both registries. The list order decides only which registry wins a version they both serve. A release that only the other registry has is stored in `myorg~nf-foo-<ver>` and loads as `myorg/nf-foo`, like an old release of a moved plugin.
- The same two registries, with a request for `myorg/nf-foo`. Only the registry that moved the plugin knows this ID. The other one returns nothing (new registry) or refuses the ID (old registry, shown as a warning). Its plain releases are not used. Nextflow does not also ask for the plain name: that is one more request per registry, only to cover a registry that did not move the plugin.
- Two registries moved `nf-foo` to different namespaces, `org1` and `org2`. The first registry in the list sets the namespace of the plain name: `nf-foo` becomes `org1/nf-foo`. Nextflow logs a warning once. The other plugin is used only when it is requested by its full ID, `org2/nf-foo`. The releases of the two plugins never mix: each registry puts its plain releases under its own `canonicalId`.
- Two registries have plain `nf-foo` for different plugins, and neither moved it. Nextflow mixes their releases, as before this change. Namespaces are how authors avoid this.

**Core plugins.** The hardcoded IDs change to `nextflow-io/...` in the same Nextflow release as the namespace support. Core zips start writing `Plugin-Namespace` only from that release. See the worked example below.

**Dev mode.** `DevPluginManager` uses `NamespaceDescriptorFinder`, but reads the MANIFEST of a local plugin from its build folders. Local and downloaded plugins get the same ID as in a release.

### Registry (summary)

- Reuses its namespace model for modules and pipelines: ID validation with an optional `<namespace>/` prefix, a `namespace_id` column on plugins, publish permission through namespace membership, `%2F`-encoded IDs in API routes, download file name `<name>-<version>.zip`, and a `plugins/` prefix in OCI storage.
- Namespaced plugins are left out of the legacy `plugins.json`. Moved plugins stay in it under their plain name.
- On publish, it checks `Plugin-Id` against the plugin name and `Plugin-Namespace` against the namespace. A moved plugin may publish without `Plugin-Namespace`, for patch releases built from older branches.
- It refuses a namespaced ID when the request's `nextflowVersion` is older than a configured minimum: the first edge release that reads `Plugin-Namespace` (planned: 26.12.0-edge). The docs name the stable release, 27.04. For such a client, `/dependencies` also leaves out namespaced plugins found in the dependency walk, and the releases that depend on one.
- Move: the owner of an unscoped plugin who is a member of the target namespace, or an admin, moves it from the plugin page. The record keeps its ID and gets the namespace. Both IDs resolve to it, and `/dependencies` returns the `canonicalId` hint.
- Revoked namespace: the namespaced ID returns 410 on metadata and download. `/dependencies` leaves it out, so Nextflow reports "Cannot find version", not "revoked". The plain name of a moved plugin keeps working.
- Enforcement is a registry flag (`REGISTRY_PLUGINS_NAMESPACE_ENFORCED`), off by default. When on, unscoped plugins that did not move cannot publish, and new unscoped plugins cannot be claimed.

See the ADR in the registry repository.

### Gradle plugin (summary)

- A new optional setting: `nextflowPlugin { namespace = 'myorg' }`.
- It writes `Plugin-Namespace: myorg` in the MANIFEST. `Plugin-Id` stays `<project.name>`.
- The release request uses the ID `myorg/<project.name>`.
- The local install path is `plugins/myorg~<project.name>-<version>`.
- The zip name stays `<project.name>-<version>.zip`.
- A refused publish shows the registry message, as today.

See the ADR in the `nextflow-plugin-gradle` repository.

### Compatibility

New namespaced plugins:

| Nextflow | Registry | Unscoped plugin | Namespaced plugin |
|---|---|---|---|
| old | old | works | not possible |
| new | old | works | fails: clear message from Nextflow, shown as a warning, then "Cannot find version" |
| old | new | works | fails: the registry message (`nextflowVersion` check) is shown as a warning, then "Cannot find version" |
| new | new | works | works |

| Gradle plugin | Registry | Unscoped plugin | Namespaced plugin |
|---|---|---|---|
| old | any | works | not possible (no `namespace` setting) |
| new | old | works | fails: the registry refuses the ID |
| new | new | works | works |

Moved plugins:

| Nextflow | Asks for `nf-foo` (moved) | Asks for `myorg/nf-foo` |
|---|---|---|
| reads `plugins.json` only | releases as today, loaded as `nf-foo` | not possible |
| `/dependencies`, no namespace support | all releases, loaded as `nf-foo` | not possible |
| namespace support | loaded as `myorg/nf-foo`, with the notice | loaded as `myorg/nf-foo` |

What older Nextflow does not get:

- new namespaced plugins, which have no plain name;
- releases of a moved plugin that depend on a namespaced ID. The registry leaves them out of `/dependencies` for older Nextflow, so it picks the newest release it can load. It also logs a warning at publish.

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

Required change in xpack, before the first build with `Plugin-Namespace`: use the name without the namespace as the product code. `core/xpack` is bundled into each plugin, so every `seqera/xpack-*` plugin must be rebuilt with this fix.

```groovy
product = wrapper.descriptor.pluginId.tokenize('/').last()
```

### Rollout

1. Registry and Gradle plugin releases, with enforcement off.
2. Nextflow release with namespace support, including the core plugin ID change.
3. Set the registry minimum Nextflow version for namespaced IDs to that release.
4. Move the core plugins to `nextflow-io` and the Seqera plugins to `seqera`. Apply the xpack product-code change before the first build with `Plugin-Namespace`.
5. Tell plugin authors when enforcement will be turned on, with a how-to for claiming a namespace and moving a plugin.
6. Turn enforcement on.

## Prototype results

A prototype in Nextflow, on top of nextflow-io/nextflow#8082, tested the Nextflow part of B5. The registry hint was simulated with a map.

- pf4j loads a zip as `<Plugin-Namespace>/<Plugin-Id>` through the `protected` descriptor hooks. Released Nextflow 26.10.4 loads a zip with `Plugin-Namespace` with no warning.
- With the hint, a moved plugin loads once, whatever ID the config, `include`, a dependency or an on-demand start uses.
- Without the hint:
  - a plain ID requested after the namespaced one fails with "There is an already loaded plugin", or starts the plugin twice for a release built without `Plugin-Namespace`;
  - a plain dependency on a moved plugin that is not loaded yet fails with `DependenciesNotFoundException`.
- Offline, the store links give the map. Released Nextflow 26.10.4 works offline on a store with links, for pinned versions and version ranges.

Other store layouts considered for the offline map: a marker file inside the plugin folder (safe, but older and newer Nextflow keep two copies) and a separate store for newer Nextflow. Links keep one copy and need no change in older Nextflow.

## Rationale

- Plugins and modules live in the same registry. One naming and ownership model is simpler for users and for the registry.
- Putting the namespace in the ID (instead of only in the registry) lets a pipeline choose the exact plugin, and lets two plugins with the same name run together.
- One record and one zip per release remove the cost that B2 and B3 put on authors and on the registry.
- `Plugin-Id` stays plain, so older Nextflow loads every release of a moved plugin without change. Because of this, the registry can make moving required.
- The price is in Nextflow: the ID mapping, which needs the registry hint, and the saved map for offline runs.

## Consequences

- Pipelines that use new namespaced plugins need a minimum Nextflow version: the first release that reads `Plugin-Namespace`. They should set `manifest.nextflowVersion`.
- With an old registry, only the namespaced plugin fails. Nextflow then retries with one plugin per request, so the unscoped plugins of the same run still resolve.
- Newer Nextflow needs the `canonicalId` hint before it loads a plugin. A registry or index without it can load a moved plugin twice.
- Offline runs need a store filled online by a Nextflow version with namespace support.
- Older Nextflow sharing a store with newer Nextflow can fail offline with a version range (see Compatibility).
- The ID shown by Nextflow (`myorg/nf-foo`) differs from `Plugin-Id` in the MANIFEST (`nf-foo`).
- Authors of moved plugins keep plain dependencies while they support older Nextflow.
- Once enforcement is on, no new unscoped plugins are added.

## Open decisions

- When to turn enforcement on, and how long after the Nextflow release.
- How an admin allows a single release that enforcement would refuse.
- Whether an admin can move a plugin whose author does not respond.
- Docs: ID format, the `Plugin-Namespace` attribute, store layout, minimum Nextflow version, the migration guide, and offline store preparation.
