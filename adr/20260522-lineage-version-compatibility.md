# Lineage record version compatibility

- Authors: Jorge Ejarque
- Status: draft
- Deciders: Paolo Di Tommaso, Ben Sherman
- Date: 2026-05-22
- Tags: lineage, serialization, versioning

## Summary

Define how lineage records evolve across Nextflow releases: how versions are stamped, when they are bumped, what changes are compatible, and how incompatible changes are migrated on read.

## Problem Statement

The `nf-lineage` module persists workflow execution metadata as JSON records under a `LinModel.VERSION` identifier (today `lineage/v1beta1`). Each record is shaped as:

```json
{ "version": "lineage/v1beta1", "kind": "WorkflowRun", "spec": { ... } }
```

`LinTypeAdapterFactory` writes the current version into every record and rejects any record whose version does not match exactly. Today there is only one version, so this works, but it is not a versioning strategy — it is a fail-fast guard.

As the lineage model evolves we need to answer three questions:

1. **When does the version string change?** Adding a nullable field is conceptually different from renaming or restructuring a field; they should not both force the same kind of code churn.
2. **How does newer Nextflow read records written by older Nextflow?** The store is append-only, but binaries upgrade; old records persist and must remain readable.
3. **How is a JSON schema associated with each version?** External consumers (Seqera Platform, third-party tooling) need a stable schema artifact per version.

Without an explicit contract, any change to the model risks breaking existing lineage stores or downstream tooling that consumed an older schema.

## Goals

- Allow older lineage records to remain readable after a Nextflow upgrade.
- Distinguish *compatible* deltas (no migration required) from *incompatible* deltas (explicit migration required).
- Make every version string map 1:1 to a published JSON schema document.
- Keep the in-memory model singular at any point in time so that downstream consumers (`LinObserver`, `LinUtils`, `LinExtensionImpl`, `LinDagRenderer`, `LinPropertyValidator`) do not need to branch on version.
- Surface unknown versions with an actionable error rather than silent failure or silent coercion.

## Non-goals

- **Backward read across binaries.** Older Nextflow installs are not required to read records written by newer Nextflow. Tolerant parsing within a MAJOR happens to provide this for compatible deltas, but it is not a guarantee across MAJOR boundaries.
- **In-place store migration.** Records on disk are never rewritten on upgrade. All version handling happens on read.
- **Per-version in-memory types.** Only the current MAJOR's model classes live in the codebase. Older records are normalised to the current shape before being handed to gson.
- **Runtime schema validation.** The runtime trusts gson. Published schemas are documentation artifacts for external consumers, not enforcement.

## Decision

Adopt **semver MAJOR.MINOR version strings** for lineage records, a **tolerant-reader contract** within a MAJOR, and a **JSON-tree migration mechanism** for cross-MAJOR records. The current `lineage/v1beta1` string is grandfathered as a legacy tier label.

## Version string scheme

Lineage records carry a `version` field with one of:

- `lineage/v1beta1` — the legacy tier label (current production format). Treated by the reader as its own opaque MAJOR.
- `lineage/<MAJOR>.<MINOR>` for all future versions, e.g. `lineage/1.0`, `lineage/1.1`, `lineage/2.0`. PATCH is omitted — bug fixes never change the wire format.

Each distinct version string maps to exactly one JSON schema document, published in the `nextflow-io/schemas` repository.

## Compatibility contract

The MAJOR component of the version string is the *compatibility tier*. Within a MAJOR, the JSON reader is required to tolerate the deltas listed below. Across MAJORs, an explicit migrator is required.

| Change | Same MAJOR (MINOR bump) | Different MAJOR |
| --- | --- | --- |
| Add an optional / nullable field with a safe default | Allowed. Old reader ignores it; new reader defaults it. | n/a |
| Add a required field | Not allowed (effectively required = no safe default). | Allowed with migrator that synthesises it. |
| Remove a field | Allowed if the field is no longer consulted by any reader. Use a two-step deprecation: stop writing in MINOR N, stop reading in MINOR N+1. | Allowed with migrator. |
| Rename a field | Not allowed. | Allowed with migrator (rewrite key). |
| Retype a field (e.g. `String` -> `List<String>`) | Not allowed. | Allowed with migrator. |
| Restructure (move field between objects, split, merge) | Not allowed. | Allowed with migrator. |
| Add a new `kind` (new `LinSerializable` subtype) | Allowed. Older readers will fail when they encounter it; that is expected. | n/a |
| Remove an existing `kind` | Not allowed within MAJOR. | Allowed with migrator. |

The rule of thumb: if gson with `withSerializeNulls(true)` and default field initialisers can round-trip the change in both directions without losing information, it is a MINOR bump. Otherwise it is a MAJOR bump.

## Read dispatch

`LinTypeAdapterFactory.read()` follows this sequence:

1. Parse the JSON object.
2. Read the `version` field. If missing, throw.
3. Extract the MAJOR component:
   - `lineage/v1beta1` -> `v1beta1` (opaque legacy major).
   - `lineage/<n>.<m>` -> `<n>` (integer prefix).
4. Look up MAJOR in a hard-coded `KNOWN_MAJORS` set. If absent, throw a `JsonParseException` that names the offending version and lists the known set, hinting that the operator may need to upgrade Nextflow.
5. If the record's version is not equal to `CURRENT_VERSION`, run the migrator chain to rewrite the JSON tree into the current MAJOR's shape. *(No migrators exist yet; this step is a no-op pending the first MAJOR bump.)*
6. Delegate to the gson `RuntimeTypeAdapterFactory` to produce the typed object.

Two pre-`spec` variants of `lineage/v1beta1` records exist in the wild:

1. **Modern layout**: payload under `spec`.
2. **Pre-spec legacy layout**: payload at root, no `spec` wrapper.

The existing inline fallback for variant 2 (currently at the tail of `LinTypeAdapterFactory.read()`) remains. It is documented in a class-level comment listing both variants so future maintainers do not remove it as dead code.

## Migration mechanism

Cross-MAJOR migration is performed by **JSON-tree transforms**. A migrator is a function `JsonObject -> JsonObject` that converts a record written under one MAJOR into a record consumable by the next MAJOR's gson adapter. Migrators chain: `v1beta1 -> 1.0 -> 2.0 -> ...`. When the third MAJOR ships, only one new edge (`2.0 -> 3.0`) needs to be written; older paths compose automatically.

The transforms are JSON-to-JSON, so each migrator is unit-testable in isolation against golden input and expected output files. Only one set of in-memory model classes (the current MAJOR) is needed; old model classes are not retained.

The first migrator is not introduced by this ADR — none is required while only `lineage/v1beta1` exists. The first concrete migrator will accompany the v1beta1 -> 1.0 graduation. Until then, `LinTypeAdapterFactory` carries the dispatch shape (KNOWN_MAJORS check) but the migration call is a documented placeholder.

## Considered Options

### Version string shape

- **Adopted: semver MAJOR.MINOR (e.g. `lineage/1.0`)** — each string maps 1:1 to a schema file; external consumers always know which schema validates a given record.
- **K8s-style (string = ABI tier only, schema rev tracked separately, e.g. `lineage/v1beta1` stable while schema grows underneath)** — sound in a live-API setting where the apiserver ships its OpenAPI alongside the binary, but awkward for an on-disk format that may be inspected by tools long after the writing Nextflow has been retired.
- **Semver with PATCH** — adds a string component that, by definition, would never change the wire format. Reader logic would never key on it.

### Reader compatibility model

- **Adopted: forward-read plus migrate-on-read.** The codebase pins to the current MAJOR; old records are normalised to that shape before reaching downstream consumers.
- **Forward-read without migration** — exposes records as the version they were written in. Pushes per-version branching into every consumer.
- **Multi-version in-memory types** — first-class K8s-style. Larger refactor than the lineage module justifies today.
- **Forward + backward read across binaries** — would require a frozen field-subset discipline and downgrade conversions. Not justified by the use case.

### Migration mechanism

- **Adopted: JSON-tree transform chain.** One transform per edge. Generalises the existing legacy-layout fallback at the tail of `LinTypeAdapterFactory.read()`. Only one model package live.
- **Per-version adapter to current model** — each old MAJOR has its own adapter targeting the current model. No chain. Rewrites grow over time.
- **Per-version model packages plus converters** — strongest type safety, biggest maintenance cost; every model class duplicated per MAJOR.
- **Declarative JSON Patch rules** — pretty for trivial cases, opaque for complex transforms, requires a patch library.

### Where migrators live

- **Adopted (for now): inline in `LinTypeAdapterFactory`, registry deferred.** No second migrator exists yet, so the registry abstraction would be premature. The first migrator (v1beta1 -> 1.0) is the trigger to extract a `serde/migration/` package and a `LinMigrationRegistry`.
- **Introduce `serde/migration/` immediately** — cleaner separation but no consumers today.

### Model package layout

- **Adopted: keep `nextflow.lineage.model.v1beta1.*` until graduation.** The name is true today. Graduation (v1beta1 -> 1.0) is the natural moment to rename, and that change rides with the model-class cleanup that graduating from beta typically requires.
- **Flatten to `nextflow.lineage.model.*` now** — pre-emptive churn across every consumer.

### JSON schema generation and hosting

- **Adopted: hosted in `nextflow-io/schemas`, no runtime validation.** Generation mechanism (auto-generated from classes, or hand-written) is an orthogonal decision tracked in that repository. The runtime trusts gson.
- **Bundle schemas in `modules/nf-lineage/src/main/resources/`** — keeps them with the code but couples external consumers to Nextflow's release cadence.
- **Validate every record on read against its schema** — duplicates gson's parse work, pays per-record cost.

### Unknown-version handling

- **Adopted: hard fail with a detailed message** naming the offending version and the known set; suggests upgrading Nextflow.
- **Best-effort tolerant parse** — risks silent coercion of unknown semantics into current shapes.
- **Warn and skip (return `null`)** — would propagate as silent data loss to every downstream consumer that currently assumes records read successfully.

## Consequences

**Positive:**

- Compatible deltas (additive nullable fields, deprecated field removal) cost a MINOR bump and a schema-file copy. No migration code, no consumer changes.
- Old records remain queryable after an upgrade. A single lineage store can carry records from multiple Nextflow versions.
- The 1:1 mapping between version string and JSON schema means any record can be validated against its declared schema using only the string in the record. External tooling does not need to know which Nextflow wrote the record.
- The migration call site is established now; the first cross-MAJOR change adds one migrator and the registry extraction, without redesigning the read path.
- Unknown-version errors are immediately diagnosable: the operator is told exactly which version was found and which are known.

**Negative:**

- Every MINOR bump (i.e. every compatible field add or removal that ships in a release) publishes a new schema file. Over many releases this accumulates; archival hygiene in `nextflow-io/schemas` becomes a real concern.
- The runtime trusts gson and does not enforce that records on disk actually conform to their declared schema. A record that lies about its version will be parsed against the wrong adapter and may produce subtly wrong objects rather than failing loudly.
- The MAJOR-only dispatch means that any future plugin that wants to teach the reader about an additional MAJOR has no SPI today — `KNOWN_MAJORS` is a static set. This is recorded in Future Work.

**Neutral:**

- `lineage/v1beta1` becomes "the legacy tier" rather than "the current version" at the moment graduation lands. Until then, the package name `model.v1beta1` and the wire-format identifier match.
- The contract above is enforced by code review, not by tooling. A reviewer who waves through a renamed field within a MINOR will break the contract silently.

## Future Work

- **Extract `serde/migration/` package and `LinMigrationRegistry`.** Triggered by the first cross-MAJOR migrator (e.g. v1beta1 -> 1.0 graduation). Each migrator becomes a class implementing `LinMigrator { String fromVersion; String toVersion; JsonObject migrate(JsonObject); }`. The registry resolves the migration chain from a record's stamped version to the current version.
- **Plugin-contributed migrators.** Expose `LinMigrationRegistry` as a plugin SPI so that a Nextflow plugin can ship a migrator for a MAJOR that the core binary does not natively understand. This effectively lets an older Nextflow read records written by a newer one as long as the appropriate migrator plugin is installed — the K8s version-skew benefit, but plugin-distributed.
- **Schema generation tooling in `nextflow-io/schemas`.** Decide between auto-generation from the Groovy classes (via `victools/jsonschema-generator` or similar) and hand-written schemas, including how PR review surfaces drift. Out of scope for this ADR.
- **Schema-validation mode for CI.** Optional runtime validation gated behind a debug flag, useful for catching drift between class and schema during development without paying the cost in production reads.
- **Graduation path for `v1beta1` -> `1.0`.** Concrete plan for the first MAJOR bump, including the model-class cleanup that graduation justifies and the introduction of the first real migrator.

## Links

- Related module: `modules/nf-lineage/src/main/nextflow/lineage/`
- Affected file: `modules/nf-lineage/src/main/nextflow/lineage/serde/LinTypeAdapterFactory.groovy`
