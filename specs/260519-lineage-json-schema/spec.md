# Feature Specification: JSON Schema generation for lineage model v1beta1

**Feature Branch**: `lineage-json-schema`
**Created**: 2026-05-19
**Status**: Draft
**Input**: User request — provide a way to produce a JSON Schema describing the JSON documents emitted by `LinEncoder` for the `nextflow.lineage.model.v1beta1` classes, suitable for checking into the repository as a versioned artifact.

## Motivation

`modules/nf-lineage/src/main/nextflow/lineage/serde/LinEncoder` serializes `LinSerializable` objects to JSON via Gson with a `RuntimeTypeAdapterFactory`. The on-the-wire shape is an envelope:

```json
{ "version": "lineage/v1beta1", "kind": "<SimpleName>", "spec": { /* subtype fields */ } }
```

There is no machine-readable description of that shape today. Consumers (validators, external tooling, future contributors editing the model) have to read the Groovy classes to understand the schema. A JSON Schema would:

- Provide a contract for the `seqera://`-style lineage JSON documents.
- Make breaking changes to the model visible in PR diffs once the generated schema is committed.
- Enable third-party validation without depending on the JVM.

## Goals

- Produce a single JSON Schema document covering all six `LinTypeAdapterFactory`-registered subtypes (`WorkflowRun`, `WorkflowOutput`, `Workflow`, `TaskRun`, `TaskOutput`, `FileOutput`) as a top-level `oneOf` over the `{version, kind, spec}` envelope.
- Generate the schema reproducibly from the compiled Groovy classes via a Gradle task.
- Check the resulting schema into `modules/nf-lineage/src/resources/schema/` so diffs surface in PRs.
- Keep all schema-generation dependencies off `nf-lineage`'s runtime/compile classpath.

## Non-Goals

- No automatic CI enforcement that the checked-in schema matches the current model (manual regeneration).
- No runtime validation of lineage JSON against the schema inside Nextflow.
- No schema for non-`LinSerializable` model classes (`Checksum`, `DataPath`, `Parameter`, `LinModel`) as top-level branches — they appear only as nested schemas referenced from the registered subtypes.
- No changes to the JSON wire format produced by `LinEncoder`.
- No changes to how `LinTypeAdapterFactory` works at runtime (single optional comment only).

## User Scenarios & Testing

### User Story 1 — Generate and commit the v1beta1 schema (Priority: P1)

A maintainer wants the first version of the schema produced and added to the repo.

**Why this priority**: Without this, the feature has zero user-visible output.

**Independent Test**: Run `./gradlew :nf-lineage:generateLineageSchema`. Confirm `modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json` exists, parses as valid JSON, declares `$schema: https://json-schema.org/draft/2020-12/schema`, and has a top-level `oneOf` of length 6 — one branch per registered subtype.

**Acceptance Scenarios**:

1. **Given** a clean checkout, **When** the task is run, **Then** the schema file is created with exactly six `oneOf` branches, each having `properties.version.const = "lineage/v1beta1"`, `properties.kind.const = "<SubtypeSimpleName>"`, and a `properties.spec` schema describing the subtype's fields.
2. **Given** the schema file already exists, **When** the task is run again with no model changes, **Then** the file content is byte-identical (or differs only in stable, deterministic ordering).
3. **Given** a sample lineage JSON document produced by `LinEncoder.encode(...)` for any subtype, **When** validated against the generated schema using a third-party validator, **Then** validation succeeds.

### User Story 2 — Refresh schema after a model change (Priority: P2)

A developer adds a new field to `WorkflowRun` (or any subtype) and needs to update the committed schema.

**Why this priority**: The on-going maintenance workflow; defines the manual-only update story chosen during brainstorming.

**Independent Test**: Add a field to one of the model classes, run the task, observe the schema file diff includes the new property under the corresponding `spec` schema.

**Acceptance Scenarios**:

1. **Given** a new `String` field on a registered subtype, **When** the task is run, **Then** the corresponding `oneOf` branch's `spec.properties` includes that field with `type: string`.
2. **Given** a new field of type `OffsetDateTime`, **When** the task is run, **Then** the field appears with `{type: string, format: date-time}`.
3. **Given** a new field of type `java.nio.file.Path`, **When** the task is run, **Then** the field appears with `{type: string}`.

### User Story 3 — Add a new `LinSerializable` subtype (Priority: P3)

A developer adds a new model class and registers it in `LinTypeAdapterFactory`.

**Why this priority**: The maintenance path most likely to drift.

**Acceptance Scenarios**:

1. **Given** a new `LinSerializable` class registered via `registerSubtype(...)`, **When** the developer also appends its fully-qualified name to the `subtypes` list in `modules/nf-lineage/build.gradle` and re-runs the task, **Then** the schema's `oneOf` grows to include a new branch for that subtype.
2. **Given** a new subtype registered in the factory but NOT added to the `subtypes` list in `build.gradle`, **When** the task is run, **Then** the schema is silently incomplete (this is the known drift trade-off, mitigated by an in-source comment — not enforced by CI).

## Design

### Component 1 — `buildSrc` module hosting the task class

New top-level `buildSrc/` directory (currently absent from the repo; the existing `buildSrc/build/` tree is just Gradle's empty-buildSrc output dir).

**`buildSrc/build.gradle`** (new):

```groovy
plugins { id 'groovy' }
repositories { mavenCentral() }
dependencies {
    implementation gradleApi()
    implementation localGroovy()
    implementation 'com.github.victools:jsonschema-generator:4.36.0'
}
```

**`buildSrc/src/main/groovy/nextflow/gradle/GenerateLineageSchemaTask.groovy`** (new): a `DefaultTask` with three inputs and one `@TaskAction`:

- `@InputFiles FileCollection classpath` — the compiled `nf-lineage` classes plus its runtime classpath
- `@Input List<String> subtypes` — fully-qualified class names to expose as top-level `oneOf` branches
- `@OutputFile File outputFile` — destination path for the schema JSON

`@TaskAction generate()` performs:

1. Build a `URLClassLoader` over `classpath`, parented to the task's own class loader. This lets the loaded model classes resolve `LinSerializable` (and JDK types) consistently.
2. Construct a `SchemaGenerator` from victools with `SchemaVersion.DRAFT_2020_12` and `OptionPreset.PLAIN_JSON`. Register two custom type definitions:
   - `OffsetDateTime` → `{"type": "string", "format": "date-time"}`
   - `java.nio.file.Path` → `{"type": "string"}`
3. Load `nextflow.lineage.model.v1beta1.LinModel` via the URLClassLoader and read its `VERSION` static field reflectively — single source of truth for the envelope's `version` const, no duplication of the literal `"lineage/v1beta1"`.
4. For each FQN in `subtypes`: load the class via the URLClassLoader, call `generator.generateSchema(cls)`, and wrap the result in an envelope object:
   ```json
   {
     "type": "object",
     "properties": {
       "version": { "const": "<LinModel.VERSION>" },
       "kind":    { "const": "<SimpleName>" },
       "spec":    <subtypeSchema>
     },
     "required": ["version", "kind", "spec"],
     "additionalProperties": false
   }
   ```
5. Compose the root document:
   ```json
   {
     "$schema": "https://json-schema.org/draft/2020-12/schema",
     "title": "Lineage v1beta1",
     "oneOf": [ /* one envelope per subtype, in input order */ ]
   }
   ```
6. Serialize with a deterministic pretty-printer (victools returns a Jackson `ObjectNode`; use Jackson's `writerWithDefaultPrettyPrinter()` to render). Write bytes to `outputFile`, creating parent dirs.

### Component 2 — Task registration in `modules/nf-lineage/build.gradle`

Append:

```groovy
import nextflow.gradle.GenerateLineageSchemaTask

tasks.register('generateLineageSchema', GenerateLineageSchemaTask) {
    description = 'Generate JSON Schema for the lineage model v1beta1'
    group = 'documentation'
    dependsOn compileGroovy
    classpath = sourceSets.main.runtimeClasspath + sourceSets.main.output
    // Keep this list in sync with LinTypeAdapterFactory.registerSubtype(...) calls
    // and re-run this task whenever the model changes.
    subtypes = [
        'nextflow.lineage.model.v1beta1.WorkflowRun',
        'nextflow.lineage.model.v1beta1.WorkflowOutput',
        'nextflow.lineage.model.v1beta1.Workflow',
        'nextflow.lineage.model.v1beta1.TaskRun',
        'nextflow.lineage.model.v1beta1.TaskOutput',
        'nextflow.lineage.model.v1beta1.FileOutput',
    ]
    outputFile = file('src/resources/schema/lineage-v1beta1.schema.json')
}
```

### Component 3 — In-source sync comment

`modules/nf-lineage/src/main/nextflow/lineage/serde/LinTypeAdapterFactory.groovy`, immediately above the `registerSubtype` chain in the constructor:

```groovy
// When adding or removing a subtype, also update the `subtypes` list in
// modules/nf-lineage/build.gradle (task `generateLineageSchema`) and re-run
// `./gradlew :nf-lineage:generateLineageSchema` to refresh
// src/resources/schema/lineage-v1beta1.schema.json.
```

This is the only modification to production source code. No behavior change.

### Component 4 — First generated schema artifact

`modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json` (new). Produced by running the task once and committed alongside the other changes.

### Testing

This is a tooling change with no runtime behavior. Coverage comes from running the task and inspecting the artifact:

- **Manual verification during implementation** (mandatory): run `./gradlew :nf-lineage:generateLineageSchema`, confirm the file is created, parse it as JSON, verify it has six `oneOf` branches each with the expected `kind` const value, and validate at least one known-good lineage JSON document (taken from an existing `LinEncoder` test fixture) against the generated schema using a standard JSON Schema 2020-12 validator.
- **No buildSrc unit test**: the task is mostly orchestration around victools. Setting up a Gradle `ProjectBuilder` harness with the correct classpath inside `buildSrc` would be more code than the task itself, and the manual verification above covers the same surface. Deferred unless the task accrues real branching logic.
- **No production-code test changes** — `LinTypeAdapterFactory` is only gaining a comment.

## Risks

- **victools schema fidelity for Groovy classes**: `@CompileStatic` + `@Canonical` produces plain JavaBean-style properties, which victools reads via the default field-introspection mode. If a field type isn't covered by the two custom definitions and isn't a primitive, `String`, `List`, `Map`, or a `LinSerializable`-related class, victools may produce a generic `{type: object}`. The current model only uses `String`, `Long`/`long`, `Map`, `List<T>`, `OffsetDateTime`, `Path`, `Checksum`, `DataPath`, `Parameter`, and nested `LinSerializable`s. Verified covered during implementation; new exotic types in the future may need new custom definitions.
- **Drift between `subtypes` list and `LinTypeAdapterFactory`**: deliberately accepted, mitigated by in-source comment. No CI enforcement per the brainstorming decision.
- **`$schema` URI choice**: draft 2020-12 is the current stable. `const` keyword used for discriminator values is supported from draft-06 onward, well within 2020-12.
- **Untyped `Map config` / `Map metadata` in `WorkflowRun`**: schema will describe them as `{type: object}` with no `properties`/`additionalProperties` constraints. This matches the intentional opacity; documented here so it's not flagged later as a generator bug.
- **Determinism**: victools and Jackson both produce stable output for the same input, but property ordering inside generated subtype schemas relies on victools' field-order behavior. If non-determinism shows up across JVMs, address by post-processing (sort `properties` keys) — not anticipated for now.

## Out of Scope

- A `checkSchema` task that fails CI on drift.
- Wiring `generateLineageSchema` into `check`, `build`, or `publish`.
- Schema generation for any model package other than `v1beta1`.
- Auto-discovery of subtypes (considered and rejected during brainstorming in favor of an explicit list with a sync comment).
- Validation of incoming JSON against the schema at runtime inside Nextflow.
- A `module-spec`–style ADR; if the v1beta1 schema becomes part of an external contract, that's a separate ADR follow-up.