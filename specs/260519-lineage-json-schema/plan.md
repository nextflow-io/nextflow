# Lineage v1beta1 JSON Schema — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a JSON Schema (draft 2020-12) describing the lineage v1beta1 wire format emitted by `LinEncoder`, generated from the compiled Groovy classes by a Gradle task, and check the resulting `lineage-v1beta1.schema.json` into `modules/nf-lineage/src/resources/schema/`.

**Architecture:** A new top-level `buildSrc/` Groovy build module hosts a custom `DefaultTask` (`GenerateLineageSchemaTask`) that depends only on `com.github.victools:jsonschema-generator` on the build-script classpath. The task takes a classpath, a hardcoded list of subtype FQNs, and an output file; it loads each subtype via a `URLClassLoader` and runs victools to produce per-subtype schemas, then wraps each in a `{version, kind, spec}` envelope and emits a single `oneOf` root document. `modules/nf-lineage/build.gradle` registers a `generateLineageSchema` task instance; no runtime dep is added to `nf-lineage`. The only production-source change is a three-line maintainer comment in `LinTypeAdapterFactory`.

**Tech Stack:** Gradle (Groovy DSL) + buildSrc, Groovy 4, victools jsonschema-generator 4.36.0 (Jackson-databind transitively). Module under work: `modules/nf-lineage` and new `buildSrc/`.

---

## File Structure

**Created:**
- `buildSrc/build.gradle` — declares groovy plugin, mavenCentral, and the victools dep on the build-script classpath only.
- `buildSrc/src/main/groovy/nextflow/gradle/GenerateLineageSchemaTask.groovy` — the `DefaultTask` subclass.
- `modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json` — generated artifact, committed.

**Modified:**
- `modules/nf-lineage/build.gradle` — register the `generateLineageSchema` task with hardcoded `subtypes` list and sync-warning comment.
- `modules/nf-lineage/src/main/nextflow/lineage/serde/LinTypeAdapterFactory.groovy` — three-line maintainer comment above the `registerSubtype` chain. No behavior change.

---

## Task 1: Scaffold `buildSrc` with the victools dependency

**Files:**
- Create: `buildSrc/build.gradle`

- [ ] **Step 1: Create `buildSrc/build.gradle`**

Create the file with this exact content:

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
plugins {
    id 'groovy'
}

repositories {
    mavenCentral()
}

dependencies {
    implementation gradleApi()
    implementation localGroovy()
    implementation 'com.github.victools:jsonschema-generator:4.36.0'
}
```

Note: `buildSrc/` already has a `build/` output directory from a prior empty buildSrc Gradle creates by default. Only the `build.gradle` file is new; do not delete the `build/` directory.

- [ ] **Step 2: Verify buildSrc compiles cleanly**

Run:
```
./gradlew help
```

Expected: succeeds with no warnings about buildSrc. Gradle compiles buildSrc once at the start of the build; if there's a syntax error in `buildSrc/build.gradle` or missing repository, this command will fail with a buildSrc compilation error. Confirm it does not.

- [ ] **Step 3: Commit**

```
git add buildSrc/build.gradle
git commit -s -m "Add buildSrc with victools jsonschema-generator dep"
```

---

## Task 2: Implement `GenerateLineageSchemaTask` — class loading + subtype enumeration

This task establishes the skeleton (inputs, classloader, subtype loading, write empty oneOf). Schema generation follows in Task 3.

**Files:**
- Create: `buildSrc/src/main/groovy/nextflow/gradle/GenerateLineageSchemaTask.groovy`

- [ ] **Step 1: Create the task class**

Create the file with this exact content:

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.gradle

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.ObjectNode
import groovy.transform.CompileStatic
import org.gradle.api.DefaultTask
import org.gradle.api.file.FileCollection
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.InputFiles
import org.gradle.api.tasks.OutputFile
import org.gradle.api.tasks.TaskAction

/**
 * Generates a JSON Schema (draft 2020-12) describing the JSON documents emitted
 * by {@code nextflow.lineage.serde.LinEncoder} for v1beta1 model classes.
 *
 * The task does not depend on the nf-lineage runtime; it loads compiled model
 * classes via a URLClassLoader and uses victools jsonschema-generator to derive
 * the per-subtype schemas. Each subtype is wrapped in a {version, kind, spec}
 * envelope matching {@code LinTypeAdapterFactory.write(...)}.
 */
@CompileStatic
class GenerateLineageSchemaTask extends DefaultTask {

    @InputFiles
    FileCollection classpath

    @Input
    List<String> subtypes

    @OutputFile
    File outputFile

    @TaskAction
    void generate() {
        final loader = buildClassLoader()
        final version = readLineageVersion(loader)
        final mapper = new ObjectMapper()
        final root = mapper.createObjectNode()
        root.put('$schema', 'https://json-schema.org/draft/2020-12/schema')
        root.put('title', 'Lineage v1beta1')
        final oneOf = root.putArray('oneOf')

        subtypes.each { String fqn ->
            final cls = loader.loadClass(fqn)
            final spec = generateSubtypeSchema(cls, mapper)
            oneOf.add(wrapEnvelope(cls.simpleName, version, spec, mapper))
        }

        outputFile.parentFile.mkdirs()
        mapper.writerWithDefaultPrettyPrinter().writeValue(outputFile, root)
        logger.lifecycle("Wrote lineage schema with ${subtypes.size()} subtypes to ${outputFile}")
    }

    private URLClassLoader buildClassLoader() {
        final urls = classpath.files.collect { it.toURI().toURL() } as URL[]
        return new URLClassLoader(urls, getClass().classLoader)
    }

    private String readLineageVersion(ClassLoader loader) {
        final linModel = loader.loadClass('nextflow.lineage.model.v1beta1.LinModel')
        return linModel.getField('VERSION').get(null) as String
    }

    // Filled in by Task 3
    private ObjectNode generateSubtypeSchema(Class<?> cls, ObjectMapper mapper) {
        return mapper.createObjectNode()
    }

    private ObjectNode wrapEnvelope(String kind, String version, ObjectNode spec, ObjectMapper mapper) {
        final env = mapper.createObjectNode()
        env.put('type', 'object')
        final props = env.putObject('properties')
        props.putObject('version').put('const', version)
        props.putObject('kind').put('const', kind)
        props.set('spec', spec)
        env.putArray('required').with { it.add('version'); it.add('kind'); it.add('spec'); return it }
        env.put('additionalProperties', false)
        return env
    }
}
```

- [ ] **Step 2: Verify buildSrc still compiles**

Run:
```
./gradlew help
```

Expected: succeeds. The task isn't wired to any project yet, so it just has to compile.

- [ ] **Step 3: Commit**

```
git add buildSrc/src/main/groovy/nextflow/gradle/GenerateLineageSchemaTask.groovy
git commit -s -m "Add GenerateLineageSchemaTask skeleton in buildSrc"
```

---

## Task 3: Implement victools subtype schema generation with custom type defs

Replace the placeholder `generateSubtypeSchema` from Task 2 with the real implementation that uses victools and registers custom type definitions for `OffsetDateTime` and `java.nio.file.Path`.

**Files:**
- Modify: `buildSrc/src/main/groovy/nextflow/gradle/GenerateLineageSchemaTask.groovy`

- [ ] **Step 1: Add imports for victools**

Add these imports near the top of the file, alongside the existing imports:

```groovy
import com.github.victools.jsonschema.generator.CustomDefinition
import com.github.victools.jsonschema.generator.OptionPreset
import com.github.victools.jsonschema.generator.SchemaGenerator
import com.github.victools.jsonschema.generator.SchemaGeneratorConfigBuilder
import com.github.victools.jsonschema.generator.SchemaVersion
import java.nio.file.Path
import java.time.OffsetDateTime
```

- [ ] **Step 2: Replace `generateSubtypeSchema` with the real implementation**

Find this placeholder in the file:

```groovy
    // Filled in by Task 3
    private ObjectNode generateSubtypeSchema(Class<?> cls, ObjectMapper mapper) {
        return mapper.createObjectNode()
    }
```

Replace it with:

```groovy
    private ObjectNode generateSubtypeSchema(Class<?> cls, ObjectMapper mapper) {
        final builder = new SchemaGeneratorConfigBuilder(SchemaVersion.DRAFT_2020_12, OptionPreset.PLAIN_JSON)
        builder.forTypesInGeneral().withCustomDefinitionProvider { javaType, context ->
            final erased = javaType.erasedType
            if (erased == OffsetDateTime) {
                final node = mapper.createObjectNode()
                node.put('type', 'string')
                node.put('format', 'date-time')
                return new CustomDefinition(node)
            }
            if (Path.isAssignableFrom(erased)) {
                final node = mapper.createObjectNode()
                node.put('type', 'string')
                return new CustomDefinition(node)
            }
            return null
        }
        final generator = new SchemaGenerator(builder.build())
        return generator.generateSchema(cls) as ObjectNode
    }
```

- [ ] **Step 3: Verify buildSrc compiles**

Run:
```
./gradlew help
```

Expected: succeeds. Confirms the victools imports resolve from the dep declared in Task 1 and the closure-based `CustomDefinitionProviderV2` SAM conversion compiles under Groovy.

- [ ] **Step 4: Commit**

```
git add buildSrc/src/main/groovy/nextflow/gradle/GenerateLineageSchemaTask.groovy
git commit -s -m "Implement victools-based subtype schema generation"
```

---

## Task 4: Register `generateLineageSchema` in `modules/nf-lineage/build.gradle`

**Files:**
- Modify: `modules/nf-lineage/build.gradle`

- [ ] **Step 1: Append task registration**

Open `modules/nf-lineage/build.gradle`. Currently the file ends with the `dependencies { ... }` block. Append a blank line, then the following block to the end of the file:

```groovy

import nextflow.gradle.GenerateLineageSchemaTask

tasks.register('generateLineageSchema', GenerateLineageSchemaTask) {
    description = 'Generate JSON Schema for the lineage model v1beta1'
    group = 'documentation'
    dependsOn compileGroovy
    classpath = sourceSets.main.runtimeClasspath
    // Keep this list in sync with LinTypeAdapterFactory.registerSubtype(...) calls
    // in src/main/nextflow/lineage/serde/LinTypeAdapterFactory.groovy.
    // After editing, re-run this task to refresh
    // src/resources/schema/lineage-v1beta1.schema.json.
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

Note: in Groovy Gradle scripts, top-level `import` statements must appear *before* the `plugins { ... }` and any other code. If Gradle complains about the import position, move the `import nextflow.gradle.GenerateLineageSchemaTask` line to the very top of the file (above `apply plugin: 'groovy'`) and leave the `tasks.register(...)` block at the end.

- [ ] **Step 2: Verify the task is discoverable**

Run:
```
./gradlew :nf-lineage:tasks --group documentation
```

Expected: output lists `generateLineageSchema` under the `Documentation tasks` group with the description "Generate JSON Schema for the lineage model v1beta1".

- [ ] **Step 3: Commit**

```
git add modules/nf-lineage/build.gradle
git commit -s -m "Register generateLineageSchema task in nf-lineage build"
```

---

## Task 5: Run the task and verify the generated schema

This is the manual verification gate mandated by the spec. No code is written; the goal is to prove the pipeline end-to-end produces a correct schema.

**Files:**
- Create (via task run): `modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json`

- [ ] **Step 1: Run the task**

```
./gradlew :nf-lineage:generateLineageSchema
```

Expected: `BUILD SUCCESSFUL`, lifecycle log line `Wrote lineage schema with 6 subtypes to .../lineage-v1beta1.schema.json`. If the build fails, debug and fix in `GenerateLineageSchemaTask.groovy` before continuing.

- [ ] **Step 2: Confirm the file exists and is valid JSON**

```
ls -la modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json
python3 -c "import json; json.load(open('modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json'))"
```

Expected: file exists and python's `json.load` runs without exception.

- [ ] **Step 3: Confirm structural correctness**

```
python3 -c "
import json
s = json.load(open('modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json'))
assert s.get('\$schema') == 'https://json-schema.org/draft/2020-12/schema', s.get('\$schema')
assert s.get('title') == 'Lineage v1beta1'
oneof = s.get('oneOf')
assert isinstance(oneof, list) and len(oneof) == 6, len(oneof) if isinstance(oneof, list) else type(oneof)
kinds = [b['properties']['kind']['const'] for b in oneof]
assert set(kinds) == {'WorkflowRun', 'WorkflowOutput', 'Workflow', 'TaskRun', 'TaskOutput', 'FileOutput'}, kinds
versions = {b['properties']['version']['const'] for b in oneof}
assert versions == {'lineage/v1beta1'}, versions
print('schema structure OK')
"
```

Expected: prints `schema structure OK` and exits 0.

- [ ] **Step 4: Validate a known-good lineage JSON document against the schema**

Locate an existing lineage encoder fixture or unit-test sample. Search for one:

```
grep -rln "lineage/v1beta1" modules/nf-lineage/src/test | head -3
```

Pick one of the listed files, find a JSON literal inside it, save it to `/tmp/sample-lineage.json`. If no fixture is convenient, hand-craft a minimal `FileOutput` document:

```json
{
  "version": "lineage/v1beta1",
  "kind": "FileOutput",
  "spec": {
    "path": "/tmp/out.txt",
    "size": 42
  }
}
```

Save to `/tmp/sample-lineage.json`, then validate:

```
pip install --user jsonschema 2>/dev/null
python3 -c "
import json
from jsonschema import Draft202012Validator
schema = json.load(open('modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json'))
doc = json.load(open('/tmp/sample-lineage.json'))
Draft202012Validator(schema).validate(doc)
print('document valid')
"
```

Expected: prints `document valid`. If validation fails, inspect the schema, identify the mismatch (likely an over-strict `required` or an unsupported type mapping), adjust `GenerateLineageSchemaTask.generateSubtypeSchema`, regenerate, and re-validate.

- [ ] **Step 5: Commit the generated schema**

```
git add modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json
git commit -s -m "Add generated lineage v1beta1 JSON Schema"
```

---

## Task 6: Add maintainer sync comment to `LinTypeAdapterFactory`

The only production-code touch. No behavior change.

**Files:**
- Modify: `modules/nf-lineage/src/main/nextflow/lineage/serde/LinTypeAdapterFactory.groovy`

- [ ] **Step 1: Insert the comment**

Open the file. Find the constructor (currently lines 47-55):

```groovy
    LinTypeAdapterFactory() {
        super(LinSerializable.class, "kind", false)
        this.registerSubtype(WorkflowRun, WorkflowRun.simpleName)
            .registerSubtype(WorkflowOutput, WorkflowOutput.simpleName)
            .registerSubtype(Workflow, Workflow.simpleName)
            .registerSubtype(TaskRun, TaskRun.simpleName)
            .registerSubtype(TaskOutput, TaskOutput.simpleName)
            .registerSubtype(FileOutput, FileOutput.simpleName)
    }
```

Replace it with:

```groovy
    LinTypeAdapterFactory() {
        super(LinSerializable.class, "kind", false)
        // When adding or removing a subtype here, also update the `subtypes` list in
        // modules/nf-lineage/build.gradle (task `generateLineageSchema`) and re-run
        // `./gradlew :nf-lineage:generateLineageSchema` to refresh
        // src/resources/schema/lineage-v1beta1.schema.json.
        this.registerSubtype(WorkflowRun, WorkflowRun.simpleName)
            .registerSubtype(WorkflowOutput, WorkflowOutput.simpleName)
            .registerSubtype(Workflow, Workflow.simpleName)
            .registerSubtype(TaskRun, TaskRun.simpleName)
            .registerSubtype(TaskOutput, TaskOutput.simpleName)
            .registerSubtype(FileOutput, FileOutput.simpleName)
    }
```

- [ ] **Step 2: Confirm no behavior change — run the module's existing tests**

```
./gradlew :nf-lineage:test
```

Expected: all tests pass. The change is comment-only; this is a sanity check.

- [ ] **Step 3: Re-run the schema task and confirm the output is unchanged**

```
./gradlew :nf-lineage:generateLineageSchema
git diff modules/nf-lineage/src/resources/schema/lineage-v1beta1.schema.json
```

Expected: `git diff` shows no changes. (The comment doesn't affect generation, so the schema must be byte-identical.)

- [ ] **Step 4: Commit**

```
git add modules/nf-lineage/src/main/nextflow/lineage/serde/LinTypeAdapterFactory.groovy
git commit -s -m "Add sync note pointing maintainers to generateLineageSchema task"
```

---

## Self-Review

Spec coverage:
- Goal 1 (single JSON Schema document, six `oneOf` branches): Tasks 2–5.
- Goal 2 (reproducible Gradle task): Task 4.
- Goal 3 (checked-in artifact): Task 5 commit.
- Goal 4 (no runtime deps on nf-lineage): Task 1 (dep lives in `buildSrc/build.gradle` only); verified implicitly by Task 6, step 2 running nf-lineage tests without any new classpath entries.
- Non-Goals respected: no CI enforcement task, no runtime validation, no auto-discovery.
- US1 acceptance (six branches, correct consts, valid JSON): Task 5, step 3.
- US1 acceptance (third-party validator passes): Task 5, step 4.
- US2 (refresh after model change): manual workflow documented via the comment added in Task 6; no separate task needed since the regeneration path is identical to Task 5.
- US3 (new subtype): the sync comment in Task 6 + the build.gradle comment in Task 4 cover the maintainer workflow.

Placeholder scan: none — all steps have concrete commands, code blocks, or assertions.

Type consistency: `classpath` (FileCollection), `subtypes` (List<String>), `outputFile` (File) are declared in Task 2 and referenced consistently in Tasks 4 and 5. `generateSubtypeSchema(Class<?>, ObjectMapper) → ObjectNode` defined in Task 2, replaced in Task 3 with the same signature.

No issues found.
