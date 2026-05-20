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

import java.nio.file.Path
import java.time.OffsetDateTime

import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.ArrayNode
import com.fasterxml.jackson.databind.node.ObjectNode
import com.github.victools.jsonschema.generator.CustomDefinition
import com.github.victools.jsonschema.generator.Option
import com.github.victools.jsonschema.generator.OptionPreset
import com.github.victools.jsonschema.generator.SchemaBuilder
import com.github.victools.jsonschema.generator.SchemaGenerator
import com.github.victools.jsonschema.generator.SchemaGeneratorConfigBuilder
import com.github.victools.jsonschema.generator.SchemaVersion
import groovy.transform.CompileStatic
import org.gradle.api.DefaultTask
import org.gradle.api.file.FileCollection
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.InputDirectory
import org.gradle.api.tasks.InputFiles
import org.gradle.api.tasks.OutputFile
import org.gradle.api.tasks.TaskAction

/**
 * Generates a JSON Schema (draft 2020-12) describing the JSON documents emitted
 * by {@code nextflow.lineage.serde.LinEncoder} for v1beta1 model classes.
 *
 * The task does not depend on the nf-lineage runtime; it loads compiled model
 * classes via a URLClassLoader and uses victools jsonschema-generator to derive
 * the per-subtype schemas. A single SchemaBuilder collects $defs across all
 * subtypes so shared types (Checksum, DataPath, Parameter) are emitted once and
 * referenced via $ref. Each subtype is wrapped in a {version, kind, spec}
 * envelope matching {@code LinTypeAdapterFactory.write(...)}.
 */
@CompileStatic
class GenerateLineageSchemaTask extends DefaultTask {

    @InputFiles
    FileCollection classpath

    @Input
    List<String> subtypes

    @InputDirectory
    File modelSourceDir

    @OutputFile
    File outputFile

    private static final String MODEL_PACKAGE = 'nextflow.lineage.model.v1beta1.'

    @TaskAction
    void generate() {
        final loader = buildClassLoader()
        final version = readLineageVersion(loader)
        final descriptions = parseClassDocs(modelSourceDir)
        final mapper = new ObjectMapper()
        final root = mapper.createObjectNode()
        root.put('$schema', 'https://json-schema.org/draft/2020-12/schema')
        root.put('title', 'Nextflow Lineage Model v1beta1')

        final schemaBuilder = newSchemaBuilder(mapper, descriptions)
        final oneOf = mapper.createArrayNode()

        subtypes.each { String fqn ->
            final cls = loader.loadClass(fqn)
            final specRef = schemaBuilder.createSchemaReference(cls) as ObjectNode
            oneOf.add(wrapEnvelope(cls.simpleName, version, specRef, mapper))
        }

        // Collect consolidated $defs from the builder; emit at root so $ref
        // targets like "#/$defs/Checksum" resolve correctly.
        final defs = schemaBuilder.collectDefinitions('$defs')
        root.set('$defs', defs)
        root.set('oneOf', oneOf)

        // Victools, when combining NULLABLE_FIELDS_BY_DEFAULT with
        // DEFINITIONS_FOR_ALL_OBJECTS, emits "<Type>-nullable" wrapper defs.
        // Inline those wrappers at the call site and drop the def entries.
        inlineNullableWrappers(root)

        outputFile.parentFile.mkdirs()
        mapper.writerWithDefaultPrettyPrinter().writeValue(outputFile, root)
        logger.lifecycle("Wrote lineage schema with ${subtypes.size()} subtypes to ${outputFile}")
    }

    private URLClassLoader buildClassLoader() {
        final urls = classpath.files.collect { it.toURI().toURL() } as URL[]
        return new URLClassLoader(urls, getClass().classLoader)
    }

    private static String readLineageVersion(ClassLoader loader) {
        final linModel = loader.loadClass('nextflow.lineage.model.v1beta1.LinModel')
        return linModel.getField('VERSION').get(null) as String
    }

    private static SchemaBuilder newSchemaBuilder(ObjectMapper mapper, Map<String, String> descriptions) {
        // LinEncoder serializes nulls (withSerializeNulls(true)) so every field
        // can appear as null in the emitted JSON — schema marks them nullable.
        final config = new SchemaGeneratorConfigBuilder(SchemaVersion.DRAFT_2020_12, OptionPreset.PLAIN_JSON)
            .with(Option.NULLABLE_FIELDS_BY_DEFAULT)
            .with(Option.NULLABLE_ARRAY_ITEMS_ALLOWED)
            .with(Option.DEFINITIONS_FOR_ALL_OBJECTS)
        config.forTypesInGeneral().withCustomDefinitionProvider({ javaType, context ->
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
        })
        config.forTypesInGeneral().withTitleResolver({ scope ->
            final cls = scope.type.erasedType
            cls.name.startsWith(MODEL_PACKAGE) ? cls.simpleName : null
        })
        config.forTypesInGeneral().withDescriptionResolver({ scope ->
            descriptions[scope.type.erasedType.name]
        })
        return new SchemaGenerator(config.build()).buildMultipleSchemaDefinitions()
    }

    private static Map<String, String> parseClassDocs(File dir) {
        final result = [:] as Map<String, String>
        if (dir == null || !dir.isDirectory()) return result
        final pattern = ~/(?s)\/\*\*(.*?)\*\/\s*(?:@\w+(?:\([^)]*\))?\s*)*(?:abstract\s+|final\s+|public\s+|static\s+)*(?:class|interface)\s+(\w+)/
        dir.eachFileRecurse { File file ->
            if (!file.name.endsWith('.groovy')) return
            final text = file.text
            final pkgMatcher = text =~ /(?m)^\s*package\s+([\w.]+)/
            if (!pkgMatcher.find()) return
            final pkg = pkgMatcher.group(1) as String
            final m = text =~ pattern
            while (m.find()) {
                final raw = m.group(1) as String
                final clsName = m.group(2) as String
                final cleaned = cleanDoc(raw)
                if (cleaned) result["${pkg}.${clsName}".toString()] = cleaned
            }
        }
        return result
    }

    private static String cleanDoc(String raw) {
        // Strip the `*` line prefix and drop everything from the first @tag onwards.
        final lines = raw.readLines().collect { line -> line.replaceFirst(/^\s*\*\s?/, '') }
        final stopIdx = lines.findIndexOf { (it =~ /^\s*@\w+/).find() }
        final kept = stopIdx >= 0 ? lines[0..<stopIdx] : lines
        final out = kept.join('\n').trim()
        return out ?: null
    }

    private static void inlineNullableWrappers(ObjectNode root) {
        final defs = root.get('$defs') as ObjectNode
        if (defs == null) return
        final wrappers = [:] as Map<String, ObjectNode>
        defs.fieldNames().each { String name ->
            if (name.endsWith('-nullable'))
                wrappers[name] = defs.get(name) as ObjectNode
        }
        if (wrappers.isEmpty()) return
        replaceNullableRefs(root, wrappers)
        wrappers.keySet().each { defs.remove(it) }
    }

    private static void replaceNullableRefs(JsonNode node, Map<String, ObjectNode> wrappers) {
        if (node instanceof ObjectNode) {
            final on = node as ObjectNode
            final ref = on.get('$ref')
            if (ref != null && on.size() == 1) {
                final target = ref.asText() - '#/$defs/'
                final body = wrappers[target]
                if (body != null) {
                    on.removeAll()
                    on.setAll(body.deepCopy() as ObjectNode)
                    return
                }
            }
            on.fields().each { replaceNullableRefs(it.value, wrappers) }
        }
        else if (node instanceof ArrayNode) {
            (node as ArrayNode).each { replaceNullableRefs(it as JsonNode, wrappers) }
        }
    }

    private static ObjectNode wrapEnvelope(String kind, String version, JsonNode spec, ObjectMapper mapper) {
        final env = mapper.createObjectNode()
        env.put('type', 'object')
        final props = env.putObject('properties')
        props.putObject('version').put('const', version)
        props.putObject('kind').put('const', kind)
        props.set('spec', spec)
        final required = env.putArray('required')
        required.add('version')
        required.add('kind')
        required.add('spec')
        env.put('additionalProperties', false)
        return env
    }
}
