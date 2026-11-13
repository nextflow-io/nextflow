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

package nextflow.module

import java.nio.file.Files
import java.nio.file.Path

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for ModuleSchemaValidator.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleSchemaValidatorTest extends Specification {

    @TempDir
    Path tempDir

    private static final String MINIMAL_SCHEMA = '''\
        {
          "$schema": "https://json-schema.org/draft/2020-12/schema",
          "type": "object",
          "properties": {
            "name": { "type": "string" },
            "description": { "type": "string" },
            "input": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "type": { "type": "string", "enum": ["string", "file", "directory"] },
                  "description": { "type": "string" }
                },
                "required": ["type", "description"]
              }
            }
          },
          "required": ["name", "description"]
        }
        '''.stripIndent()

    private Path writeSchema(String text = MINIMAL_SCHEMA) {
        final p = tempDir.resolve('schema.json')
        Files.writeString(p, text)
        return p
    }

    private Path writeMeta(String yaml) {
        final p = tempDir.resolve('meta.yml')
        Files.writeString(p, yaml)
        return p
    }

    def 'should pass validation when meta.yml satisfies the schema' () {
        given:
        def schema = writeSchema()
        def meta = writeMeta('''\
            name: nf-core/fastqc
            description: Run FastQC
            '''.stripIndent())

        when:
        def errors = ModuleSchemaValidator.validate(meta, schema.toString())

        then:
        errors.isEmpty()
    }

    def 'should report missing required fields against the schema' () {
        given:
        def schema = writeSchema()
        def meta = writeMeta('''\
            name: nf-core/fastqc
            '''.stripIndent())

        when:
        def errors = ModuleSchemaValidator.validate(meta, schema.toString())

        then:
        !errors.isEmpty()
        errors.any { it.contains('description') }
    }

    def 'should report invalid input type against the schema enum' () {
        given:
        def schema = writeSchema()
        def meta = writeMeta('''\
            name: nf-core/fastqc
            description: Run FastQC
            input:
              - name: reads
                type: bogus
                description: input reads
            '''.stripIndent())

        when:
        def errors = ModuleSchemaValidator.validate(meta, schema.toString())

        then:
        !errors.isEmpty()
        errors.any { it.toLowerCase().contains('type') || it.toLowerCase().contains('enum') }
    }

    def 'should accept a file: URI for the schema location' () {
        given:
        def schema = writeSchema()
        def meta = writeMeta('''\
            name: nf-core/fastqc
            description: Run FastQC
            '''.stripIndent())

        when:
        def errors = ModuleSchemaValidator.validate(meta, schema.toUri().toString())

        then:
        errors.isEmpty()
    }

    def 'should hard-fail when the schema cannot be loaded' () {
        given:
        def meta = writeMeta('name: x\ndescription: y\n')

        when:
        ModuleSchemaValidator.validate(meta, tempDir.resolve('does-not-exist.json').toString())

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Failed to load module schema')
    }

    def 'should hard-fail when the schema does not declare a supported draft' () {
        given:
        def schema = writeSchema('''\
            {
              "type": "object",
              "properties": {
                "name": { "type": "string" }
              }
            }
            '''.stripIndent())
        def meta = writeMeta('name: x\n')

        when:
        ModuleSchemaValidator.validate(meta, schema.toString())

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Cannot determine JSON Schema draft')
    }

    def 'resolveSchemaLocation: explicit value wins over env, meta and default' () {
        given:
        def meta = writeMeta('name: x\n$schema: /from/meta.json\n')
        SysEnv.push([(ModuleSchemaValidator.SCHEMA_ENV_VAR): '/from/env.json'])

        when:
        def loc = ModuleSchemaValidator.resolveSchemaLocation(meta, '/explicit.json')

        then:
        loc == '/explicit.json'

        cleanup:
        SysEnv.pop()
    }

    def 'resolveSchemaLocation: env var wins over meta $schema and default' () {
        given:
        def meta = writeMeta('name: x\n$schema: /from/meta.json\n')
        SysEnv.push([(ModuleSchemaValidator.SCHEMA_ENV_VAR): '/from/env.json'])

        when:
        def loc = ModuleSchemaValidator.resolveSchemaLocation(meta, null)

        then:
        loc == '/from/env.json'

        cleanup:
        SysEnv.pop()
    }

    def 'resolveSchemaLocation: meta $schema used when no explicit or env override' () {
        given:
        def meta = writeMeta('name: x\n$schema: /from/meta.json\n')

        expect:
        ModuleSchemaValidator.resolveSchemaLocation(meta, null) == '/from/meta.json'
    }

    def 'resolveSchemaLocation: falls back to the default schema url' () {
        given:
        def meta = writeMeta('name: x\n')

        expect:
        ModuleSchemaValidator.resolveSchemaLocation(meta, null) == ModuleSchemaValidator.DEFAULT_SCHEMA_URL
    }
}
