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

package nextflow.config

import java.nio.file.Files
import java.nio.file.Path

import nextflow.SysEnv
import spock.lang.Specification

/**
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class SchemaParamsHelperTest extends Specification {

    Path tempDir

    def setup() {
        tempDir = Files.createTempDirectory('nf-schema-test')
        SysEnv.push([:])
    }

    def cleanup() {
        SysEnv.pop()
        tempDir.toFile().deleteDir()
    }

    private Path writeSchema(String content) {
        final f = tempDir.resolve('nextflow_schema.json')
        f.text = content
        return f
    }

    def 'reads top-level types from a flat schema'() {
        given:
        writeSchema('''
        {
            "properties": {
                "input":  { "type": "string"  },
                "cpus":   { "type": "integer" },
                "ratio":  { "type": "number"  },
                "skip":   { "type": "boolean" }
            }
        }
        '''.stripIndent())

        when:
        def types = SchemaParamsHelper.readSchemaTypes(tempDir.resolve('nextflow_schema.json'))

        then:
        types == [input: 'string', cpus: 'integer', ratio: 'number', skip: 'boolean']
    }

    def 'reads types nested under definitions / $defs / allOf'() {
        given:
        writeSchema('''
        {
            "allOf": [
                { "$ref": "#/$defs/group_one" }
            ],
            "definitions": {
                "group_one": {
                    "properties": {
                        "input":   { "type": "string"  },
                        "max_cpus":{ "type": "integer" }
                    }
                }
            },
            "$defs": {
                "group_two": {
                    "properties": {
                        "skip_qc": { "type": "boolean" }
                    }
                }
            }
        }
        '''.stripIndent())

        when:
        def types = SchemaParamsHelper.readSchemaTypes(tempDir.resolve('nextflow_schema.json'))

        then:
        types == [input: 'string', max_cpus: 'integer', skip_qc: 'boolean']
    }

    def 'coerces CLI string params to schema-declared types'() {
        given:
        writeSchema('''
        {
            "properties": {
                "input":   { "type": "string"  },
                "cpus":    { "type": "integer" },
                "ratio":   { "type": "number"  },
                "skip":    { "type": "boolean" }
            }
        }
        '''.stripIndent())
        def cli = [input: 'data.csv', cpus: '8', ratio: '0.5', skip: 'true', extra: 'hello']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.input == 'data.csv'
        cli.cpus == 8
        cli.cpus instanceof Integer
        cli.ratio == 0.5f
        cli.ratio instanceof Float
        cli.skip == Boolean.TRUE
        cli.extra == 'hello' // not in schema -- left untouched
    }

    def 'leaves un-coercible values as strings'() {
        given:
        writeSchema('''
        {
            "properties": {
                "cpus": { "type": "integer" },
                "skip": { "type": "boolean" }
            }
        }
        '''.stripIndent())
        def cli = [cpus: 'abc', skip: 'maybe']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.cpus == 'abc'
        cli.skip == 'maybe'
    }

    def 'does nothing when schema file is missing'() {
        given:
        def cli = [cpus: '8']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.cpus == '8'
    }

    def 'does nothing when schema is malformed'() {
        given:
        writeSchema('{ broken')
        def cli = [cpus: '8']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.cpus == '8'
    }

    def 'is a no-op when params type detection is disabled'() {
        given:
        writeSchema('{ "properties": { "cpus": { "type": "integer" } } }')
        SysEnv.push(NXF_DISABLE_PARAMS_TYPE_DETECTION: 'true')
        def cli = [cpus: '8']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.cpus == '8'

        cleanup:
        SysEnv.pop()
    }

    def 'is a no-op under syntax parser v1'() {
        given:
        writeSchema('{ "properties": { "cpus": { "type": "integer" } } }')
        SysEnv.push(NXF_SYNTAX_PARSER: 'v1')
        def cli = [cpus: '8']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.cpus == '8'

        cleanup:
        SysEnv.pop()
    }

    def 'ignores schema properties that have no type field'() {
        given:
        writeSchema('''
        {
            "properties": {
                "input":  { "description": "no type given here" },
                "cpus":   { "type": "integer" }
            }
        }
        '''.stripIndent())
        def cli = [input: '42', cpus: '8']

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.input == '42'   // schema didn't declare a type -- left as string
        cli.cpus == 8
    }

    def 'first declaration wins when a property is repeated across definitions'() {
        given:
        writeSchema('''
        {
            "definitions": {
                "first": {
                    "properties": { "size": { "type": "integer" } }
                },
                "second": {
                    "properties": { "size": { "type": "string" } }
                }
            }
        }
        '''.stripIndent())

        when:
        def types = SchemaParamsHelper.readSchemaTypes(tempDir.resolve('nextflow_schema.json'))

        then:
        types == [size: 'integer']
    }

    def 'coerces booleans regardless of letter casing'() {
        given:
        writeSchema('{ "properties": { "flag": { "type": "boolean" } } }')

        expect:
        coerceFlag(input) == expected

        where:
        input   | expected
        'TRUE'  | Boolean.TRUE
        'False' | Boolean.FALSE
        'true'  | Boolean.TRUE
        'false' | Boolean.FALSE
    }

    private Object coerceFlag(String input) {
        def cli = [flag: input]
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)
        return cli.flag
    }

    def 'preserves non-string values already typed by upstream parsing'() {
        given:
        writeSchema('{ "properties": { "cpus": { "type": "integer" } } }')
        // value supplied as Integer (e.g. from a JSON params file) is left alone
        def cli = [cpus: 16] as Map<String,Object>

        when:
        SchemaParamsHelper.applySchemaTypes(tempDir, cli)

        then:
        cli.cpus == 16
        cli.cpus instanceof Integer
    }
}
