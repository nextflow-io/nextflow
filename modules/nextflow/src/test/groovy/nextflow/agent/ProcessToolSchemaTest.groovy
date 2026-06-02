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
package nextflow.agent

import java.nio.file.Files

import nextflow.Session
import nextflow.script.ProcessDef
import nextflow.script.ScriptMeta
import nextflow.script.parser.v2.ScriptLoaderV2
import test.Dsl2Spec

/**
 * Verifies that {@link ProcessToolSchema} derives portable JSON-schema maps from
 * a typed process's declared inputs/outputs, and fails loudly for kinds not yet
 * supported as agent tools (tuple, path, ...).
 *
 * Schemas are obtained by loading a real Nextflow script so the actual compiled
 * {@code ProcessConfigV2} / typed-I/O metadata is exercised.
 */
class ProcessToolSchemaTest extends Dsl2Spec {

    private ProcessDef loadProcess(String processDecl) {
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = """
            nextflow.enable.types = true

            ${processDecl}

            workflow {
            }
            """.stripIndent()
        parser.parse(file)
        parser.runScript()
        return ScriptMeta.get(parser.script).getProcess('greet')
    }

    def 'should derive the input schema from a scalar typed input'() {
        given:
        def proc = loadProcess('''
            process greet {
                input:
                name: String

                output:
                result: String

                exec:
                result = "Hello ${name}!"
            }
            ''')

        expect:
        ProcessToolSchema.inputSchema(proc) == [
            type: 'object',
            properties: [name: [type: 'string']],
            required: ['name'],
            additionalProperties: false,
        ]
    }

    def 'should mark an optional input as not required and map number/integer/boolean'() {
        given:
        def proc = loadProcess('''
            process greet {
                input:
                name: String
                count: Integer?
                score: Double
                flag: Boolean

                output:
                result: String

                exec:
                result = "x"
            }
            ''')

        when:
        def schema = ProcessToolSchema.inputSchema(proc)

        then:
        schema.properties == [
            name : [type: 'string'],
            count: [type: 'integer'],
            score: [type: 'number'],
            flag : [type: 'boolean'],
        ]
        schema.required == ['name', 'score', 'flag']
        schema.additionalProperties == false
    }

    def 'should derive the output schema from a named scalar output'() {
        given:
        def proc = loadProcess('''
            process greet {
                input:
                name: String

                output:
                answer: String

                exec:
                answer = "Hello ${name}!"
            }
            ''')

        expect:
        ProcessToolSchema.outputSchema(proc) == [
            type: 'object',
            properties: [answer: [type: 'string']],
            required: ['answer'],
            additionalProperties: false,
        ]
    }

    def 'should throw a loud error for a tuple input'() {
        given:
        def proc = loadProcess('''
            process greet {
                input:
                tuple(id: String, num: Integer)

                output:
                result: String

                exec:
                result = "x"
            }
            ''')

        when:
        ProcessToolSchema.inputSchema(proc)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('greet')
        e.message.contains('not yet supported as an agent tool')
    }

    def 'should throw a loud error for a path input'() {
        given:
        def proc = loadProcess('''
            process greet {
                input:
                sample: Path

                output:
                result: String

                exec:
                result = "x"
            }
            ''')

        when:
        ProcessToolSchema.inputSchema(proc)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('sample')
        e.message.contains('not yet supported as an agent tool')
    }
}
