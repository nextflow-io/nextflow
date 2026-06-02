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
import nextflow.script.AgentDef
import nextflow.script.ScriptMeta
import nextflow.script.parser.v2.ScriptLoaderV2
import test.Dsl2Spec

/**
 * Verifies that {@link RecordSchema#of} reflects a record output class into a
 * portable JSON-schema map, rejecting unsupported field types.
 *
 * The record classes are obtained by loading a real Nextflow script (so the
 * actual compiled record classes + {@code @Nullable} annotations are exercised),
 * mirroring how an agent output type is resolved at runtime.
 */
class RecordSchemaTest extends Dsl2Spec {

    private Class loadOutputType(String recordDecls, String outputDecl) {
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = """
            nextflow.enable.types = true

            ${recordDecls}

            agent probe_agent {
                input:
                    q: Question

                output:
                    ${outputDecl}

                prompt:
                "hello"
            }

            workflow {
            }
            """.stripIndent()
        parser.parse(file)
        parser.runScript()
        def agent = ScriptMeta.get(parser.script).getDefinitions().find { it instanceof AgentDef } as AgentDef
        return agent.outputs[0].type
    }

    def 'should derive a JSON schema from a record output type'() {
        given:
        def cls = loadOutputType('''
            record Question { text: String }
            record Answer { answer: String; confidence: Double }
            ''', 'a: Answer')

        when:
        def schema = RecordSchema.of(cls)

        then:
        schema == [
            type: 'object',
            properties: [
                answer    : [type: 'string'],
                confidence: [type: 'number'],
            ],
            required: ['answer', 'confidence'],
            additionalProperties: false,
        ]
    }

    def 'should mark optional fields as not required'() {
        given:
        def cls = loadOutputType('''
            record Question { text: String }
            record Answer { answer: String; note: String? }
            ''', 'a: Answer')

        when:
        def schema = RecordSchema.of(cls)

        then:
        schema.properties.keySet() == ['answer', 'note'] as Set
        schema.required == ['answer']
    }

    def 'should reject a Path output field'() {
        given:
        def cls = loadOutputType('''
            record Question { text: String }
            record WithPath { p: Path }
            ''', 'a: WithPath')

        when:
        RecordSchema.of(cls)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('p')
    }
}
