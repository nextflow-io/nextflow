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
package nextflow.script.parser.v2

import java.nio.file.Files

import nextflow.Session
import nextflow.script.AgentDef
import nextflow.script.ScriptMeta
import test.Dsl2Spec

/**
 * End-to-end smoke test: parsing a script with an `agent` block must
 * register an {@link AgentDef} on the script's {@link ScriptMeta}.
 */
class AgentScriptLoadingTest extends Dsl2Spec {

    def 'should load a script with an agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            agent hello_agent {
                prompt:
                "hello"
            }

            workflow {
            }
            '''.stripIndent()

        when:
        parser.parse(file)
        parser.runScript()

        then:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        definitions.any { it instanceof AgentDef && it.name == 'hello_agent' }

        cleanup:
        file.parent.deleteDir()
    }

    def 'should load a script with a directive-rich agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    question: String

                output:
                    plan: String

                prompt:
                """
                Question: ${question}
                """
            }

            workflow {
            }
            '''.stripIndent()

        when:
        parser.parse(file)
        parser.runScript()

        then:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        def agent = definitions.find { it instanceof AgentDef && it.name == 'eval_agent' } as AgentDef
        agent != null
        agent.model == 'openai/gpt-5-mini'
        agent.instruction == 'You are helpful.'
        agent.maxIterations == 20
        agent.tools == []
        agent.inputs*.name == ['question']
        agent.outputs*.name == ['plan']
        agent.prompt != null
        agent.prompt.source.contains('Question:')

        cleanup:
        file.parent.deleteDir()
    }
}
