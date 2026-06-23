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

    def 'should load a script with a minimal record-typed agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            record Question { text: String }
            record Answer { answer: String }

            agent hello_agent {
                input:
                    q: Question

                output:
                    a: Answer

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

    def 'should load a script with a val-typed agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'

                input:
                    question: String

                output:
                    answer: String

                prompt:
                """
                Answer: ${question}
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
        def agent = definitions.find { it instanceof AgentDef && it.name == 'qa' } as AgentDef
        agent != null
        agent.inputs*.name == ['question']
        agent.outputs*.name == ['answer']
        and:
        // val I/O resolve to scalar (String) types, not record classes
        (agent.inputs[0].type as Class) == String
        (agent.outputs[0].type as Class) == String

        cleanup:
        file.parent.deleteDir()
    }

    def 'should load a script with a directive-rich record-typed agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            record Question { text: String; context: String? }
            record Answer { answer: String; confidence: Double }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    q: Question

                output:
                    a: Answer

                prompt:
                """
                Question: ${q.text}
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
        agent.inputs*.name == ['q']
        agent.outputs*.name == ['a']
        agent.prompt != null
        agent.prompt.source.contains('Question:')
        and:
        // I/O are named record types: the resolved output type is the compiled
        // record class (e.g. `Answer`) with the declared fields
        def inputType = agent.inputs[0].type as Class
        def outputType = agent.outputs[0].type as Class
        inputType.name.endsWith('Question')
        outputType.name.endsWith('Answer')
        (outputType.declaredFields*.name as Set).containsAll(['answer', 'confidence'])

        cleanup:
        file.parent.deleteDir()
    }

    def 'should load a script with capability tools (module_run + filesystem) and expose them via getTools()'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            process word_stats {
                input:  text: String
                output: stats: String
                exec:   stats = "{}"
            }

            agent analyst {
                model 'openai/gpt-5-mini'
                instruction 'Analyse text.'

                tools 'module_run', 'filesystem'

                input:
                    text: String
                output:
                    summary: String

                prompt:
                """
                Analyse: ${text}
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
        def agent = definitions.find { it instanceof AgentDef && it.name == 'analyst' } as AgentDef
        agent != null
        agent.tools.contains('module_run')
        agent.tools.contains('filesystem')
        agent.tools.size() == 2

        cleanup:
        file.parent.deleteDir()
    }

    def 'should load a script with a goal directive and expose it via getGoal()'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-5-mini'
                instruction 'be concise'
                goal 'answer the question accurately'

                input:
                    question: String

                output:
                    answer: String

                prompt:
                """
                ${question}
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
        def agent = definitions.find { it instanceof AgentDef && it.name == 'qa' } as AgentDef
        agent != null
        agent.goal == 'answer the question accurately'
        agent.instruction == 'be concise'

        cleanup:
        file.parent.deleteDir()
    }
}
