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
package nextflow.script.parser

import nextflow.script.ast.AgentNode
import nextflow.script.ast.ScriptNode
import nextflow.script.control.ScriptParser
import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * @see nextflow.script.parser.ScriptAstBuilder
 */
class AgentParserTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    List<SyntaxException> check(String contents) {
        return TestUtils.check(scriptParser, contents)
    }

    ScriptNode parse(String contents) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents.stripIndent())
        scriptParser.analyze()
        assert !TestUtils.hasSyntaxErrors(source)
        return source.getAST() as ScriptNode
    }

    def 'should parse a minimal agent definition'() {
        when:
        def script = parse('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    val question

                output:
                    val plan

                prompt:
                """
                Question: ${question}
                """
            }
            ''')

        then:
        script.agents.size() == 1
        def node = script.agents[0] as AgentNode
        node.name == 'eval_agent'
        node.inputs.length == 1
        node.inputs[0].name == 'question'
    }
}
