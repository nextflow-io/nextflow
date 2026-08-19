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
import nextflow.script.control.ScriptToGroovyHelper
import org.codehaus.groovy.control.SourceUnit
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

    SourceUnit parseSource(String contents) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents.stripIndent())
        scriptParser.analyze()
        assert !TestUtils.hasSyntaxErrors(source)
        return source
    }

    // -- T3 (design §7.2/D3): the prompt closure's free-variable refs (params.*, task.ext.*)
    //    are captured via the same collector that populates process-body BodyDef.valRefs, so
    //    AgentToGroovyVisitor can fold them into the synthetic PromptDef/BodyDef cache key.

    def 'should capture params.* prompt globals as prompt valRefs (excluding declared inputs)'() {
        given:
        def source = parseSource('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                tools()

                input:
                question: String

                output:
                answer: String

                prompt:
                """
                ${question} threshold=${params.threshold} args=${task.ext.args}
                """
            }
            ''')

        when:
        def node = (source.getAST() as ScriptNode).agents[0] as AgentNode
        def refs = new ScriptToGroovyHelper(source).getVariableRefs(node.prompt)
        def names = refs.expressions
            .collect { expr -> expr.arguments.expressions[0].value }
            .sort()

        then: 'params.* and task.ext.* are captured; the declared input `question` is NOT'
        names == ['params.threshold', 'task.ext.args']
    }

    def 'should capture no prompt valRefs when the prompt references only declared inputs'() {
        given:
        def source = parseSource('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                tools()

                input:
                question: String

                output:
                answer: String

                prompt:
                """
                Question: ${question}
                """
            }
            ''')

        when:
        def node = (source.getAST() as ScriptNode).agents[0] as AgentNode
        def refs = new ScriptToGroovyHelper(source).getVariableRefs(node.prompt)

        then:
        refs.expressions.isEmpty()
    }

    def 'should parse a minimal agent definition'() {
        when:
        def script = parse('''\
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                plan: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                question: Question

                output:
                plan: Answer

                prompt:
                """
                Question: ${question.text}
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

    def 'should report an error for agent without prompt section'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            record Question {
                text: String
            }

            agent broken {
                model 'openai/gpt-5-mini'
                instruction 'x'
                tools()

                input:
                q: Question
            }
            ''')

        then:
        errors.size() == 1
        errors[0].getOriginalMessage() == 'Invalid agent definition -- check for missing or out-of-order section labels'
    }

    def 'should resolve an agent reference from a workflow'() {
        when:
        def script = parse('''\
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                text: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'x'
                tools()

                input:
                q: Question

                output:
                r: Answer

                prompt:
                """
                ${q.text}
                """
            }

            workflow {
                channel.of('hi') | eval_agent | view
            }
            ''')

        then:
        script.agents.size() == 1
        script.workflows.size() == 1
        // The `parse` helper asserts no syntax errors. If `eval_agent` failed
        // to resolve in the workflow body, that assertion would fail.
    }

    def 'should resolve directives and prompt variables in an agent body'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                plan: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                question: Question

                output:
                plan: Answer

                prompt:
                """
                Question: ${question.text}
                """
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should resolve the repeatable label directive in an agent body'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent eval_agent {
                label 'reasoning'
                label 'fast'
                model 'openai/gpt-5-mini'

                input:
                question: String

                output:
                answer: String

                prompt:
                """
                Question: ${question}
                """
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should accept record-typed agent I/O'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                answer: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()

                input:
                q: Question

                output:
                a: Answer

                prompt:
                """
                ${q.text}
                """
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should accept val agent I/O'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()

                input:
                question: String

                output:
                answer: String

                prompt:
                """
                ${question}
                """
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should reject destructured record agent I/O'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            record Answer {
                answer: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()

                input:
                record(text: String)

                output:
                a: Answer

                prompt:
                """
                ${text}
                """
            }
            ''')

        then:
        !errors.isEmpty()
        errors.any { it.getOriginalMessage().contains('named record type') }
    }

    def 'should reject a tuple agent input'() {
        when:
        // a tuple input declares no context slot for its components, so an agent would
        // half-ignore it: nothing in the input JSON and nothing staged
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                tuple(a: Integer, b: Path)

                output:
                answer: String

                prompt:
                "go"
            }
            ''')

        then:
        errors.any { it.getOriginalMessage().contains('tuple inputs are not supported') }
        and: 'the message identifies WHICH input, by its components -- a tuple parameter has no name'
        errors.any { it.getOriginalMessage().contains('Agent input `tuple(a, b)`') }
        and: 'the message says what to do instead'
        errors.any { it.getOriginalMessage().contains('separate input') }
    }

    def 'should reject a bare expression as an agent output'() {
        when:
        // the shared `processOutput` rule admits a bare expression (a process lowers it to `$out`),
        // and an agent has no such thing -- so it must be a diagnostic, not a dropped statement
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                file('report.md')

                prompt:
                "go"
            }
            ''')

        then:
        errors.any { it.getOriginalMessage().contains('Agent output must be declared as `name: Type`') }
    }

    def 'should reject more than one agent output'() {
        when:
        // unlike a typed process (a warning there), an agent MUST have a single output: the model
        // answers one value, so multiple outputs have no lowering -- combine them into a record
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                answer: String
                score: Integer

                prompt:
                "go"
            }
            ''')

        then:
        errors.any { it.getOriginalMessage().contains('Agent should have only one output') }
    }

    def 'should resolve file/files in an agent output but not the process-only directives'() {
        when: 'file() is in scope'
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                report: Path = file('report.md')

                prompt:
                "go"
            }
            ''')

        then:
        errors.isEmpty()

        when: 'files() is in scope'
        errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                notes: Set<Path> = files('*.txt')

                prompt:
                "go"
            }
            ''')

        then:
        errors.isEmpty()

        when: 'stdout() is process-only, so it must not resolve in an agent output'
        errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                answer: String = stdout()

                prompt:
                "go"
            }
            ''')

        then:
        errors.any { it.getOriginalMessage().contains('stdout') }
    }
}
