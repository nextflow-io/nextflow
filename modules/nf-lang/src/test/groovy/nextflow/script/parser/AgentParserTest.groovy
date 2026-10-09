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
                stdout()

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
                stdout()

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
        errors[0].getStartLine() == 14
        errors[0].getOriginalMessage() == "Unexpected input: '}'"
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
                stdout()

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
                stdout()

                prompt:
                """
                ${question}
                """
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should accept destructured agent inputs'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                record(id: String, reads: Path)
                tuple(a: Integer, b: Path)

                output:
                record(id: id, a: a, summary: stdout())

                prompt:
                "Summarize ${id} ${reads} ${a} ${b}"
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should accept a composed agent output'() {
        when:
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
        errors.isEmpty()
    }

    def 'should reject multiple agent outputs'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                count: Integer
                report: Path = file('report.md')

                prompt:
                "go"
            }
            ''')

        then:
        errors.any { it.getOriginalMessage() == 'Agent should have only one output -- combine outputs into a record' }
    }

    def 'should allow a named output when it is the only agent output'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                summary: String = stdout()

                prompt:
                "go"
            }
            ''')

        then:
        errors.isEmpty()
    }

    def 'should check the types of model-answered agent outputs'() {
        when:
        def errors = check("""\
            nextflow.enable.types = true

            record Answer {
                text: String
            }

            record BadAnswer {
                text: String
                meta: Map
            }

            agent qa {
                input:
                q: String

                output:
                ${declaration}

                prompt:
                "go"
            }
            """)

        then:
        errors.collect { it.getOriginalMessage() }.findAll { it.startsWith('Agent output') } == expected

        where:
        declaration                     || expected
        'a: Answer'                     || []
        'n: Integer'                    || []
        'x: Float'                      || []
        'ok: Boolean'                   || []
        'f: Path'                       || []
        's: String'                     || []
        'items: List<String>'           || []
        'answers: List<Answer>'         || []
        'nested: List<List<Path>>'      || []
        'items: List'                   || ['Agent output `items` has unsupported type List -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type']
        'items: List<Map>'              || ['Agent output `items` has unsupported type List<Map> -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type']
        'd: Double'                     || ['Agent output `d` has unsupported type Double -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type']
        'm: Map'                        || ['Agent output `m` has unsupported type Map -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type']
        'b: BadAnswer'                  || ['Agent output `b` has unsupported field `meta` with type Map -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type']
        'bs: List<BadAnswer>'           || ['Agent output `bs` has unsupported field `meta` with type Map -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type']
        'q'                             || ['Agent output `q` should declare a type -- typed outputs are answered by the model']
    }

    def 'should resolve file/files/stdout in an agent output but not the process-only directives'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                record(
                    summary: stdout(),
                    report: file('report.md'),
                    notes: files('*.txt')
                )

                prompt:
                "go"
            }
            ''')

        then:
        errors.isEmpty()

        when: 'eval() is process-only, so it must not resolve in an agent output'
        errors = check('''\
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                eval('date')

                prompt:
                "go"
            }
            ''')

        then:
        errors.any { it.getOriginalMessage().contains('eval') }
    }
}
