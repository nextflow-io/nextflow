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

package nextflow.script.control

import nextflow.script.ast.ASTNodeMarker
import nextflow.script.ast.AssignmentExpression
import nextflow.script.ast.ProcessNodeV2
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * An agent's typed I/O lowers through the SAME compiler units a process's does, so the implicit
 * stagers and the output unstagers must come out identical for identical declarations.
 *
 * @see nextflow.script.control.AgentToGroovyVisitor
 * @see nextflow.script.control.ImplicitStagers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentToGroovyTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    def 'should generate an implicit stager for each Path-typed agent input'() {
        when:
        final blocks = lowerAgent("""
            agent qa {
                input:
                ${declaration}
                output:
                answer: String
                prompt: "go"
            }
            """)

        then:
        stagerTargets(blocks.stagers) == expected

        where:
        declaration                     || expected
        'contigs: Path'                 || ['contigs']
        'contigs: Path?'                || ['contigs']
        'reads: Set<Path>'              || ['reads']
        'reads: List<Path>'             || ['reads']
        'n: Integer'                    || []
        's: String'                     || []
        'contigs: Path\nn: Integer'     || ['contigs']
    }

    def 'should generate the same stagers an equivalent process generates'() {
        given:
        final text = '''
            nextflow.enable.types = true

            record Sample {
                id: String
                seq: Path
                index: Path
            }

            process CHECK {
                input:
                sample: Sample
                n: Integer
                output:
                answer: String = 'x'
                script:
                'true'
            }

            agent qa {
                input:
                sample: Sample
                n: Integer
                output:
                answer: String
                prompt: "go"
            }
            '''
        final source = analyze(text)

        when: 'the record recursion runs for both'
        final agent = new AgentToGroovyVisitor(source).transform(source.getAST().getAgents().first())
        final process = new ProcessToGroovyVisitorV2(source)
            .transform((ProcessNodeV2) source.getAST().getProcesses().first())

        then:
        stagerTargets(blockAt(agent, 1)) == ['sample.seq', 'sample.index']
        and: 'the process stagers block is byte-identical in shape'
        stagerTargets(blockAt(process, 1)) == stagerTargets(blockAt(agent, 1))
    }

    def 'should carry the nullable marker into the input declaration'() {
        when:
        final blocks = lowerAgent('''
            agent qa {
                input:
                a: Path
                b: Path?
                output:
                answer: String
                prompt: "go"
            }
            ''')

        then: 'the third argument of _input_ is the optional flag'
        callArgs(blocks.inputs) == [['a', 'java.nio.file.Path', false], ['b', 'java.nio.file.Path', true]]
    }

    def 'should lower a file output into an unstager plus a value-carrying output'() {
        when:
        final blocks = lowerAgent('''
            agent qa {
                input:
                q: String
                output:
                answer: String
                report: Path = file('report.md')
                notes: Set<Path> = files('*.txt')
                prompt: "go"
            }
            ''')

        then: 'one unstager per file()/files() call, keyed per agent'
        callNames(blocks.unstagers) == ['_unstage_files', '_unstage_files']
        unstagerKeys(blocks.unstagers) == ['$path0', '$path1']

        and: 'a bare output stays 2-arg (the model answers it); an RHS output carries its closure'
        blocks.outputs.statements.collect { arity(it) } == [2, 3, 3]
    }

    def 'should not lower env or eval in an agent output'() {
        when:
        final blocks = lowerAgent('''
            agent qa {
                input:
                q: String
                output:
                answer: String = env('HOME')
                prompt: "go"
            }
            ''')

        then: 'no unstager is generated: an agent has no task environment to read back'
        blocks.unstagers.statements.isEmpty()
    }

    def 'should place the stagers and unstagers before the prompt statement'() {
        when:
        final blocks = lowerAgent('''
            agent qa {
                input:
                contigs: Path
                output:
                report: Path = file('report.md')
                prompt: "go"
            }
            ''')

        then: 'the prompt must stay LAST: its value is the closure return value'
        blocks.order == ['directives', 'stagers', 'unstagers', 'inputs', 'outputs', 'prompt']
        and:
        stagerTargets(blocks.stagers) == ['contigs']
        callNames(blocks.unstagers) == ['_unstage_files']
    }

    def 'should lower every file()/files() call form a process output accepts'() {
        given: """the one-arg and two-arg opts forms, singular and plural. DSL resolution matches by
                  NAME only, so this covers the LOWERING rather than the declared signature"""
        final blocks = lowerAgent('''
            nextflow.enable.types = true

            agent qa {
                input:
                x: String
                output:
                one:   Path      = file('report.md')
                typed: Path      = file(type: 'file', 'report.md')
                many:  Set<Path> = files('*.tsv')
                opts:  Set<Path> = files(hidden: true, '*.log')
                prompt: "go"
            }
            ''')

        expect: 'all four resolve, and each lowers to its own unstager'
        callNames(blocks.unstagers) == ['_unstage_files'] * 4
    }

    def 'should resolve file() in an agent output to the same scope a process output resolves it to'() {
        given:
        final source = analyze('''
            nextflow.enable.types = true

            process CHECK {
                input:
                x: String
                output:
                report: Path = file('report.md')
                script:
                "true"
            }

            agent qa {
                input:
                x: String
                output:
                report: Path = file('report.md')
                prompt: "go"
            }
            ''')
        final ast = source.getAST()

        when:
        final agentTarget = fileCallTarget(ast.getAgents().first().outputs)
        final processTarget = fileCallTarget(((ProcessNodeV2) ast.getProcesses().first()).outputs)

        then: 'the shared declaration -- NOT the driver-side global `ScriptDsl.file`, which also matches by name'
        agentTarget == 'nextflow.script.dsl.FileOutputDsl'
        and: 'and provably one declaration, not two that may drift'
        processTarget == agentTarget
    }

    // -----------------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------------

    private analyze(String contents) {
        scriptParser.compiler().getSources().clear()
        final source = scriptParser.parse('main.nf', contents.stripIndent())
        scriptParser.analyze()
        final errors = TestUtils.getErrors(source)
        assert errors.isEmpty() : errors.collect { it.getOriginalMessage() + ' @' + it.getStartLine() }.join('; ')
        return source
    }

    /** Lower the (single) agent and name the six blocks of the generated closure body. */
    private Map lowerAgent(String contents) {
        final source = analyze(contents)
        final stmt = new AgentToGroovyVisitor(source).transform(source.getAST().getAgents().first())
        return [
            order    : ['directives', 'stagers', 'unstagers', 'inputs', 'outputs', 'prompt'],
            stagers  : blockAt(stmt, 1),
            unstagers: blockAt(stmt, 2),
            inputs   : blockAt(stmt, 3),
            outputs  : blockAt(stmt, 4),
        ]
    }

    /** The i-th statement of the generated `agent('name') { ... }` closure body. */
    private static BlockStatement blockAt(org.codehaus.groovy.ast.stmt.Statement stmt, int index) {
        final call = (MethodCallExpression) ((ExpressionStatement) stmt).getExpression()
        final closure = (ClosureExpression) call.getArguments().getExpression(1)
        return (BlockStatement) ((BlockStatement) closure.getCode()).getStatements().get(index)
    }

    /** The `stageAs({ <target> })` targets, rendered as source-like text. */
    private static List<String> stagerTargets(BlockStatement block) {
        return block.getStatements().collect { st ->
            final call = (MethodCallExpression) ((ExpressionStatement) st).getExpression()
            assert call.getMethodAsString() == 'stageAs'
            final closure = (ClosureExpression) call.getArguments().getExpression(0)
            return render(closureBody(closure))
        }
    }

    /** The single expression of a one-statement closure, wrapped in a block or not. */
    private static closureBody(ClosureExpression closure) {
        final code = closure.getCode()
        final st = code instanceof BlockStatement ? code.getStatements().get(0) : code
        return ((ExpressionStatement) st).getExpression()
    }

    private static String render(expr) {
        if( expr instanceof VariableExpression )
            return expr.getName()
        if( expr instanceof PropertyExpression )
            return render(expr.getObjectExpression()) + '.' + expr.getPropertyAsString()
        return expr.getText()
    }

    private static List<String> callNames(BlockStatement block) {
        return block.getStatements().collect { st ->
            ((MethodCallExpression) ((ExpressionStatement) st).getExpression()).getMethodAsString()
        }
    }

    private static List<String> unstagerKeys(BlockStatement block) {
        return block.getStatements().collect { st ->
            final call = (MethodCallExpression) ((ExpressionStatement) st).getExpression()
            return ((ConstantExpression) call.getArguments().getExpression(0)).getValue()
        }
    }

    private static int arity(st) {
        final call = (MethodCallExpression) ((ExpressionStatement) st).getExpression()
        return call.getArguments().getExpressions().size()
    }

    /** The (name, type, optional) triples of the generated `_input_` calls. */
    private static List callArgs(BlockStatement block) {
        return block.getStatements().collect { st ->
            final call = (MethodCallExpression) ((ExpressionStatement) st).getExpression()
            final args = call.getArguments().getExpressions()
            return [
                ((ConstantExpression) args[0]).getValue(),
                args[1].getType().getName(),
                ((ConstantExpression) args[2]).getValue() ]
        }
    }

    /**
     * The declaring class the `file(...)` call in the FIRST output declaration resolved to.
     * `file` is overloaded, so the scope visitor records METHOD_OVERLOADS rather than a single
     * METHOD_TARGET; both spellings are read so the assertion does not depend on the arity used.
     */
    private static String fileCallTarget(Statement outputs) {
        final st = ((BlockStatement) outputs).getStatements().first()
        final assign = (AssignmentExpression) ((ExpressionStatement) st).getExpression()
        final call = (MethodCallExpression) assign.getRightExpression()
        final overloads = (List<MethodNode>) call.getNodeMetaData(ASTNodeMarker.METHOD_OVERLOADS)
        final mn = overloads ? overloads.first() : (MethodNode) call.getNodeMetaData(ASTNodeMarker.METHOD_TARGET)
        return mn?.getDeclaringClass()?.getName()
    }

}
