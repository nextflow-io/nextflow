/*
 * Copyright 2024-2025, Seqera Labs
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

package nextflow.script.types

import nextflow.script.control.ScriptParser
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.SourceUnit
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckerTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    SourceUnit parse(String contents) {
        def source = scriptParser.parse('main.nf', contents)
        scriptParser.analyze()
        return source
    }

    Statement parseStatement(String contents) {
        def source = parse(contents)
        assert !source.getErrorCollector().hasErrors()
        def entry = source.getAST().getEntry()
        assert entry.main instanceof BlockStatement
        return entry.main.statements.last()
    }

    Expression parseExpression(String contents) {
        def stmt = parseStatement(contents)
        assert stmt instanceof ExpressionStatement
        return stmt.expression
    }

    Class type(Expression exp) {
        return TypeChecker.getType(exp).getTypeClass()
    }

    def 'should infer the type of a literal expression' () {
        when:
        def exp = parseExpression(
            '''
            true
            '''
        )
        then:
        type(exp) == Boolean

        when:
        exp = parseExpression(
            '''
            42
            '''
        )
        then:
        type(exp) == int

        when:
        exp = parseExpression(
            '''
            3.14
            '''
        )
        then:
        type(exp) == BigDecimal

        when:
        exp = parseExpression(
            '''
            'hello'
            '''
        )
        then:
        type(exp) == nextflow.script.types.shim.String

        when:
        exp = parseExpression(
            '''
            "hello ${'world!'}"
            '''
        )
        then:
        type(exp) == GString

        when:
        exp = parseExpression(
            '''
            [1, 2, 3]
            '''
        )
        then:
        type(exp) == nextflow.script.types.shim.List

        when:
        exp = parseExpression(
            '''
            [a: 1, b: 2, c: 3]
            '''
        )
        then:
        type(exp) == nextflow.script.types.shim.Map
    }

    def 'should infer the return type of a function call' () {
        when:
        def exp = parseExpression(
            '''
            env('HELLO')
            '''
        )
        then:
        type(exp) == nextflow.script.types.shim.String
    }

    def 'should infer the type of a property access' () {
        when:
        def exp = parseExpression(
            '''
            workflow.projectDir
            '''
        )
        then:
        type(exp) == nextflow.script.types.shim.Path
    }

    def 'should infer the type of a variable declaration' () {
        when:
        def exp = parseExpression(
            '''
            def x = 'hello'
            '''
        )
        def target = exp.getLeftExpression()
        then:
        type(target) == nextflow.script.types.shim.String
    }

    def 'should infer the return type of a process invocation' () {
        when:
        def exp = parseExpression(
            '''
            process hello {
                input:
                val(target)

                output:
                val('hello'), emit: x
                val(target), emit: y

                script:
                """
                """
            }

            workflow {
                hello('world')
            }
            '''
        )
        def cn = TypeChecker.getType(exp)
        then:
        cn.getField('x') != null
        cn.getField('y') != null
    }

    def 'should infer the type of a workflow output' () {
        when:
        def exp = parseExpression(
            '''
            workflow hello {
                take:
                target

                emit:
                x = 'hello'
                y = target
            }

            workflow {
                hello('world')
            }
            '''
        )
        def cn = TypeChecker.getType(exp)
        then:
        cn.getField('x') != null
        cn.getField('y') != null
    }

}
