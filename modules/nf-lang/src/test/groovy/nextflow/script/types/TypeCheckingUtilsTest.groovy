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

import java.nio.file.Path

import nextflow.script.control.ScriptParser
import nextflow.script.dsl.Types
import nextflow.script.namespaces.ChannelNamespace
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.GenericsType
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.ast.tools.GenericsUtils
import org.codehaus.groovy.control.SourceUnit
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll

import static nextflow.script.types.TypeCheckingUtils.*

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingUtilsTest extends Specification {

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

    boolean checkType(Expression exp, Class type) {
        assert Types.isEqual(getType(exp), ClassHelper.makeCached(type))
        return true
    }

    @Unroll
    def 'should resolve the type of a literal expression' () {
        expect:
        checkType(parseExpression(EXPR), TYPE)

        where:
        EXPR                    | TYPE
        'true'                  | Boolean
        '42'                    | Integer
        '3.14'                  | Float
        "'hello'"               | String
        '"hello ${"world!"}"'   | String
        '[1, 2, 3]'             | List
        '[a: 1, b: 2, c: 3]'    | Map
    }

    @Unroll
    def 'should resolve the return type of a function call' () {
        expect:
        checkType(parseExpression(EXPR), TYPE)

        where:
        EXPR            | TYPE
        "env('HELLO')"  | String
    }

    @Unroll
    def 'should resolve the type of a property access' () {
        expect:
        checkType(parseExpression(EXPR), TYPE)

        where:
        EXPR                        | TYPE
        'workflow.projectDir'       | Path
        'workflow.projectDir.name'  | String
    }

    @Unroll
    def 'should resolve the return type of a method call' () {
        expect:
        checkType(parseExpression(EXPR), TYPE)

        where:
        EXPR                            | TYPE
        'workflow.projectDir.exists()'  | Boolean
    }

    def 'should resolve the return type of a process invocation' () {
        when:
        def exp = parseExpression(
            '''
            nextflow.enable.types = true

            process hello {
                input:
                target: String

                output:
                x = 'hello'
                y = target

                script:
                """
                """
            }

            workflow {
                hello('world')
            }
            '''
        )
        def cn = getType(exp)
        then:
        cn.getField('x') != null
        cn.getField('y') != null
    }

    def 'should resolve the type of a workflow output' () {
        when:
        def exp = parseExpression(
            '''
            nextflow.enable.types = true

            workflow hello {
                take:
                target: String

                emit:
                x = 'hello'
                y = target
            }

            workflow {
                hello('world')
            }
            '''
        )
        def cn = getType(exp)
        then:
        cn.getField('x') != null
        cn.getField('y') != null
    }

    def 'should extract placeholder mapping from a parameterized type' () {
        when:
        def cn = CHANNEL_TYPE
        def spec = GenericsUtils.extractPlaceholders(cn)
        then:
        spec == [:]

        when:
        cn = makeType(CHANNEL_TYPE, ClassHelper.Integer_TYPE)
        spec = GenericsUtils.extractPlaceholders(cn)
        then:
        spec[new GenericsType.GenericsTypeName('E')]?.getType() == ClassHelper.Integer_TYPE
    }

    def 'should apply placeholder mapping to a parameterized type' () {
        given:
        def declaringType = makeType(CHANNEL_TYPE, ClassHelper.Integer_TYPE)

        when:
        def cn = resolveGenericType(CHANNEL_TYPE, declaringType)
        def spec = GenericsUtils.extractPlaceholders(cn)
        then:
        spec[new GenericsType.GenericsTypeName('E')]?.getType() == ClassHelper.Integer_TYPE

        when:
        def mn = ClassHelper.makeCached(ChannelNamespace).getDeclaredMethods('of').first()
        def returnType = resolveGenericType(mn.getReturnType(), declaringType)
        spec = GenericsUtils.extractPlaceholders(returnType)
        then:
        spec[new GenericsType.GenericsTypeName('E')]?.getType() == ClassHelper.Integer_TYPE
    }

    private static final ClassNode CHANNEL_TYPE = ClassHelper.makeCached(Channel)

}
