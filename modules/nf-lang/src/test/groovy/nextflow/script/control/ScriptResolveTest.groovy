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

import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * @see nextflow.script.control.ScriptResolveVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptResolveTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    List<SyntaxException> check(String contents) {
        return TestUtils.check(scriptParser, contents)
    }

    def 'should report an error for conflicting declarations' () {
        when:
        def errors = check(
            '''\
            def hello() {
            }

            workflow hello {
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 4
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == '`hello` is already declared'

        when:
        errors = check(
            '''\
            def hello(x) {
                def x = 'world'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == '`x` is already declared'

        when:
        errors = check(
            '''\
            workflow hello {
                take:
                x

                main:
                def x = 'world'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 6
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == '`x` is already declared'
    }

    def 'should report an error for an unrecognized process directive' () {
        when:
        def errors = check(
            '''\
            process hello {
                cpu 8

                input:
                stdout

                output:
                stdin

                script:
                ""
            }
            '''
        )
        then:
        errors.size() == 3
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Unrecognized process directive `cpu`'
        errors[1].getStartLine() == 5
        errors[1].getStartColumn() == 5
        errors[1].getOriginalMessage() == 'Unrecognized process input qualifier `stdout`'
        errors[2].getStartLine() == 8
        errors[2].getStartColumn() == 5
        errors[2].getOriginalMessage() == 'Unrecognized process output qualifier `stdin`'
    }

    def 'should allow stdin and stdout to be used without parens' () {
        when:
        def errors = check(
            '''\
            process hello {
                input:
                tuple val(id), stdin

                output:
                tuple val(id), stdout

                script:
                ""
            }
            '''
        )
        then:
        errors.size() == 0
    }

    def 'should report an error for an undefined variable' () {
        when:
        def errors = check(
            '''\
            println(x)
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == '`x` is not defined'

        when:
        errors = check(
            '''\
            def hello() {
                x = 'hello'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == '`x` was assigned but not declared'

        when:
        errors = check(
            '''\
            workflow {
                [1, 2, 3].each { i ->
                    x = 'hello'
                }
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == 'Variables in a closure should be declared with `def`'
    }

    def 'should report an error for an undefined function' () {
        when:
        def errors = check(
            '''\
            hello('world')
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == '`hello` is not defined'
    }

    def 'should report an error for an undefined type' () {
        when:
        def errors = check(
            '''\
            'hello' as FooBar
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == '`FooBar` is not defined'

        when:
        errors = check(
            '''\
            [] as List<FooBar>
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == '`FooBar` is not defined'

        when:
        errors = check(
            '''\
            [:] as Map<Foo,Bar>
            '''
        )
        then:
        errors.size() == 2
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == '`Foo` is not defined'
        errors[1].getStartLine() == 1
        errors[1].getStartColumn() == 1
        errors[1].getOriginalMessage() == '`Bar` is not defined'

        when:
        errors = check(
            '''\
            x = new FooBar()
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == '`FooBar` is not defined'
    }

    def 'should report an error for an invalid process invocation' () {
        when:
        def errors = check(
            '''\
            process hello {
                script:
                """
                echo hello
                """
            }

            def greet() {
                hello()
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 9
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Processes can only be called from a workflow'

        when:
        errors = check(
            '''\
            process hello {
                script:
                """
                echo hello
                """
            }

            workflow {
                [1, 2, 3].each { i ->
                    hello()
                }
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 10
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == 'Processes cannot be called from within a closure'
    }

    def 'should report an error for an invalid workflow invocation' () {
        when:
        def errors = check(
            '''\
            workflow hello {
            }

            def greet() {
                hello()
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 5
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Workflows can only be called from a workflow'

        when:
        errors = check(
            '''\
            workflow hello {
            }

            workflow {
                [1, 2, 3].each { i ->
                    hello()
                }
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 6
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == 'Workflows cannot be called from within a closure'
    }

    def 'should recognize fully-qualified class name' () {
        when:
        def errors = check(
            '''\
            new groovy.json.JsonSlurper()
            '''
        )
        then:
        errors.size() == 0

        when:
        errors = check(
            '''\
            new groovy.json.JsonGenerator.Options()
            '''
        )
        then:
        errors.size() == 0
    }

}
