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

package nextflow.script.parser

import nextflow.script.control.ScriptParser
import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * @see nextflow.script.parser.ScriptAstBuilder
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptAstBuilderTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    List<SyntaxException> check(String contents) {
        return TestUtils.check(scriptParser, contents)
    }

    def 'should report an error for invalid syntax' () {
        when:
        def errors = check(
            '''\
            workflow {
                foo(1, 2 3, 4)
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 14
        errors[0].getOriginalMessage() == "Unexpected input: '3'"

        when:
        errors = check(
            '''\
            process hello {
                debug true:

                script:
                """
                echo hello
                """
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 16
        errors[0].getOriginalMessage() == "Unexpected input: '\\n'"

        when:
        errors = check(
            '''\
            process hello {
                output:
                tuple val("hello") val("goodbye")

                script:
                """
                echo hello
                """
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 24
        errors[0].getOriginalMessage() == "Unexpected input: 'val'"

        when:
        errors = check(
            '''\
            process hello {
                output:
                tuple val('hello'), val('goodbye') emit: message

                script:
                """
                echo hello
                """
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 40
        errors[0].getOriginalMessage() == "Unexpected input: 'emit'"

        when:
        errors = check(
            '''\
            workflow foo {
                take:
                x
                y,
                z

                main:
                println 'hello'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 4
        errors[0].getStartColumn() == 6
        errors[0].getOriginalMessage() == "Unexpected input: ','"
    }

    def 'should report an error for mixing script declarations with statements' () {
        when:
        def errors = check(
            '''\
            println 'Hello world!'

            workflow {
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage().contains "Statements cannot be mixed with script declarations"
    }

    def 'should report an error for implicit process script section' () {
        when:
        def errors = check(
            '''\
            process hello {
                input:
                val greeting

                """
                echo '${greeting}!'
                """
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 5
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "The `script:` or `exec:` label is required when other sections are present"
    }

    def 'should report an error for invalid process directive' () {
        when:
        def errors = check(
            '''\
            process hello {
                container = 'ubuntu'

                script:
                ""
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Invalid process directive"

        when:
        errors = check(
            '''\
            process hello {
                if( true )
                    container = 'ubuntu'
                else
                    container = 'centos'

                script:
                ""
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Invalid process directive"
    }

    def 'should report an error for invalid workflow emits' () {
        when:
        def errors = check(
            '''\
            workflow hello {
                emit:
                'hello'
                target = 'world'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Every emit must be assigned to a name when there are multiple emits"

        when:
        errors = check(
            '''\
            workflow hello {
                emit:
                if( true )
                    'hello'
                else
                    'world'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Invalid workflow emit -- must be a name, assignment, or expression"
    }

    def 'should report error for defining a typed process without preview flag' () {
        when:
        def errors = check(
            '''\
            process hello {
                input:
                message: String

                output:
                result: String

                exec:
                result = message
            }
            '''
        )
        then:
        errors.size() >= 2
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Typed input declaration is not allowed in legacy process -- set `nextflow.preview.types = true` to use typed processes in this script"
        errors[1].getStartLine() == 6
        errors[1].getStartColumn() == 5
        errors[1].getOriginalMessage() == "Typed output declaration is not allowed in legacy process -- set `nextflow.preview.types = true` to use typed processes in this script"

        when:
        errors = check(
            '''\
            nextflow.preview.types = true

            process hello {
                input:
                message: String

                output:
                result: String

                exec:
                result = message
            }
            '''
        )
        then:
        errors.size() == 0
    }

    def 'should report error for defining a legacy process with preview flag enabled' () {
        when:
        def errors = check(
            '''\
            nextflow.preview.types = true

            process hello {
                input:
                val message

                output:
                val result

                exec:
                result = message
            }
            '''
        )
        then:
        errors.size() >= 1
        errors[0].getStartLine() == 5
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Invalid input declaration in typed process"

        when:
        errors = check(
            '''\
            process hello {
                input:
                val message

                output:
                val result

                exec:
                result = message
            }
            '''
        )
        then:
        errors.size() == 0
    }

    def 'should report error for invalid topic statement' () {
        when:
        def errors = check(
            '''\
            nextflow.preview.types = true

            process hello {
                topic:
                versions = stdout()

                script:
                ""
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 5
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Invalid process topic statement"

        when:
        errors = check(
            '''\
            nextflow.preview.types = true

            process hello {
                topic:
                stdout() >> 'versions'

                script:
                ""
            }
            '''
        )
        then:
        errors.size() == 0
    }

    def 'should report error for invalid record definition' () {
        when:
        def errors = check(
            '''\
            record Sample {}
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == "Missing record body"

        when:
        errors = check(
            '''\
            record Sample {
                id
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Missing field type"
    }

}
