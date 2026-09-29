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
 * @see nextflow.script.control.CallArityVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class CallArityTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    List<SyntaxException> check(String contents) {
        return TestUtils.check(scriptParser, contents)
    }

    def 'should report an error for an invalid method call' () {
        when:
        def errors = check(
            '''\
            process hello {
                script:
                """
                echo hello
                """
            }

            workflow {
                hello('world')
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 9
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Incorrect number of call arguments, expected 0 but received 1'

        when:
        errors = check(
            '''\
            workflow hello {
                take:
                x
                y

                main:
                println 'hello'
            }

            workflow {
                hello('world')
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 11
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Incorrect number of call arguments, expected 2 but received 1'
    }

    def 'should report an error for an agent call with the wrong number of arguments' () {
        when:
        def errors = check(
            '''\
            agent reporter {
                model 'openai/gpt-4o'
                instruction 'i'

                input:
                sample: String

                output:
                report: String

                prompt:
                """
                Report on ${sample}.
                """
            }

            workflow {
                reporter('a', 'b')
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Incorrect number of call arguments, expected 1 but received 2'

        when: 'the call arity matches the declared inputs'
        errors = check(
            '''\
            agent reporter {
                model 'openai/gpt-4o'
                instruction 'i'

                input:
                sample: String

                output:
                report: String

                prompt:
                """
                Report on ${sample}.
                """
            }

            workflow {
                reporter('a')
            }
            '''
        )
        then:
        errors.size() == 0
    }

}
