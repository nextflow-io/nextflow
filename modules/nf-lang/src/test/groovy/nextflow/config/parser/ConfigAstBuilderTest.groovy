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

package nextflow.config.parser

import nextflow.config.control.ConfigParser
import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * @see nextflow.config.parser.ConfigAstBuilder
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigAstBuilderTest extends Specification {

    @Shared
    ConfigParser parser

    def setupSpec() {
        parser = new ConfigParser()
    }

    List<SyntaxException> check(String contents) {
        return TestUtils.check(parser, contents)
    }

    def 'should report an error for invalid syntax' () {
        when:
        def errors = check(
            '''\
            process {
                cpus = 4
                memory 8.GB
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 12
        errors[0].getOriginalMessage() == "Unexpected input: '8'"

        when:
        errors = check(
            '''\
            if( true ) {
                process.cpus = 8
            }
            else {
                process.cpus = 4
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == "If statements cannot be mixed with config statements"

        when:
        errors = check(
            '''\
            process {
                withName: 'HELLO' {
                    if( true ) {
                        cpus = 8
                    }
                    else {
                        cpus = 4
                    }
                }
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 3
        errors[0].getStartColumn() == 9
        errors[0].getOriginalMessage() == "If statements cannot be mixed with config statements"

        when:
        errors = check(
            '''\
            try {
                process.cpus = 4
            }
            catch( Exception e ) {
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == "Try-catch blocks cannot be mixed with config statements"

        when:
        errors = check(
            '''\
            def trace_timestamp = new Date().format('yyyy-MM-dd_HH-mm-ss')
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == "Variable declarations cannot be mixed with config statements"
    }

}
