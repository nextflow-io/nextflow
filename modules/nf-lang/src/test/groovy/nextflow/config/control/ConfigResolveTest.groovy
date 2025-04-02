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

package nextflow.config.control

import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Ignore
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * @see nextflow.config.control.ConfigResolveVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigResolveTest extends Specification {

    @Shared
    ConfigParser parser

    def setupSpec() {
        parser = new ConfigParser()
    }

    List<SyntaxException> check(String contents) {
        return TestUtils.check(parser, contents)
    }

    def 'should report an error for an undefined variable' () {
        when:
        def errors = check(
            '''\
            process.clusterOptions = "--cpus ${process.cpus}"
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 36
        errors[0].getOriginalMessage() == '`process` is not defined'
    }

    def 'should report an error for an invalid dynamic config option' () {
        when:
        def errors = check(
            '''\
            report.file = { "report.html" }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 1
        errors[0].getOriginalMessage() == 'Dynamic config options are only allowed in the `process` scope'

        when:
        errors = check(
            '''\
            process.clusterOptions = { "--cpus ${task.cpus}" }
            '''
        )
        then:
        errors.size() == 0
    }

    @Ignore("config includes aren't working in tests due to FileSystemNotFoundException")
    def 'should report an error for an invalid config include' () {
        when:
        def errors = check(
            '''\
            profiles {
                includeConfig 'foo.config'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == 'Config includes are only allowed at the top-level or in a profile'

        when:
        errors = check(
            '''\
            profiles {
                foo {
                    includeConfig 'foo.config'
                }
            }
            '''
        )
        then:
        errors.size() == 0
    }

}
