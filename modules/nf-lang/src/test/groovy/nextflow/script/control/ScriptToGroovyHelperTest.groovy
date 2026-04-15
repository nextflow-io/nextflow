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

import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptToGroovyHelperTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    def parse(String contents) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents.stripIndent())
        assert !TestUtils.hasSyntaxErrors(source)
        return source
    }

    def 'should get source text of process body' () {
        given:
        def source = parse('''\
            process hello {
                script:
                """
                echo 'hello!'
                """
            }
            ''')
        def sgh = new ScriptToGroovyHelper(source)

        when:
        def process = source.getAST().getProcesses().first()
        then:
        sgh.getSourceText(process.exec) == '    """\n    echo \'hello!\'\n    """\n'
    }

}
