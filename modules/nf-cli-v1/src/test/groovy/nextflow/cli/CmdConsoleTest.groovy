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

package nextflow.cli

import nextflow.exception.AbortOperationException
import spock.lang.Specification

class CmdConsoleTest extends Specification {

    def 'should report that the console command was removed'() {
        when:
        new CmdConsole().run()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'The `console` command has been removed'
    }

    def 'should report removal when usage is requested'() {
        when:
        new CmdConsole().usage()
        then:
        thrown(AbortOperationException)

        when:
        new CmdConsole().usage([])
        then:
        thrown(AbortOperationException)
    }

    def 'should report removal through the launcher'() {
        expect:
        new Launcher().command(args as String[]).run() == 1

        where:
        args                    | _
        ['console']             | _
        ['console', 'foo.nf']   | _
        ['console', '-h']       | _
        ['help', 'console']     | _
    }
}
