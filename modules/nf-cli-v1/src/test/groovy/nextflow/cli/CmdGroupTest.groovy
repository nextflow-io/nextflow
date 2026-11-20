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
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

class CmdGroupTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should list the sub-commands in the usage'() {
        when:
        new CmdPipeline().usage([])
        then:
        def out = capture.toString()
        out.contains('Usage: nextflow pipeline <command> [options]')
        out.contains('  install     Install a pipeline from the registry')
        out.contains('  view        Show pipeline information')
    }

    def 'should fail on an unknown sub-command'() {
        when:
        new CmdModule().usage(['foo'])
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Unknown module sub-command: foo'
    }

    def 'should reject unknown sub-command arguments with a hint'() {
        given:
        def cmd = new CmdPipeline(args: ['foo'], unknownOptions: [])
        cmd.launcher = new Launcher()
        when:
        cmd.run()
        then:
        def e = thrown(AbortOperationException)
        e.message.endsWith("-- Check the available commands and options and syntax with 'nextflow pipeline -h'")
    }
}
