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

import java.nio.file.Files

import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class CmdLintTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should report compilation errors' () {

        given:
        def dir = Files.createTempDirectory('test')

        dir.resolve('main.nf').text = '''\
            process HELLO {

                script:
                """
                ${
                    params.is_paired_end
                        ? "..."
                        : "..."
                }
                """
            }
            '''

        dir.resolve('nextflow.config').text = '''\
            process {
                withLabel:
                'bambino' {
                    container = "..."
                }
            }
            '''

        when:
        def cmd = new CmdLint()
        cmd.args = [dir.toFile().toString()]
        cmd.launcher = Mock(Launcher) {
            getOptions() >> Mock(CliOptions)
        }
        cmd.run()

        then:
        thrown(AbortOperationException)
        and:
        capture.toString().contains("Error $dir/main.nf:5:19: Unexpected input: '\\n'")
        capture.toString().contains("│   5 |                 \${")
        capture.toString().contains("╰     |                   ^")
        capture.toString().contains("Error $dir/nextflow.config:2:27: Unexpected input: '\\n'")
        capture.toString().contains("│   2 |                 withLabel:")
        capture.toString().contains("╰     |                           ^")

        cleanup:
        dir?.deleteDir()
    }

}
