/*
 * Copyright 2020-2021, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import spock.lang.Specification

import java.nio.file.Files

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdTestTest extends Specification {

    def 'should create reports' () {

        given:
        def pipeline = Files.createTempDirectory('test')
        def workdir = Files.createTempDirectory('workdir')
        def file = Files.createFile(pipeline.resolve('dummy_test.nf'))
        file.text = """
        nextflow.enable.dsl=2
        
        process foo {
            tag "\$id"
            input: val(id)
            output: stdout
            script:
            "echo -n hello \$id"
        }
        
        testflow test1 {
          when:
            foo( channel.of('A','B') )
          then:
            assert foo.outputsCount() == 1
            assert foo.emissionsCount() == 2
        
            foo.emissionWith(tag:'A') {
                assert out[0] == 'hello A'
            }
        
            foo.emissionWith(tag:'B') {
                assert out[0] == 'hello B'
            }
        }
        """.stripIndent()
        def cmd = new CmdTest(launcher: new Launcher())
        cmd.args = [pipeline.toAbsolutePath().toString()]
        cmd.workDir = workdir.toAbsolutePath().toString()

        when: 'tests are run'
        cmd.run()

        then: 'should create XML reports'
        workdir.resolve('test-results/TEST-dummy_test.xml').exists()

        and: 'should create HTML reports'
        workdir.resolve('test-reports/index.html').exists()
        workdir.resolve('test-reports/dummy_test.html').exists()

        cleanup:
        pipeline?.deleteDir()
        workdir?.deleteDir()

    }

}
