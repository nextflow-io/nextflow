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

package nextflow.cli.pipeline

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import nextflow.pipeline.PipelineInfo
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture
import org.junit.Rule

class CmdPipelineListTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path base

    Path seed() {
        def dir = base.resolve('pipelines/nf-core/rnaseq')
        Files.createDirectories(dir)
        Files.writeString(dir.resolve('main.nf'), 'x')
        new PipelineInfo(name: 'nf-core/rnaseq', version: '3.14.0', checksum: PipelineInfo.checksum(dir)).save(dir)
        return dir
    }

    def 'should list installed pipelines as text'() {
        given:
        seed()
        def cmd = new CmdPipelineList(root: base)

        when:
        cmd.run()

        then:
        def lines = capture.toString().readLines()
        lines[0] == 'Pipeline'.padRight(40) + 'Version'.padRight(15) + 'Status'
        lines[1] == '-' * 65
        lines[2] == 'nf-core/rnaseq'.padRight(40) + '3.14.0'.padRight(15) + 'ok'
    }

    def 'should flag modified pipelines in json'() {
        given:
        Files.writeString(seed().resolve('main.nf'), 'changed')
        def cmd = new CmdPipelineList(root: base, output: 'json')

        when:
        cmd.run()
        def json = new JsonSlurper().parseText(capture.toString()) as Map

        then:
        json.pipelines == [[name: 'nf-core/rnaseq', version: '3.14.0', status: 'modified']]
    }

    def 'should list a pipeline with a damaged .pipeline-info'() {
        given:
        def dir = base.resolve('pipelines/acme/bare')
        Files.createDirectories(dir)
        Files.writeString(dir.resolve('.pipeline-info'), 'registryUrl=x\n')
        def cmd = new CmdPipelineList(root: base)

        when:
        cmd.run()

        then:
        capture.toString().readLines()[2] == '-'.padRight(40) + '-'.padRight(15) + 'modified'
    }
}
