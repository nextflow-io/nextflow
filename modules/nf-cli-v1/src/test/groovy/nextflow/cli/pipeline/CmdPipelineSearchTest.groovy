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

import groovy.json.JsonSlurper
import io.seqera.npr.api.schema.v1.ListPipelinesResponse
import io.seqera.npr.api.schema.v1.PipelineMetadata
import io.seqera.npr.api.schema.v1.PipelineSummary
import io.seqera.npr.client.RegistryClient
import nextflow.cli.Launcher
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

class CmdPipelineSearchTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    RegistryClient client = Mock(RegistryClient)

    static ListPipelinesResponse response() {
        new ListPipelinesResponse().pipelines([
            new PipelineSummary().name('nf-core/rnaseq').latestVersion('3.14.0')
                .metadata(new PipelineMetadata().description('RNA-seq analysis'))
        ])
    }

    CmdPipelineSearch cmd(String output = 'simple') {
        def cmd = new CmdPipelineSearch()
        cmd.args = ['rnaseq']
        cmd.limit = 5
        cmd.output = output
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.client = client
        return cmd
    }

    def 'should print name, version and description'() {
        when:
        cmd().run()

        then:
        1 * client.listPipelines('rnaseq', 1, 5) >> response()
        capture.toString().contains("Searching for 'rnaseq'...")
        capture.toString().contains('Top 1 matching pipeline(s):')
        capture.toString().contains('  nf-core/rnaseq\n    Version: 3.14.0\n    Description: RNA-seq analysis\n')
    }

    def 'should print json'() {
        when:
        cmd('json').run()
        def json = new JsonSlurper().parseText(capture.toString()) as List

        then:
        1 * client.listPipelines('rnaseq', 1, 5) >> response()
        json[0].name == 'nf-core/rnaseq'
        json[0].latestVersion == '3.14.0'
        json[0].description == 'RNA-seq analysis'
    }

    def 'should report no results'() {
        when:
        cmd().run()

        then:
        1 * client.listPipelines('rnaseq', 1, 5) >> new ListPipelinesResponse()
        capture.toString().contains('No pipelines found')
    }
}
