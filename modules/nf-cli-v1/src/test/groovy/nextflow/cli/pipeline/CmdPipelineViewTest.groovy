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
import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineMetadata
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.client.RegistryClient
import nextflow.cli.Launcher
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

class CmdPipelineViewTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    RegistryClient client = Mock(RegistryClient)

    static final Map SCHEMA = [
        '$defs'   : [
            input_output: [properties: [
                input : [type: 'string', description: 'Samplesheet.\nCSV.'],
                outdir: [description: 'Output dir (no type)']
            ]]
        ],
        properties: [flag: [type: 'boolean', description: 'A top-level flag']]
    ]

    static PipelineRelease release(String version, Map schema = SCHEMA) {
        def r = new PipelineRelease().version(version).gitRevision(version).commitId('abc123')
        r.schema = schema
        return r
    }

    static Pipeline pipeline() {
        new Pipeline().latest(release('3.14.0')).repositoryUrl('https://github.com/nf-core/rnaseq')
            .metadata(new PipelineMetadata().keywords(['rna', 'seq']))
    }

    CmdPipelineView cmd(String output = 'text') {
        def cmd = new CmdPipelineView()
        cmd.args = ['nf-core/rnaseq']
        cmd.output = output
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.client = client
        return cmd
    }

    def 'should print pipeline info, params and no usage'() {
        when:
        cmd().run()
        def out = capture.toString()

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline()
        0 * client.getPipelineRelease(_, _)
        out.contains('Name:        nf-core/rnaseq')
        out.contains('Version:     3.14.0')
        out.contains('Repository:  https://github.com/nf-core/rnaseq')
        out.contains('Revision:    3.14.0')
        out.contains('Commit:      abc123')
        out.contains('Keywords:    rna, seq')
        out.contains('/pipelines/nf-core/rnaseq@3.14.0')
        out.contains('--input (string)')
        out.contains('Samplesheet. CSV.')
        out.contains('--outdir\n')
        out.contains('--flag (boolean)')
        !out.contains('Usage:')
    }

    def 'should fetch a given version'() {
        given:
        def c = cmd()
        c.version = '3.13.0'

        when:
        c.run()

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline()
        1 * client.getPipelineRelease('nf-core/rnaseq', '3.13.0') >> release('3.13.0', null)
        capture.toString().contains('Version:     3.13.0')
        !capture.toString().contains('Parameters:')
    }

    def 'should print json'() {
        when:
        cmd('json').run()
        def lines = capture.toString().readLines()
        def json = new JsonSlurper().parseText(lines[lines.findIndexOf { it.trim().startsWith('{') }..-1].join('\n')) as Map

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline()
        json.version == '3.14.0'
        json.commitId == 'abc123'
        json.repositoryUrl == 'https://github.com/nf-core/rnaseq'
        json.params*.name == ['input', 'outdir', 'flag']
        !json.containsKey('usage')
    }

}
