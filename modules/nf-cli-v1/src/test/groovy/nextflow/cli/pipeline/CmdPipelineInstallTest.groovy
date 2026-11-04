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

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.client.RegistryClient
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.pipeline.PipelineGitFetcher
import spock.lang.Specification
import org.junit.Rule
import spock.lang.TempDir
import test.OutputCapture

class CmdPipelineInstallTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path base

    RegistryClient client = Mock(RegistryClient)
    PipelineGitFetcher fetcher = Mock(PipelineGitFetcher)

    static PipelineRelease release(String version) {
        new PipelineRelease().version(version).gitRevision(version).commitId('abc')
    }

    static Pipeline pipeline(String latest = '3.14.0') {
        new Pipeline().latest(release(latest)).repositoryUrl('https://github.com/nf-core/rnaseq').workflowPath('main.nf')
    }

    CmdPipelineInstall cmd(List<String> args, boolean force = false) {
        def cmd = new CmdPipelineInstall()
        cmd.args = args
        cmd.force = force
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.client = client
        cmd.fetcher = fetcher
        cmd.root = base
        return cmd
    }

    def 'should install a specific version'() {
        given:
        def c = cmd(['nf-core/rnaseq'])
        c.version = '3.14.0'

        when:
        c.run()

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline()
        1 * client.getPipelineRelease('nf-core/rnaseq', '3.14.0') >> release('3.14.0')
        1 * fetcher.fetch(_, '3.14.0', 'abc', _) >> { u, r, c2, Path t -> Files.writeString(t.resolve('main.nf'), 'x') }
        Files.isDirectory(base.resolve('pipelines/nf-core/rnaseq'))
        def out = capture.toString()
        def dir = base.resolve('pipelines/nf-core/rnaseq').toAbsolutePath()
        out.indexOf('Installing pipeline nf-core/rnaseq@3.14.0...') >= 0
        out.indexOf('Installing pipeline nf-core/rnaseq@3.14.0...') < out.indexOf("Pipeline nf-core/rnaseq@3.14.0 installed in ${dir}")
    }

    def 'should install the latest version by default'() {
        when:
        cmd(['nf-core/rnaseq']).run()

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline()
        0 * client.getPipelineRelease(_, _)
        1 * fetcher.fetch(_, '3.14.0', 'abc', _) >> { u, r, c2, Path t -> Files.writeString(t.resolve('main.nf'), 'x') }
        Files.isDirectory(base.resolve('pipelines/nf-core/rnaseq'))
    }

    def 'should refuse to overwrite a modified copy without -force'() {
        given:
        client.getPipeline(_) >>> [pipeline('1.0.0'), pipeline('2.0.0')]
        fetcher.fetch(*_) >> { u, r, c2, Path t -> Files.writeString(t.resolve('main.nf'), 'x') }
        cmd(['nf-core/rnaseq']).run()
        Files.writeString(base.resolve('pipelines/nf-core/rnaseq/main.nf'), 'changed')

        when:
        cmd(['nf-core/rnaseq']).run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('-force')
    }

    def 'should report an up to date pipeline instead of installing it again'() {
        given:
        client.getPipeline(_) >> pipeline()
        fetcher.fetch(*_) >> { u, r, c2, Path t -> Files.writeString(t.resolve('main.nf'), 'x') }
        cmd(['nf-core/rnaseq']).run()

        when:
        cmd(['nf-core/rnaseq']).run()

        then:
        capture.toString().contains('Pipeline nf-core/rnaseq@3.14.0 is already installed and up to date')
        capture.toString().count('installed in') == 1
    }
}
