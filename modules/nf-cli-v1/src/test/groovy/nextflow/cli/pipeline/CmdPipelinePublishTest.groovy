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

import io.seqera.npr.api.schema.v1.CreatePipelineReleaseResponse
import io.seqera.npr.api.schema.v1.UploadPipelineReleaseResponse
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.cli.Launcher
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.Git
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

class CmdPipelinePublishTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path base

    RegistryClient client = Mock(RegistryClient)

    CmdPipelinePublish cmd(List<String> args = [base.toString()]) {
        def cmd = new CmdPipelinePublish()
        cmd.args = args
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.client = client
        return cmd
    }

    private void pipeline(String manifest = "name = 'myorg/mypipe'; version = '1.2.0'") {
        Files.writeString(base.resolve('nextflow.config'), "manifest { ${manifest} }")
        Files.writeString(base.resolve('main.nf'), 'workflow { }')
    }

    private Git gitCheckout(String origin = 'https://github.com/myorg/mypipe.git') {
        def git = Git.init().setDirectory(base.toFile()).setInitialBranch('main').call()
        git.repository.config.setString('remote', 'origin', 'url', origin)
        git.repository.config.save()
        git.add().addFilepattern('.').call()
        git.commit().setMessage('init').setSign(false).call()
        git.tag().setName('v1.2.0').setAnnotated(false).call()
        return git
    }

    private static CreatePipelineReleaseResponse created() {
        new CreatePipelineReleaseResponse().releaseId(7L)
    }

    def 'should publish bundle and provenance'() {
        given:
        pipeline()
        gitCheckout().close()

        when:
        cmd().run()

        then:
        1 * client.createPipelineRelease('myorg/mypipe', {
            it.version == '1.2.0' && it.checksum ==~ /sha256:.{64}/ && it.size > 0 &&
            it.repositoryUrl == 'https://github.com/myorg/mypipe' && it.gitRevision == 'v1.2.0' && it.commitId
        }) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> new UploadPipelineReleaseResponse()
        def out = capture.toString()
        out.indexOf('Publishing pipeline myorg/mypipe@1.2.0...') >= 0
        out.indexOf('Publishing pipeline myorg/mypipe@1.2.0...') < out.indexOf('Published pipeline myorg/mypipe@1.2.0')
    }

    def 'should let -version override the manifest'() {
        given:
        pipeline()
        def c = cmd()
        c.version = '2.0.0'

        when:
        c.run()

        then:
        1 * client.createPipelineRelease('myorg/mypipe', { it.version == '2.0.0' }) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> new UploadPipelineReleaseResponse()
        capture.toString().contains('Published pipeline myorg/mypipe@2.0.0')
    }

    def 'should publish without git provenance outside git'() {
        given:
        pipeline()

        when:
        cmd().run()

        then:
        1 * client.createPipelineRelease('myorg/mypipe', {
            it.repositoryUrl == null && it.gitRevision == null && it.commitId == null
        }) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> new UploadPipelineReleaseResponse()
    }

    def 'should publish bundle only when the origin has no public url'() {
        given:
        pipeline()
        gitCheckout('/some/local/path').close()

        when:
        cmd().run()

        then:
        1 * client.createPipelineRelease('myorg/mypipe', {
            it.repositoryUrl == null && it.gitRevision == null && it.commitId == null
        }) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> new UploadPipelineReleaseResponse()
        capture.toString().contains('without Git provenance')
    }

    def 'should refuse a dirty checkout without -force'() {
        given:
        pipeline()
        gitCheckout().close()
        Files.writeString(base.resolve('main.nf'), 'workflow { println 1 }')

        when:
        cmd().run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('-force')
        0 * client._

        when:
        def c = cmd()
        c.force = true
        c.run()

        then:
        1 * client.createPipelineRelease(_, _) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> new UploadPipelineReleaseResponse()
    }

    def 'should refuse missing or invalid name or version'() {
        given:
        pipeline(manifest)

        when:
        cmd().run()

        then:
        thrown(AbortOperationException)
        0 * client._

        where:
        manifest << ["version = '1.0.0'", "name = 'mypipe'; version = '1.0.0'", "name = 'myorg/mypipe'",
                    "name = 'MyOrg/Pipe'; version = '1.0.0'", "name = 'org/_p'; version = '1.0.0'", "name = 'org/1p'; version = '1.0.0'"]
    }

    def 'should refuse an empty bundle'() {
        given:
        pipeline()
        def repo = Git.init().setDirectory(base.toFile()).call()
        // everything untracked: nothing would be published
        repo.close()

        when:
        cmd().run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No files')
        e.message.contains('untracked')
        0 * client._
    }

    def 'should refuse a bundle over 50 MB before upload'() {
        given:
        pipeline()
        def bytes = new byte[51 * 1024 * 1024]
        new Random().nextBytes(bytes)
        Files.write(base.resolve('big.bin'), bytes)

        when:
        cmd().run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('51.0 MB')
        e.message.contains('50 MB')
        0 * client._
    }

    def 'dry run should print the release without calling the registry'() {
        given:
        pipeline()
        gitCheckout().close()
        def c = cmd()
        c.dryRun = true

        when:
        c.run()

        then:
        0 * client._
        def out = capture.toString()
        out.contains('myorg/mypipe')
        out.contains('1.2.0')
        out.contains('sha256:')
        out.contains('Size:')
        out.contains('https://github.com/myorg/mypipe')
        out.contains('v1.2.0')
    }

    def 'should report a failed upload with a retry hint'() {
        given:
        pipeline()

        when:
        cmd().run()

        then:
        1 * client.createPipelineRelease(_, _) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> { throw new RegistryException('upload failed: 500 boom') }
        def e = thrown(AbortOperationException)
        e.message.contains('500 boom')
        e.message.contains('the draft expires in 1 hour')
    }

    def 'should show the server reason on #status'() {
        given:
        pipeline()

        when:
        cmd().run()

        then:
        1 * client.createPipelineRelease(_, _) >> { throw new RegistryException("${status} ${reason}".toString()) }
        def e = thrown(AbortOperationException)
        e.message.contains(reason)
        0 * client.uploadPipelineRelease(_, _)

        where:
        status | reason
        409    | 'version 1.2.0 already exists'
        413    | 'bundle too large'
    }

    def 'should delete the temporary bundle'() {
        given:
        pipeline()
        Path sent

        when:
        cmd().run()

        then:
        1 * client.createPipelineRelease(_, _) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> { long id, Path p ->
            sent = p
            assert Files.exists(p)
            new UploadPipelineReleaseResponse()
        }
        !Files.exists(sent)
    }

    def 'should use the -registry url'() {
        given:
        pipeline()
        RegistryConfig used
        def mock = client
        def c = new CmdPipelinePublish() {
            @Override
            protected RegistryClient newClient(Map config) {
                used = registryConfig(config)
                return mock
            }
        }
        c.args = [base.toString()]
        c.registryUrl = 'https://reg.example.com/'

        when:
        c.run()

        then:
        1 * client.createPipelineRelease(_, _) >> created()
        1 * client.uploadPipelineRelease(7L, _) >> new UploadPipelineReleaseResponse()
        used.url == 'https://reg.example.com'
    }
}
