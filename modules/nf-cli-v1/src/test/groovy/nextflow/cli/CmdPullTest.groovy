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

import static nextflow.scm.MultiRevisionRepositoryStrategy.REPOS_SUBDIR
import static nextflow.scm.MultiRevisionRepositoryStrategy.REVISION_SUBDIR

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.client.RegistryClient
import nextflow.plugin.Plugins
import org.eclipse.jgit.api.Git
import spock.lang.IgnoreIf

import java.nio.file.Files

import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class CmdPullTest extends Specification {

    def cleanup() {
        Plugins.stop()
    }

    @Requires({ System.getenv('NXF_GITHUB_ACCESS_TOKEN') })
    def 'should pull the github repository in the local folder'() {

        given:
        def accessToken = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def dir = Files.createTempDirectory('test')
        def cmd = new CmdPull(args: ['nextflow-io/hello'], root: dir.toFile(), revision: '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8', hubUser: accessToken)

        when:
        cmd.run()
        then:
        dir.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/' + '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' + '/.git').exists()
        dir.resolve(REPOS_SUBDIR+ '/nextflow-io/hello/' + REVISION_SUBDIR + '/' + '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' + '/README.md').exists()

        cleanup:
        dir?.deleteDir()
        Plugins.stop()
    }

    def 'should pull a registry release at its pinned commit under its registry name'() {
        given: 'a bare repo with two commits'
        def folder = Files.createTempDirectory('test')
        def sourceDir = folder.resolve('source-repo')
        Files.createDirectories(sourceDir)
        def git = Git.init().setDirectory(sourceDir.toFile()).call()
        sourceDir.resolve('main.nf').text = "println 'v1'"
        git.add().addFilepattern('.').call()
        def v1 = git.commit().setSign(false).setMessage('v1').call().name()
        sourceDir.resolve('main.nf').text = "println 'v2'"
        git.add().addFilepattern('.').call()
        git.commit().setSign(false).setMessage('v2').call()
        git.close()
        def bareDir = folder.resolve('hello.git')
        Git.cloneRepository().setURI(sourceDir.toUri().toString()).setDirectory(bareDir.toFile()).setBare(true).call().close()

        and: 'a registry release 1.0.0 pinned to the first commit'
        def client = Mock(RegistryClient)
        def pipeline = new Pipeline().repositoryUrl("file:${bareDir}".toString()).latest(new PipelineRelease().version('2.0.0'))
        client.getPipeline('acme/hello') >> pipeline
        client.getPipelineRelease('acme/hello', '1.0.0') >> new PipelineRelease().version('1.0.0').gitRevision('v1.0.0').commitId(v1)

        and:
        def root = folder.resolve('assets')
        def cmd = new CmdPull(args: ['acme/hello'], root: root.toFile(), revision: '1.0.0', hubProvider: 'registry.nextflow.io', registryClient: client)
        cmd.launcher = Mock(Launcher) { getOptions() >> null }

        when:
        cmd.run()

        then:
        root.resolve("$REPOS_SUBDIR/acme/hello/$REVISION_SUBDIR/$v1/main.nf").text == "println 'v1'"

        cleanup:
        folder?.deleteDir()
    }
}
