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

import java.nio.file.Files
import java.nio.file.Path

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.client.RegistryClient
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.transport.URIish
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Project commands given a registry address of a pipeline
 */
class RegistryPipelineCmdTest extends Specification {

    static final String PAGE_URL = 'https://registry.nextflow.io/pipelines/acme/hello'

    @TempDir
    Path folder

    private File savedRoot
    private Path root
    private String bareUrl
    private String v1
    private RegistryClient client

    def setup() {
        savedRoot = AssetManager.root
        root = folder.resolve('assets')
        AssetManager.root = root.toFile()

        // a bare repo with two commits; release 1.0.0 is pinned to the first
        def source = folder.resolve('source')
        def git = Git.init().setDirectory(source.toFile()).call()
        source.resolve('main.nf').text = "println 'v1'"
        git.add().addFilepattern('.').call()
        v1 = git.commit().setSign(false).setMessage('v1').call().name()
        source.resolve('main.nf').text = "println 'v2'"
        git.add().addFilepattern('.').call()
        git.commit().setSign(false).setMessage('v2').call()
        git.close()
        def bare = folder.resolve('hello.git')
        Git.cloneRepository().setURI(source.toUri().toString()).setDirectory(bare.toFile()).setBare(true).call().close()
        bareUrl = "file:${bare}".toString()

        client = Mock(RegistryClient)
        client.getPipeline('acme/hello') >> new Pipeline().repositoryUrl(bareUrl).latest(new PipelineRelease().version('2.0.0'))
        client.getPipelineRelease('acme/hello', '1.0.0') >> new PipelineRelease().version('1.0.0').gitRevision('v1.0.0').commitId(v1)
    }

    def cleanup() {
        AssetManager.root = savedRoot
        Plugins.stop()
    }

    private void pullRelease() {
        def cmd = new CmdPull(args: [PAGE_URL], root: root.toFile(), revision: '1.0.0', registryClient: client)
        cmd.launcher = launcher()
        cmd.run()
        Plugins.stop()
    }

    private Launcher launcher() {
        Mock(Launcher) { getOptions() >> null }
    }

    private static String stdout(Closure action) {
        final saved = System.out
        final buffer = new ByteArrayOutputStream()
        System.out = new PrintStream(buffer)
        try {
            action.call()
        }
        finally {
            System.out = saved
        }
        return buffer.toString()
    }

    def 'should give the info of a cached registry pipeline by its page URL' () {
        given:
        pullRelease()
        def buffer = new ByteArrayOutputStream()
        def cmd = new CmdInfo(args: [PAGE_URL], out: new PrintStream(buffer))
        cmd.launcher = launcher()

        when:
        cmd.run()

        then:
        buffer.toString().contains(' project name: acme/hello')
        buffer.toString().contains(" local path  : $root/$REPOS_SUBDIR/acme/hello")
    }

    def 'should view a cached registry pipeline by its page URL' () {
        given:
        pullRelease()
        def cmd = new CmdView(args: [PAGE_URL], revision: v1, quiet: true)
        cmd.launcher = launcher()

        when:
        def lines = stdout { cmd.run() }.readLines()

        then:
        lines.contains("println 'v1'")
        !lines.contains("println 'v2'")
    }

    def 'should drop a cached registry pipeline by its page URL' () {
        given:
        pullRelease()
        def cmd = new CmdDrop(args: [PAGE_URL], force: true)
        cmd.launcher = launcher()

        when:
        cmd.run()

        then:
        !AssetManager.list().contains('acme/hello')
    }

    def 'should not read the registry config for a plain name' () {
        given:
        pullRelease()

        when: 'there is no launcher to read the config with'
        new CmdInfo(args: ['acme/other']).run()

        then:
        def e = thrown(AbortOperationException)
        e.message == 'Unknown project `acme/other`'
    }

    def 'should fail clearly for an uncached registry pipeline when NXF_HUB is a registry host' () {
        given:
        SysEnv.push(NXF_HUB: 'registry.nextflow.io')
        CMD.launcher = launcher()

        when:
        CMD.run()

        then:
        def e = thrown(AbortOperationException)
        e.message == 'Pipeline `acme/other` is not downloaded -- use `nextflow pull` to download it'

        cleanup:
        SysEnv.pop()

        where:
        CMD << [new CmdInfo(args: ['acme/other']), new CmdView(args: ['acme/other']), new CmdDrop(args: ['acme/other'])]
    }

    def 'should list a cached registry pipeline under its registry name' () {
        given:
        pullRelease()

        expect:
        stdout { new CmdList().run() }.readLines().findAll { !it.contains(' DEBUG ') } == ['acme/hello']
    }

    def 'should clone a registry release at its pinned commit' () {
        given:
        def target = folder.resolve('clone')
        def cmd = new CmdClone(args: [PAGE_URL, target.toString()], revision: '1.0.0', registryClient: client)
        cmd.launcher = launcher()

        when:
        cmd.run()

        then:
        target.resolve('main.nf').text == "println 'v1'"
        def git = Git.open(target.toFile())
        git.repository.resolve('HEAD').name() == v1
        new URIish(git.repository.config.getString('remote', 'origin', 'url')) == new URIish(bareUrl)
        and: 'a full clone, as for a Git project'
        git.repository.refDatabase.getRefsByPrefix('refs/remotes/origin/')

        cleanup:
        git?.close()
    }
}
