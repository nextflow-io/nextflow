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
import nextflow.cli.CmdDrop
import nextflow.cli.CmdInfo
import nextflow.cli.CmdList
import nextflow.cli.CmdPipeline
import nextflow.cli.CmdPull
import nextflow.cli.CmdView
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.scm.MultiRevisionRepositoryStrategy
import org.eclipse.jgit.api.Git
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

/**
 * Tests for the `-global` option of the pipeline commands, which use the user-level store,
 * and for the deprecated project commands
 */
class CmdPipelineUserLevelTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path base

    File savedRoot

    def setup() {
        savedRoot = AssetManager.root
        AssetManager.root = Files.createDirectories(base.resolve('assets')).toFile()
    }

    def cleanup() {
        Plugins.stop()
        AssetManager.root = savedRoot
    }

    private static Path assets() {
        return AssetManager.root.toPath()
    }

    private static void seedDirs(String... names) {
        names.each { Files.createDirectories(assets().resolve(it)) }
    }

    /** A Git repository with one commit and the given remote, as downloaded by the legacy strategy */
    private static Path seedRepo(String name, String remote) {
        final dir = Files.createDirectories(assets().resolve(name))
        Files.writeString(dir.resolve('main.nf'), 'workflow { }')
        try (final git = Git.init().setDirectory(dir.toFile()).setInitialBranch('master').call()) {
            git.add().addFilepattern('main.nf').call()
            git.commit().setMessage('first').setSign(false).call()
            final config = git.repository.config
            config.setString('remote', 'origin', 'url', remote)
            config.save()
        }
        return dir
    }

    // the printed lines, without the log
    private List<String> stdout() {
        return capture.toString().readLines().findAll { String it -> it && !(it =~ /^\d\d:\d\d:\d\d\.\d+ \[/) }
    }

    private static <T> T parse(String... args) {
        return new Launcher().parseMainArgs(args).command as T
    }

    def 'should list the pipelines in the user-level store'() {
        given:
        seedDirs('nextflow-io/hello', 'acme/foo')

        when:
        new CmdPipelineList(global: true).run()
        then:
        stdout() == ['acme/foo', 'nextflow-io/hello']
    }

    def 'should list the pipelines in the user-level store as json'() {
        given:
        seedDirs('nextflow-io/hello', 'acme/foo')

        when:
        new CmdPipelineList(global: true, output: 'json').run()
        then:
        (new JsonSlurper().parseText(stdout().join('\n')) as Map).pipelines == [[name: 'acme/foo'], [name: 'nextflow-io/hello']]
    }

    def 'should report an empty user-level store'() {
        when:
        new CmdPipelineList(global: true).run()
        then:
        stdout() == ['No pipelines installed in the user-level store']
    }

    def 'should hint at the user-level store when no pipeline is installed in the project'() {
        given:
        seedDirs(*USER_LEVEL)

        when:
        new CmdPipelineList(root: base.resolve('project')).run()
        then:
        stdout() == EXPECTED

        where:
        USER_LEVEL                          | EXPECTED
        []                                  | ['No pipelines installed']
        ['acme/foo']                        | ['No pipelines installed', '1 pipeline is installed in the user-level store -- use `nextflow pipeline list -g` to list them']
        ['acme/foo', 'nextflow-io/hello']   | ['No pipelines installed', '2 pipelines are installed in the user-level store -- use `nextflow pipeline list -g` to list them']
    }

    def 'should view a pipeline in the user-level store'() {
        given:
        final dir = seedRepo('acme/hello', 'https://github.com/acme/hello.git')

        when:
        new CmdPipelineView(args: ['acme/hello'], global: true).run()
        then:
        final out = stdout()
        out.contains(' project name: acme/hello')
        out.contains(" local path  : ${dir}".toString())
        out.contains(' main script : main.nf')
    }

    def 'should view a pipeline in the user-level store as json'() {
        given:
        final dir = seedRepo('acme/hello', 'https://github.com/acme/hello.git')

        when:
        new CmdPipelineView(args: ['acme/hello'], global: true, output: 'json').run()
        final json = new JsonSlurper().parseText(stdout().join('\n')) as Map
        then:
        json.name == 'acme/hello'
        json.repositoryUrl == 'https://github.com/acme/hello'
        json.localPath == dir.toString()
        json.mainScript == 'main.nf'
        json.revisions instanceof Map
    }

    def 'should reject an unknown pipeline or a version when viewing the user-level store'() {
        when:
        new CmdPipelineView(args: ['acme/hello'], global: true).run()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Pipeline acme/hello is not installed in the user-level store'

        when:
        new CmdPipelineView(args: ['acme/hello'], global: true, version: '1.0').run()
        then:
        e = thrown(AbortOperationException)
        e.message == 'Option `-version` is not supported with `-global`'
    }

    def 'should map a registry page URL to the registry name in the user-level store'() {
        given:
        final launcher = new Launcher().parseMainArgs('pipeline')
        final cmd = new CmdPipelineView(launcher: launcher)

        expect:
        PipelineCmdSupport.userLevelName(cmd, 'https://registry.nextflow.io/pipelines/nf-core/demo') == 'nf-core/demo'
        PipelineCmdSupport.userLevelName(cmd, 'nf-core/demo') == 'nf-core/demo'
        PipelineCmdSupport.userLevelName(cmd, 'https://github.com/nf-core/demo') == 'https://github.com/nf-core/demo'
    }

    def 'should remove a pipeline from the user-level store'() {
        given:
        final dir = seedRepo('acme/hello', 'https://github.com/acme/hello.git')

        when:
        new CmdPipelineRemove(args: ['acme/hello'], global: true).run()
        then:
        !Files.exists(dir)
        stdout() == ['Pipeline acme/hello removed']
    }

    def 'should report a pipeline that is not in the user-level store'() {
        when:
        new CmdPipelineRemove(args: ['acme/hello'], global: true).run()
        then:
        stdout() == ['Pipeline acme/hello is not installed in the user-level store']
    }

    def 'should find the provider of a pipeline whose default branch is not checked out'() {
        given: 'a bare repository from GitLab, as a registry release leaves it, with only its pinned commit checked out'
        final source = base.resolve('source')
        Files.createDirectories(source)
        Files.writeString(source.resolve('main.nf'), 'workflow { }')
        try (final git = Git.init().setDirectory(source.toFile()).setInitialBranch('master').call()) {
            git.add().addFilepattern('main.nf').call()
            git.commit().setMessage('first').setSign(false).call()
        }
        final bare = assets().resolve("${MultiRevisionRepositoryStrategy.REPOS_SUBDIR}/acme/hello/${MultiRevisionRepositoryStrategy.BARE_REPO}")
        Git.cloneRepository().setURI(source.toUri().toString()).setDirectory(bare.toFile()).setBare(true).call().close()
        try (final git = Git.open(bare.toFile())) {
            final config = git.repository.config
            config.setString('remote', 'origin', 'url', 'https://gitlab.com/acme/hello.git')
            config.save()
        }

        when:
        final manager = new AssetManager('acme/hello')
        then:
        !manager.isLocal()
        manager.hub == 'gitlab'

        cleanup:
        manager?.close()
    }

    def 'should hide the deprecated project commands'() {
        when:
        new Launcher().printCommands([new CmdDrop(), new CmdInfo(), new CmdList(), new CmdPipeline(), new CmdPull(), new CmdView()])
        then:
        capture.toString() == '''
                Commands:
                  info       Print system runtime information
                  pipeline   Manage registry pipelines
                  pull       Download or update a project

                '''
                .stripIndent()
    }

    def 'should point to the replacement in the help of a deprecated command'() {
        when:
        new Launcher().usage(NAME)
        then:
        capture.toString().contains(DESCRIPTION)

        where:
        NAME    | DESCRIPTION
        'drop'  | '(deprecated -- use `nextflow pipeline remove -g`)'
        'list'  | '(deprecated -- use `nextflow pipeline list -g`)'
        'view'  | '(deprecated)'
    }

    def 'should not show the project options in the help of info'() {
        when:
        new Launcher().usage('info')
        then:
        final out = capture.toString()
        out.contains('-d')
        !out.contains('project name')
        !out.contains('-check-updates')
        !out.contains('yaml')
    }

    def 'should warn and still run the deprecated project commands'() {
        given:
        seedRepo('acme/hello', 'https://github.com/acme/hello.git')

        when:
        parse(*ARGS).run()
        then:
        capture.toString().contains(WARNING)

        where:
        ARGS                                | WARNING
        ['list']                            | 'The `list` command is deprecated -- use `nextflow pipeline list -g` instead'
        ['info', 'acme/hello']              | 'The `info <project>` command is deprecated -- use `nextflow pipeline view -g <project>` instead'
        ['view', 'acme/hello']              | 'The `view` command is deprecated -- use `nextflow pipeline view -g` to find the local path of a pipeline'
        ['drop', 'acme/hello']              | 'The `drop` command is deprecated -- use `nextflow pipeline remove -g` instead'
    }
}
