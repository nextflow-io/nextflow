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

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import groovy.json.JsonSlurper
import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.cli.CmdBase
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

    private int mark

    // only check the output printed from now on
    private void resetOutput() {
        mark = capture.toString().length()
    }

    // the printed lines, without the log; the jgit log may end a message with a lone quote on its own line
    private List<String> stdout() {
        return capture.toString().substring(mark).readLines().findAll { String it -> it && it != "'" && !(it =~ /^\d\d:\d\d:\d\d\.\d+ \[/) }
    }

    // each command initializes the plugins, as on the command line
    private static void exec(CmdBase cmd) {
        Plugins.stop()
        cmd.run()
    }

    private static <T> T parse(String... args) {
        return new Launcher().parseMainArgs(args).command as T
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

    def 'should not accept -g when viewing a pipeline'() {
        when:
        new JCommander(new CmdPipelineView()).parse('-g', 'acme/hello')
        then:
        thrown(ParameterException)
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
                  pipeline   Manage pipelines

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
        'pull'  | '(deprecated -- use `nextflow pipeline install -g`)'
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
        ['info', 'acme/hello']              | 'The `info <project>` command is deprecated -- use `nextflow pipeline list -g` instead'
        ['view', 'acme/hello']              | 'The `view` command is deprecated -- use `nextflow pipeline list -g` to find the local path of a pipeline'
        ['drop', 'acme/hello']              | 'The `drop` command is deprecated -- use `nextflow pipeline remove -g` instead'
    }

    def 'should warn and still run the deprecated pull command'() {
        given:
        final repo = gitRepo('hello')

        when:
        parse('pull', repo.url).run()
        then:
        capture.toString().contains('The `pull` command is deprecated -- use `nextflow pipeline install -g` instead')
        Files.exists(assets().resolve("${MultiRevisionRepositoryStrategy.REPOS_SUBDIR}/local/hello/${MultiRevisionRepositoryStrategy.REVISION_SUBDIR}/${repo.v2}/main.nf"))
    }

    /** A Git repository with two commits, with lightweight tags v1.0.0 and v2.0.0, on the master branch; a `file:` URL is stored under `local/<name>` */
    private Map<String,String> gitRepo(String name) {
        final dir = Files.createDirectories(base.resolve('remote').resolve(name))
        try (final git = Git.init().setDirectory(dir.toFile()).setInitialBranch('master').call()) {
            Files.writeString(dir.resolve('main.nf'), "println 'v1'")
            git.add().addFilepattern('main.nf').call()
            final v1 = git.commit().setMessage('v1').setSign(false).call().name()
            git.tag().setName('v1.0.0').setAnnotated(false).call()
            Files.writeString(dir.resolve('main.nf'), "println 'v2'")
            git.add().addFilepattern('main.nf').call()
            final v2 = git.commit().setMessage('v2').setSign(false).call().name()
            git.tag().setName('v2.0.0').setAnnotated(false).call()
            return [url: dir.toUri().toString(), v1: v1, v2: v2]
        }
    }

    private static Path clone(String name, String commit) {
        return assets().resolve("${MultiRevisionRepositoryStrategy.REPOS_SUBDIR}/${name}/${MultiRevisionRepositoryStrategy.REVISION_SUBDIR}/${commit}")
    }

    private CmdPipelineInstall installCmd(Map opts) {
        // each install initializes the plugins, as on the command line
        Plugins.stop()
        final cmd = new CmdPipelineInstall(opts + [global: true])
        cmd.launcher = new Launcher().parseMainArgs('pipeline')
        return cmd
    }

    private RegistryClient registryOf(Map<String,String> repo) {
        final client = Mock(RegistryClient)
        client.getPipeline('acme/hello') >> new Pipeline().repositoryUrl(repo.url).latest(new PipelineRelease().version('2.0.0').gitRevision('v2.0.0').commitId(repo.v2))
        client.getPipelineRelease('acme/hello', '1.0.0') >> new PipelineRelease().version('1.0.0').gitRevision('v1.0.0').commitId(repo.v1)
        return client
    }

    def 'should install a registry release in the user-level store at its pinned commit'() {
        given:
        final repo = gitRepo('hello')
        final client = registryOf(repo)

        when: 'the latest release'
        installCmd(args: [NAME], client: client).run()
        then:
        clone('acme/hello', repo.v2).resolve('main.nf').text == "println 'v2'"
        stdout() == ["Pipeline acme/hello@2.0.0 installed in ${clone('acme/hello', repo.v2)}".toString()]

        when: 'again'
        resetOutput()
        installCmd(args: [NAME], client: client).run()
        then:
        stdout() == ['Pipeline acme/hello@2.0.0 is already installed']

        when: 'an older release, kept next to the latest one'
        resetOutput()
        installCmd(args: [NAME], version: '1.0.0', client: client).run()
        then:
        clone('acme/hello', repo.v1).resolve('main.nf').text == "println 'v1'"
        Files.exists(clone('acme/hello', repo.v2))
        stdout() == ["Pipeline acme/hello@1.0.0 installed in ${clone('acme/hello', repo.v1)}".toString()]

        where:
        NAME << ['acme/hello', 'https://registry.nextflow.io/pipelines/acme/hello']
    }

    def 'should install a Git repository in the user-level store'() {
        given:
        final repo = gitRepo('hello')
        final name = 'local/hello'

        when: 'the default branch'
        installCmd(args: [repo.url]).run()
        then:
        clone(name, repo.v2).resolve('main.nf').text == "println 'v2'"
        stdout() == ["Pipeline ${name}@master installed in ${clone(name, repo.v2)}".toString()]

        when: 'again'
        resetOutput()
        installCmd(args: [repo.url]).run()
        then:
        stdout() == ["Pipeline ${name}@master is already installed".toString()]

        when: 'a tag'
        resetOutput()
        installCmd(args: [repo.url], version: 'v1.0.0').run()
        then:
        clone(name, repo.v1).resolve('main.nf').text == "println 'v1'"
        stdout() == ["Pipeline ${name}@v1.0.0 installed in ${clone(name, repo.v1)}".toString()]
    }

    def 'should fetch a branch again when installing it in the user-level store'() {
        given:
        final repo = gitRepo('hello')
        final name = 'local/hello'
        installCmd(args: [repo.url]).run()

        and: 'a new commit on the branch'
        final v3 = Git.open(base.resolve('remote/hello').toFile()).withCloseable { git ->
            Files.writeString(base.resolve('remote/hello/main.nf'), "println 'v3'")
            git.add().addFilepattern('main.nf').call()
            git.commit().setMessage('v3').setSign(false).call().name()
        }

        when:
        resetOutput()
        installCmd(args: [repo.url], version: 'master').run()
        then:
        clone(name, v3).resolve('main.nf').text == "println 'v3'"
        Files.exists(clone(name, repo.v2))
        stdout() == ["Pipeline ${name}@master installed in ${clone(name, v3)}".toString()]
    }

    def 'should tell how to install a Git pipeline when a name is not in the registry'() {
        given:
        final client = Mock(RegistryClient)
        client.getPipeline('nextflow-io/hello') >> { throw new RegistryException('Not found', 404) }
        client.getPipelineRelease('nextflow-io/hello', _) >> { throw new RegistryException('Not found', 404) }

        when:
        installCmd(args: ['nextflow-io/hello'], version: VERSION, client: client).run()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Pipeline nextflow-io/hello is not in the registry -- a Git pipeline is installed by its repository URL, e.g. `nextflow pipeline install -g https://github.com/nextflow-io/hello`'

        where:
        VERSION << [null, '1.0.0']
    }

    def 'should report a missing version of a registry pipeline'() {
        given:
        final client = registryOf(gitRepo('hello'))
        client.getPipelineRelease('acme/hello', '9.9.9') >> { throw new RegistryException('Not found', 404) }

        when:
        installCmd(args: ['acme/hello'], version: '9.9.9', client: client).run()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Pipeline acme/hello@9.9.9 not found in the registry'
    }

    def 'should reject -force when installing in the user-level store'() {
        when:
        installCmd(args: ['acme/hello'], force: true).run()
        then:
        def e = thrown(AbortOperationException)
        e.message.startsWith('Option `-force` is not supported with `-global`')
    }

    def 'should accept -r as a hidden alias of -version'() {
        given:
        final client = registryOf(gitRepo('hello'))

        when:
        final cmd = installCmd(client: client)
        new JCommander(cmd).parse('-g', '-r', '1.0.0', 'acme/hello')
        cmd.run()
        then:
        cmd.version == '1.0.0'
        stdout().last().startsWith('Pipeline acme/hello@1.0.0 installed in ')

        when:
        final view = new CmdPipelineView(client: client)
        view.launcher = new Launcher().parseMainArgs('pipeline')
        new JCommander(view).parse('-r', '1.0.0', 'acme/hello')
        then:
        view.version == null
        view.revision == '1.0.0'

        expect:
        !usage(new CmdPipelineInstall()).contains('-r')
        usage(new CmdPipelineInstall()).contains('-version')
        !usage(new CmdPipelineView()).contains('-r')
    }

    /**
     * A repository with lightweight tags v1.0.0 and v2.0.0 on master, and a `dev` branch from v2.0.0
     * with an annotated tag v3.0.0 and an untagged head
     */
    private Map<String,String> branchedRepo() {
        final repo = gitRepo('hello')
        final dir = base.resolve('remote/hello')
        try (final git = Git.open(dir.toFile())) {
            git.checkout().setCreateBranch(true).setName('dev').call()
            Files.writeString(dir.resolve('main.nf'), "println 'v3'")
            git.add().addFilepattern('main.nf').call()
            repo.v3 = git.commit().setMessage('v3').setSign(false).call().name()
            repo.v3tag = git.tag().setName('v3.0.0').setAnnotated(true).setMessage('v3').setSigned(false).call().objectId.name()
            Files.writeString(dir.resolve('main.nf'), "println 'v4'")
            git.add().addFilepattern('main.nf').call()
            repo.v4 = git.commit().setMessage('v4').setSign(false).call().name()
            git.checkout().setName('master').call()
        }
        return repo
    }

    /** Install the default branch, a tag, an annotated tag and a branch of the {@link #branchedRepo} */
    private Map<String,String> installRevisions() {
        final repo = branchedRepo()
        installCmd(args: [repo.url]).run()
        installCmd(args: [repo.url], version: 'v1.0.0').run()
        installCmd(args: [repo.url], version: 'v3.0.0').run()
        installCmd(args: [repo.url], version: 'dev').run()
        resetOutput()
        return repo
    }

    private static String line(String name, String version, String commit, String status) {
        return name.padRight(13) + version.padRight(9) + commit.padRight(12) + status
    }

    def 'should list each downloaded revision in the user-level store'() {
        given:
        final repo = installRevisions()

        when:
        exec(new CmdPipelineList(global: true))
        then: 'a tag is preferred to a branch, and an annotated tag shows its commit'
        stdout() == [
            line('NAME', 'VERSION', 'COMMIT', 'STATUS'),
            line('local/hello', 'dev', repo.v4.take(10), 'ok'),
            line('local/hello', 'v1.0.0', repo.v1.take(10), 'ok'),
            line('local/hello', 'v2.0.0', repo.v2.take(10), 'ok'),
            line('local/hello', 'v3.0.0', repo.v3.take(10), 'ok'),
        ]

        when: 'a clone is changed'
        resetOutput()
        Files.writeString(clone('local/hello', repo.v1).resolve('main.nf'), 'changed')
        exec(new CmdPipelineList(global: true))
        then:
        stdout()[2] == line('local/hello', 'v1.0.0', repo.v1.take(10), 'modified')
    }

    def 'should list each downloaded revision in the user-level store as json'() {
        given:
        final repo = installRevisions()
        Files.writeString(clone('local/hello', repo.v4).resolve('main.nf'), 'changed')

        when:
        exec(new CmdPipelineList(global: true, output: 'json'))
        final json = new JsonSlurper().parseText(stdout().join('\n')) as Map
        final url = new AssetManager('local/hello').withCloseable { it.repositoryUrl }
        then:
        json.pipelines == [
            [name: 'local/hello', version: 'dev', commit: repo.v4, status: 'modified', path: clone('local/hello', repo.v4).toString(), repositoryUrl: url],
            [name: 'local/hello', version: 'v1.0.0', commit: repo.v1, status: 'ok', path: clone('local/hello', repo.v1).toString(), repositoryUrl: url],
            [name: 'local/hello', version: 'v2.0.0', commit: repo.v2, status: 'ok', path: clone('local/hello', repo.v2).toString(), repositoryUrl: url],
            [name: 'local/hello', version: 'v3.0.0', commit: repo.v3, status: 'ok', path: clone('local/hello', repo.v3tag).toString(), repositoryUrl: url],
        ]
    }

    def 'should list a legacy download in the user-level store, skipping a broken one'() {
        given:
        final dir = seedRepo('acme/hello', 'https://github.com/acme/hello.git')
        final commit = Git.open(dir.toFile()).withCloseable { it.repository.resolve('HEAD').name() }
        and: 'a multi-revision download whose bare repository is broken'
        Files.createDirectories(assets().resolve("${MultiRevisionRepositoryStrategy.REPOS_SUBDIR}/acme/broken/${MultiRevisionRepositoryStrategy.BARE_REPO}"))
        Files.createDirectories(clone('acme/broken', 'a' * 40))

        when:
        exec(new CmdPipelineList(global: true))
        then:
        stdout() == [
            'NAME        VERSION  COMMIT      STATUS',
            "acme/hello  master   ${commit.take(10)}  ok".toString()
        ]
        capture.toString().contains('Unable to read pipeline acme/broken in the user-level store')

        when:
        resetOutput()
        exec(new CmdPipelineList(global: true, output: 'json'))
        then:
        (new JsonSlurper().parseText(stdout().join('\n')) as Map).pipelines == [
            [name: 'acme/hello', version: 'master', commit: commit, status: 'ok', path: dir.toString(), repositoryUrl: 'https://github.com/acme/hello']
        ]
    }

    def 'should remove one revision from the user-level store'() {
        given:
        final repo = installRevisions()

        when: 'by tag'
        exec(new CmdPipelineRemove(args: ['local/hello'], global: true, version: 'v1.0.0'))
        then:
        !Files.exists(clone('local/hello', repo.v1))
        Files.exists(clone('local/hello', repo.v2))
        stdout() == ['Pipeline local/hello@v1.0.0 removed']

        when: 'by annotated tag'
        resetOutput()
        exec(new CmdPipelineRemove(args: ['local/hello'], global: true, version: 'v3.0.0'))
        then:
        !Files.exists(clone('local/hello', repo.v3tag))
        stdout() == ['Pipeline local/hello@v3.0.0 removed']

        when: 'by short commit, with the -r alias'
        resetOutput()
        final cmd = new CmdPipelineRemove()
        new JCommander(cmd).parse('-g', '-r', repo.v4.take(7), 'local/hello')
        exec(cmd)
        then:
        !Files.exists(clone('local/hello', repo.v4))
        stdout() == ["Pipeline local/hello@${repo.v4.take(7)} removed".toString()]

        when: 'the last revision, by branch'
        resetOutput()
        exec(new CmdPipelineRemove(args: ['local/hello'], global: true, version: 'master'))
        then: 'removes the pipeline'
        !Files.exists(assets().resolve("${MultiRevisionRepositoryStrategy.REPOS_SUBDIR}/local/hello"))
        stdout() == ['Pipeline local/hello@master removed']
    }

    def 'should not remove an unknown or modified revision from the user-level store'() {
        given:
        final repo = installRevisions()

        when:
        exec(new CmdPipelineRemove(args: ['local/hello'], global: true, version: '9.9.9'))
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Pipeline local/hello@9.9.9 is not installed in the user-level store -- installed revisions: dev, v1.0.0, v2.0.0, v3.0.0'

        when:
        Files.writeString(clone('local/hello', repo.v1).resolve('main.nf'), 'changed')
        exec(new CmdPipelineRemove(args: ['local/hello'], global: true, version: 'v1.0.0'))
        then:
        e = thrown(AbortOperationException)
        e.message == 'Pipeline local/hello@v1.0.0 has local modifications -- use -force to remove it'
        Files.exists(clone('local/hello', repo.v1))

        when:
        exec(new CmdPipelineRemove(args: ['local/hello'], global: true, version: 'v1.0.0', force: true))
        then:
        !Files.exists(clone('local/hello', repo.v1))
        Files.exists(clone('local/hello', repo.v2))
    }

    def 'should remove a legacy download by its version'() {
        given:
        final dir = seedRepo('acme/hello', 'https://github.com/acme/hello.git')

        when:
        exec(new CmdPipelineRemove(args: ['acme/hello'], global: true, version: 'master'))
        then:
        !Files.exists(dir)
        stdout() == ['Pipeline acme/hello@master removed']
    }

    def 'should reject -version when removing from the project'() {
        when:
        exec(new CmdPipelineRemove(args: ['acme/hello'], version: '1.0.0', root: base.resolve('project')))
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Option `-version` requires `-global` -- the project-level store holds a single version of each pipeline'
    }

    def 'should hint at the user-level store when removing a pipeline that is not in the project'() {
        given:
        seedRepo('acme/hello', 'https://github.com/acme/hello.git')

        when:
        exec(new CmdPipelineRemove(args: [NAME], root: base.resolve('project')))
        then:
        stdout() == EXPECTED

        where:
        NAME            | EXPECTED
        'acme/hello'    | ['Pipeline acme/hello is not installed', 'Pipeline acme/hello is installed in the user-level store -- use `nextflow pipeline remove -g acme/hello` to remove it']
        'acme/other'    | ['Pipeline acme/other is not installed']
    }

    private static String usage(Object cmd) {
        final out = new StringBuilder()
        new JCommander(cmd).usage(out)
        return out.toString()
    }
}
