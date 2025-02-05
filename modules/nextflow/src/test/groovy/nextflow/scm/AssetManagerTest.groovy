/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.scm

import spock.lang.IgnoreIf
import spock.lang.Tag
import spock.lang.Unroll

import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Config
import org.junit.Rule
import spock.lang.Requires
import spock.lang.Specification
import test.TemporaryPath
import java.nio.file.Path
import java.nio.file.Paths
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class AssetManagerTest extends Specification {

    static String GIT_CONFIG_TEXT = '''
            [remote "origin"]
                url = https://github.com/nextflow-io/nextflow.git
                fetch = +refs/heads/*:refs/remotes/origin/*
            [branch "master"]
                remote = origin
                merge = refs/heads/master
            '''
            .stripIndent()



    static final GIT_CONFIG_LONG = '''
        [core]
            repositoryformatversion = 0
            filemode = true
            bare = false
            logallrefupdates = true
            ignorecase = true
        [remote "origin"]
            url = git@github.com:nextflow-io/nextflow.git
            fetch = +refs/heads/*:refs/remotes/origin/*
        [branch "master"]
            remote = origin
            merge = refs/heads/master
        [submodule "tests"]
            url = git@github.com:nextflow-io/tests.git
        '''
        .stripIndent()

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def setup() {
        AssetManager.root = tempDir.root.toFile()
    }

    // Helper method to grab the default brasnch if set in ~/.gitconfig
    String getLocalDefaultBranch() {
        def defaultBranch = 'master'
        def gitconfig = Paths.get(System.getProperty('user.home'),'.gitconfig');
        if(gitconfig.exists()) {
            def config = new Config()
            config.fromText(gitconfig.text)
            defaultBranch = config.getString('init', null, 'defaultBranch') ?: 'master'
        }
        return defaultBranch
    }

    @Tag("core")
    def testList() {

        given:
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        when:
        def list = AssetManager.list()
        then:
        list.sort() == ['cbcrg/pipe1','cbcrg/pipe2','ncbi/blast']

        expect:
        AssetManager.find('blast') == 'ncbi/blast'
        AssetManager.find('pipe1') == 'cbcrg/pipe1'
        AssetManager.find('pipe') as Set == ['cbcrg/pipe1', 'cbcrg/pipe2'] as Set

    }


    @Tag("core")
    @Unroll
    def 'should resolve name #input to #expected'() {
        given:
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()
        def manager = new AssetManager()

        expect:
        manager.resolveName(input) == expected

        where:
        input                   | expected
        'x/y'                  | 'x/y'
        'blast'                | 'ncbi/blast'
        'ncbi/blast/script.nf' | 'ncbi/blast'
        'blast/script.nf'      | 'ncbi/blast'
    }

    @Tag("core")
    @Unroll
    def 'should throw exception for invalid name #input'() {
        given:
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()
        def manager = new AssetManager()

        when:
        manager.resolveName(input)

        then:
        thrown(AbortOperationException)

        where:
        input << ['pipe', 'pipe/alpha/beta', '../blast/script.nf', './blast/script.nf']
    }


    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPull() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when:
        manager.download()
        then:
        noExceptionThrown()

    }


    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullTagTwice() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("v1.2")
        then:
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when:
        manager.download("v1.2")
        then:
        noExceptionThrown()
    }

    // The hashes used here are NOT associated with tags.
    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullHashTwice() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")
        then:
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when:
        def result = manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")
        then:
        noExceptionThrown()
        result == "Already-up-to-date"
    }


    // Downloading a branch first and then pulling the branch
    // should work fine, unlike with tags.
    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullBranchTwice() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("mybranch")
        then:
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when:
        manager.download("mybranch")
        then:
        noExceptionThrown()
    }

    // First clone a repo with a tag, then forget to include the -r argument
    // when you execute nextflow.
    // Note that while the download will work, execution will fail subsequently
    // at a separate check - this just tests that we don't fail because of a detached head.
    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullTagThenBranch() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("v1.2")
        then:
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when:
        manager.download()
        then:
        noExceptionThrown()
    }


    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testClone() {

        given:
        def dir = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers:[github: [auth: token]]])

        when:
        manager.clone(dir.toFile())

        then:
        dir.resolve('README.md').exists()
        dir.resolve('.git').isDirectory()

    }

    @Tag("config")
    def testGetScriptName() {

        given:
        def dir = tempDir.getRoot()
        dir.resolve('sub1').mkdir()
        dir.resolve('sub1/nextflow.config').text = "manifest.mainScript = 'pippo.nf'"
        dir.resolve('sub2').mkdir()

        when:
        def holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.getMainScriptName() == 'pippo.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub2').toFile()
        then:
        holder.getMainScriptName() == 'main.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        holder.resolveName('nextflow/hello')
        then:
        holder.getMainScriptName() == 'pippo.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        holder.resolveName('nextflow/hello/my-script.nf')
        then:
        holder.getMainScriptName() == 'my-script.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.resolveName('nextflow-io/hello/x/y/z/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'x/y/z/my-script.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.resolveName('nextflow-io/hello/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'my-script.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.resolveName('hello/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'my-script.nf'

    }

    @Tag("core")
    @Unroll
    def 'should create hub provider for #providerName'() {
        when:
        def manager = new AssetManager()
        def repo = manager.createHubProvider(providerName)

        then:
        repo instanceof RepositoryProvider

        where:
        providerName << ['github', 'bitbucket', 'gitlab']
    }


    @Tag("config")
    def 'should read manifest file' () {

        given:
        def config =
                '''
                manifest {
                    homePage = 'http://foo.com'
                    mainScript = 'hello.nf'
                    defaultBranch = 'super-stuff'
                    description = 'This pipeline do this and that'
                    author = 'Hi Dude'
                }
                '''
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = config
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')
        then:
        holder.getMainScriptName() == 'hello.nf'
        holder.manifest.getDefaultBranch() == 'super-stuff'
        holder.manifest.getHomePage() == 'http://foo.com'
        holder.manifest.getDescription() == 'This pipeline do this and that'
        holder.manifest.getAuthor() == 'Hi Dude'

    }

    @Tag("config")
    def 'should return default main script file' () {

        given:
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = 'empty: 1'
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')

        then:
        holder.getMainScriptName() == 'main.nf'
        holder.getHomePage() == 'https://github.com/foo/bar'
        holder.manifest.getDefaultBranch() == null
        holder.manifest.getDescription() == null

    }

    @Tag("git")
    @Unroll
    def 'should parse Git URL #url'() {
        given:
        def manager = new AssetManager()
        // Add provider configurations before testing
        manager.providerConfigs.add(new ProviderConfig('local-scm', [platform: 'github', server: 'http://foo.bar.com']))
        manager.providerConfigs.add(new ProviderConfig('gitea', [platform: 'gitea', server: 'http://my-server.org/sub1']))

        when:
        def result = manager.resolveNameFromGitUrl(url)

        then:
        result == expected
        manager.hub == hub

        where:
        url                                         | expected            | hub
        'nextflow/pipe'                            | null                | null
        'https://gitlab.com/pditommaso/hello.git'  | 'pditommaso/hello'  | 'gitlab'
        'file:/user/repo/projects/hello.git'       | 'local/hello'       | 'file:/user/repo/projects'
        'https://foo.bar.com/project/xyz.git'      | 'project/xyz'       | 'local-scm'
        'http://my-server.org/sub1/foo/bar.git'    | 'foo/bar'          | 'gitea'
    }

    @Tag("git")
    @Unroll
    def 'should filter remote branches #refName'() {
        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.download()
        def branches = manager.getBranchList()

        expect:
        def remote_branch = branches.find { it.name == refName }
        remote_branch != null
        AssetManager.isRemoteBranch(remote_branch) == isRemote

        where:
        refName                         | isRemote
        'refs/remotes/origin/HEAD'      | false
        'refs/remotes/origin/master'    | true
        'refs/heads/master'             | false
    }

    @Tag("git")
    @Unroll
    def "should parse Git config for #testCase"() {
        given:
        def dir = tempDir.root
        // Initialize git repo with proper configuration
        def init = Git.init()
        init.setInitialBranch('master')
        def repo = init.setDirectory(dir.toFile()).call()
        dir.resolve('.git/config').text = configText
        repo.add().addFilepattern('.').call()
        repo.commit().setSign(false).setAll(true).setMessage('First commit').call()
        repo.close()

        when:
        def manager = new AssetManager().setLocalPath(dir.toFile())

        then:
        manager.getGitConfigRemoteUrl() == expectedUrl

        where:
        testCase                        | configText                | expectedUrl
        'simple origin URL'             | GIT_CONFIG_TEXT           | 'https://github.com/nextflow-io/nextflow.git'
        'SSH origin URL'                | GIT_CONFIG_LONG           | 'git@github.com:nextflow-io/nextflow.git'
        'multiple remotes'             | '''
            [remote "origin"]
                url = https://github.com/foo/bar.git
            [remote "upstream"]
                url = https://git.example.com/group/project.git
            [branch "master"]
                remote = origin
                merge = refs/heads/master
        '''.stripIndent()               | 'https://github.com/foo/bar.git'
        'non-default branch'            | '''
            [remote "origin"]
                url = https://gitlab.com/user/repo.git
            [branch "dev"]
                remote = origin
                merge = refs/heads/dev
        '''.stripIndent()               | 'https://gitlab.com/user/repo.git'
        'custom default branch'         | '''
            [remote "origin"]
                url = https://github.com/org/repo.git
            [branch "main"]
                remote = origin
                merge = refs/heads/main
        '''.stripIndent()               | 'https://github.com/org/repo.git'
    }

    @Tag("script")
    def "should create script file object correctly"() {
        given:
        def dir = tempDir.root
        // create the repo dir
        dir.resolve('main.nf').text = "println 'Hello world'"
        dir.resolve('nextflow.config').text = 'manifest {  }'
        dir.resolve('foo.nf').text = 'this is foo content'

        // Initialize git repo with proper configuration
        def init = Git.init()
        init.setInitialBranch('master')
        def repo = init.setDirectory(dir.toFile()).call()
        repo.add().addFilepattern('.').call()
        def commit = repo.commit().setSign(false).setAll(true).setMessage('First commit').call()
        repo.close()

        when:
        def p = Mock(RepositoryProvider) { getRepositoryUrl() >> 'https://github.com/nextflow-io/nextflow' }
        and:
        def manager = new AssetManager(provider: p)
                .setLocalPath(dir.toFile())
                .setProject('nextflow-io/nextflow')
        and:
        def script = manager.getScriptFile()
        then:
        script.localPath == dir
        script.commitId == commit.name()
        script.revision == getLocalDefaultBranch()
        script.parent == dir
        script.text == "println 'Hello world'"
        script.repository == 'https://github.com/nextflow-io/nextflow'
        script.projectName == 'nextflow-io/nextflow'

        when:
        p = Mock(RepositoryProvider) { getRepositoryUrl() >> 'https://github.com/nextflow-io/nextflow' }
        and:
        manager = new AssetManager(provider: p)
                .setLocalPath(dir.toFile())
                .setProject('nextflow-io/nextflow')
        and:
        script = manager.getScriptFile('foo.nf')
        then:
        script.localPath == dir
        script.commitId == commit.name()
        script.revision == getLocalDefaultBranch()
        script.parent == dir
        script.text == "this is foo content"
        script.repository == 'https://github.com/nextflow-io/nextflow'
        script.projectName == 'nextflow-io/nextflow'
    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should download branch specified'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])

        when:
        manager.download("dev")
        then:
        folder.resolve('nextflow-io/nf-test-branch/.git').isDirectory()
        and:
        folder.resolve('nextflow-io/nf-test-branch/workflow.nf').text == "println 'Hello'\n"

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should fetch main script from branch specified'() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])

        expect:
        manager.checkValidRemoteRepo('dev')
        and:
        manager.getMainScriptName() == 'workflow.nf'

    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should download tag specified'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])

        when:
        manager.download("v0.1")
        then:
        folder.resolve('nextflow-io/nf-test-branch/.git').isDirectory()
        and:
        folder.resolve('nextflow-io/nf-test-branch/workflow.nf').text == "println 'Hello'\n"

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should identify default branch when downloading repo'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/socks', [providers: [github: [auth: token]]])

        when:
        // simulate calling `nextflow run nextflow-io/socks` without specifying a revision
        manager.download()
        manager.checkout(null)
        then:
        folder.resolve('nextflow-io/socks/.git').isDirectory()
        manager.getCurrentRevision() == 'main'

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'can filter remote branches'() {
        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.download()
        def branches = manager.getBranchList()

        when:
        def remote_head = branches.find { it.name == 'refs/remotes/origin/HEAD' }
        then:
        remote_head != null
        !AssetManager.isRemoteBranch(remote_head)

        when:
        def remote_master = branches.find { it.name == 'refs/remotes/origin/master' }
        then:
        remote_master != null
        AssetManager.isRemoteBranch(remote_master)

        when:
        def local_master = branches.find { it.name == 'refs/heads/master' }
        then:
        local_master != null
        !AssetManager.isRemoteBranch(local_master)
    }

    @Tag("config")
    def 'should handle manifest with gitmodules configuration'() {
        given:
        def config = '''
            manifest {
                gitmodules = ['module1', 'module2']
                recurseSubmodules = true
            }
            '''
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = config
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT
        dir.resolve('foo/bar/.gitmodules').text = '''
            [submodule "module1"]
                path = module1
                url = https://github.com/org/module1.git
            [submodule "module2"]
                path = module2
                url = https://github.com/org/module2.git
            '''

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')

        then:
        holder.manifest.gitmodules == ['module1', 'module2']
        holder.manifest.recurseSubmodules == true
    }

    @Tag("config")
    def 'should handle manifest with string gitmodules configuration'() {
        given:
        def config = '''
            manifest {
                gitmodules = 'module1, module2'
                recurseSubmodules = true
            }
            '''
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = config
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')

        then:
        holder.manifest.gitmodules.tokenize(', ') == ['module1', 'module2']
        holder.manifest.recurseSubmodules == true
    }

    @Tag("git")
    def 'should handle malformed manifest file'() {
        given:
        def config = '''
            manifest {
                malformed syntax here
            }
            '''
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = config
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')
        and:
        holder.getManifest()

        then:
        thrown(AbortOperationException)
    }

    @Tag("git")
    def 'should handle missing manifest file'() {
        given:
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')
        def manifest = holder.getManifest()

        then:
        manifest != null
        manifest.getMainScript() == 'main.nf'
        manifest.getDefaultBranch() == null
        manifest.getHomePage() == null
    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should handle submodule operations'() {
        given:
        def dir = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        dir.resolve('foo/bar').mkdirs()
        def mainDir = dir.resolve('foo/bar').toFile()

        // Initialize main repo with proper configuration
        def init = Git.init()
        init.setInitialBranch('master')
        def repo = init.setDirectory(mainDir).call()

        // Set up Git config with remote URL
        mainDir.resolve('.git/config').text = """[core]
        repositoryformatversion = 0
        filemode = true
        bare = false
        logallrefupdates = true
[remote "origin"]
        url = https://github.com/nextflow-io/hello.git
        fetch = +refs/heads/*:refs/remotes/origin/*
[branch "master"]
        remote = origin
        merge = refs/heads/master"""

        // Create and add a submodule
        def submoduleDir = dir.resolve('foo/submodule').toFile()
        submoduleDir.mkdirs()
        def subInit = Git.init()
        subInit.setInitialBranch('master')
        def subRepo = subInit.setDirectory(submoduleDir).call()
        new File(submoduleDir, 'test.txt').text = 'test content'
        subRepo.add().addFilepattern('.').call()
        subRepo.commit().setSign(false).setMessage('First commit in submodule').call()

        // Add submodule to main repo
        repo.submoduleAdd()
            .setPath('submodule')
            .setURI(submoduleDir.toURI().toString())
            .call()
        repo.add().addFilepattern('.').call()
        repo.commit().setSign(false).setMessage('Added submodule').call()

        // Create manifest with submodule config
        new File(mainDir, 'nextflow.config').text = '''
            manifest {
                gitmodules = ['submodule']
                recurseSubmodules = true
            }
            '''

        when:
        def holder = new AssetManager()
        holder.build('foo/bar', [providers: [github: [auth: token]]])
        holder.updateModules()

        then:
        new File(mainDir, '.gitmodules').exists()
        new File(mainDir, 'submodule/test.txt').exists()
    }

    @Tag("git")
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should handle repository errors'() {
        given:
        def dir = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        dir.resolve('foo/bar').mkdirs()
        def mainDir = dir.resolve('foo/bar').toFile()

        // Set up Git config with non-existent repository
        mainDir.resolve('.git').mkdirs()
        mainDir.resolve('.git/config').text = """[core]
        repositoryformatversion = 0
        filemode = true
        bare = false
        logallrefupdates = true
[remote "origin"]
        url = https://github.com/nextflow-io/non-existent-repo.git
        fetch = +refs/heads/*:refs/remotes/origin/*
[branch "master"]
        remote = origin
        merge = refs/heads/master"""

        when:
        def holder = new AssetManager()
        holder.build('foo/bar', [providers: [github: [auth: token]]])
        holder.getGit()

        then:
        thrown(RepositoryNotFoundException)
    }

}
