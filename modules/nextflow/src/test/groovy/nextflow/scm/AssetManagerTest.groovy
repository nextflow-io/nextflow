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
 * Tests for the AssetManager class which handles operations on remote and local installed pipelines.
 * Tests are organized into logical groups:
 * - Project name resolution and validation
 * - SCM provider handling
 * - Git URL parsing and repository configuration
 * - Script and manifest handling
 * - Branch and reference management
 * - Integration tests with real repositories
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

    // Helper method to grab the default branch if set in ~/.gitconfig
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

    // ----------------------------------------------------------------
    // Project name resolution and validation tests
    // ----------------------------------------------------------------

    @Unroll
    def "resolving project name '#name' should return '#expected'"() {
        given: 'an asset manager instance'
        def manager = new AssetManager()

        expect: 'the name is correctly resolved'
        manager.resolveName(name) == expected

        where:
        name                                        | expected
        'nextflow-io/hello'                        | 'nextflow-io/hello'
        'hello'                                    | 'nextflow-io/hello'
        'https://github.com/nextflow-io/hello.git' | 'nextflow-io/hello'
        'https://gitlab.com/user/repo.git'         | 'user/repo'
        'file:///path/to/repo.git'                 | 'local/repo'
    }

    @Unroll
    def "invalid project name '#name' should throw AbortOperationException"() {
        given: 'an asset manager instance'
        def manager = new AssetManager()

        when: 'resolving an invalid project name'
        manager.resolveName(name)

        then: 'an exception is thrown'
        thrown(AbortOperationException)

        where: 'invalid names include'
        name << ['./invalid', '../invalid', '/invalid', 'a/b/c/d']
    }

    // ----------------------------------------------------------------
    // SCM provider handling tests
    // ----------------------------------------------------------------

    @Unroll
    def "creating hub provider for '#providerName' should return #expectedClass.simpleName"() {
        given: 'an asset manager instance'
        def manager = new AssetManager()

        when: 'creating a provider'
        def provider = manager.createHubProvider(providerName)

        then: 'the correct provider type is returned'
        provider.class == expectedClass

        where:
        providerName  | expectedClass
        'github'      | GithubRepositoryProvider
        'gitlab'      | GitlabRepositoryProvider
        'bitbucket'   | BitbucketRepositoryProvider
        'gitea'       | GiteaRepositoryProvider
    }

    // ----------------------------------------------------------------
    // Git URL parsing and repository configuration tests
    // ----------------------------------------------------------------

    @Unroll
    def "parsing Git URL '#url' should resolve to project '#expected' with hub '#expectedHub'"() {
        given: 'an asset manager instance'
        def manager = new AssetManager()

        expect: 'URL is correctly parsed'
        manager.resolveNameFromGitUrl(url) == expected
        and: 'hub is correctly identified'
        manager.hub == expectedHub

        where:
        url                                         | expected           | expectedHub
        'https://github.com/foo/bar.git'           | 'foo/bar'         | 'github'
        'https://gitlab.com/foo/bar.git'           | 'foo/bar'         | 'gitlab'
        'file:///path/to/repo.git'                 | 'local/repo'      | 'file:/path/to'
        'https://custom.gitea.org/foo/bar.git'     | 'foo/bar'         | 'gitea'
    }

    @Unroll
    def "parsing Git config should extract correct remote URL and domain"() {
        given: 'a temporary git directory with config'
        def dir = tempDir.root
        dir.resolve('.git').mkdir()
        dir.resolve('.git/config').text = configText

        and: 'an asset manager instance'
        def manager = new AssetManager().setLocalPath(dir.toFile())

        expect: 'remote URL and domain are correctly parsed'
        manager.getGitConfigRemoteUrl() == expectedUrl
        and: 'domain is correctly extracted'
        manager.getGitConfigRemoteDomain() == expectedDomain

        where:
        configText                                                           | expectedUrl                                    | expectedDomain
        GIT_CONFIG_TEXT                                                     | 'https://github.com/nextflow-io/nextflow.git' | 'github.com'
        GIT_CONFIG_LONG                                                     | 'git@github.com:nextflow-io/nextflow.git'     | 'github.com'
        '''[remote "origin"]
            url = https://gitlab.com/user/repo.git'''                       | 'https://gitlab.com/user/repo.git'            | 'gitlab.com'
        '''[remote "origin"]
            url = file:///local/path/repo.git'''                           | 'file:///local/path/repo.git'                 | '/local/path'
    }

    // ----------------------------------------------------------------
    // Script and manifest handling tests
    // ----------------------------------------------------------------

    @Unroll
    def "resolving script name with #description"() {
        given: 'a test repository directory'
        def dir = tempDir.getRoot()
        dir.resolve('test/repo').mkdirs()

        and: 'optional manifest configuration'
        if (configContent) {
            dir.resolve('test/repo/nextflow.config').text = configContent
        }

        and: 'git configuration'
        dir.resolve('test/repo/.git').mkdir()
        dir.resolve('test/repo/.git/config').text = GIT_CONFIG_TEXT

        and: 'an asset manager instance'
        def manager = new AssetManager()
        manager.setLocalPath(dir.resolve('test/repo').toFile())

        expect: 'script name is correctly resolved'
        manager.getMainScriptName() == expected

        where:
        description                 | scriptName | configContent                            | expected
        'default settings'         | null       | null                                    | 'main.nf'
        'manifest script setting'  | null       | 'manifest { mainScript = "custom.nf" }' | 'custom.nf'
        'manifest overrides input' | 'test.nf'  | 'manifest { mainScript = "custom.nf" }' | 'custom.nf'
    }

    // ----------------------------------------------------------------
    // Branch and reference management tests
    // ----------------------------------------------------------------

    @Unroll
    def "branch reference '#refName' should #description"() {
        given: 'a mock Git reference'
        def ref = Mock(Ref)
        ref.name >> refName

        expect: 'reference is correctly classified'
        AssetManager.isRemoteBranch(ref) == expected

        where:
        refName                       | expected | description
        'refs/remotes/origin/master' | true     | 'be identified as remote branch'
        'refs/remotes/origin/HEAD'   | false    | 'not be identified as remote branch'
        'refs/heads/master'          | false    | 'not be identified as remote branch'
        'refs/tags/v1.0'             | false    | 'not be identified as remote branch'
    }

    // ----------------------------------------------------------------
    // Integration tests with real repositories
    // ----------------------------------------------------------------

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should download and pull repository"() {
        given: 'a GitHub token and asset manager'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'downloading for the first time'
        manager.download()

        then: 'repository is cloned successfully'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'downloading again'
        manager.download()

        then: 'no errors occur'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should handle tag downloads correctly"() {
        given: 'a GitHub token and asset manager'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'downloading a specific tag'
        manager.download("v1.2")

        then: 'repository is cloned successfully'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'downloading the same tag again'
        manager.download("v1.2")

        then: 'no errors occur'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should handle commit hash downloads correctly"() {
        given: 'a GitHub token and asset manager'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'downloading a specific commit'
        manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")

        then: 'repository is cloned successfully'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'downloading the same commit again'
        def result = manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")

        then: 'repository is up-to-date'
        noExceptionThrown()
        result == "Already-up-to-date"
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should handle branch downloads correctly"() {
        given: 'a GitHub token and asset manager'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'downloading a specific branch'
        manager.download("mybranch")

        then: 'repository is cloned successfully'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'downloading the same branch again'
        manager.download("mybranch")

        then: 'no errors occur'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should handle tag-to-branch transitions"() {
        given: 'a GitHub token and asset manager'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'downloading a tag first'
        manager.download("v1.2")

        then: 'repository is cloned successfully'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'downloading without revision (defaults to branch)'
        manager.download()

        then: 'no errors occur'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should clone repository to specified directory"() {
        given: 'a GitHub token and asset manager'
        def dir = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers:[github: [auth: token]]])

        when: 'cloning the repository'
        manager.clone(dir.toFile())

        then: 'repository is cloned with expected files'
        dir.resolve('README.md').exists()
        dir.resolve('.git').isDirectory()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should identify default branch correctly"() {
        given: 'a GitHub token and asset manager'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/socks', [providers: [github: [auth: token]]])

        when: 'downloading without revision'
        manager.download()
        manager.checkout(null)

        then: 'correct default branch is identified'
        folder.resolve('nextflow-io/socks/.git').isDirectory()
        manager.getCurrentRevision() == 'main'

        when: 'downloading again'
        manager.download()

        then: 'no errors occur'
        noExceptionThrown()
    }

    // ----------------------------------------------------------------
    // Manifest handling tests
    // ----------------------------------------------------------------

    def "should read manifest file correctly"() {
        given: 'a test repository with manifest'
        def config = '''
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

        when: 'building the asset manager'
        def holder = new AssetManager()
        holder.build('foo/bar')

        then: 'manifest values are correctly read'
        holder.getMainScriptName() == 'hello.nf'
        holder.manifest.getDefaultBranch() == 'super-stuff'
        holder.manifest.getHomePage() == 'http://foo.com'
        holder.manifest.getDescription() == 'This pipeline do this and that'
        holder.manifest.getAuthor() == 'Hi Dude'
    }

    def "should handle default manifest values"() {
        given: 'a test repository with minimal config'
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = 'empty: 1'
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when: 'building the asset manager'
        def holder = new AssetManager()
        holder.build('foo/bar')

        then: 'default values are used'
        holder.getMainScriptName() == 'main.nf'
        holder.getHomePage() == 'https://github.com/foo/bar'
        holder.manifest.getDefaultBranch() == null
        holder.manifest.getDescription() == null
    }

    def "should create script file object correctly"() {
        given: 'a test repository with script files'
        def dir = tempDir.root
        dir.resolve('main.nf').text = "println 'Hello world'"
        dir.resolve('nextflow.config').text = 'manifest {  }'
        dir.resolve('foo.nf').text = 'this is foo content'

        and: 'a git repository'
        def init = Git.init()
        def repo = init.setDirectory(dir.toFile()).call()
        repo.add().addFilepattern('.').call()
        def commit = repo.commit().setSign(false).setAll(true).setMessage('First commit').call()
        repo.close()

        and: 'git configuration'
        dir.resolve('.git/config').text = GIT_CONFIG_TEXT

        when: 'creating script file for main script'
        def p = Mock(RepositoryProvider) { getRepositoryUrl() >> 'https://github.com/nextflow-io/nextflow' }
        def manager = new AssetManager(provider: p)
                .setLocalPath(dir.toFile())
                .setProject('nextflow-io/nextflow')
        def script = manager.getScriptFile()

        then: 'script file object is correctly configured'
        script.localPath == dir
        script.commitId == commit.name()
        script.revision == getLocalDefaultBranch()
        script.parent == dir
        script.text == "println 'Hello world'"
        script.repository == 'https://github.com/nextflow-io/nextflow'
        script.projectName == 'nextflow-io/nextflow'

        when: 'creating script file for secondary script'
        p = Mock(RepositoryProvider) { getRepositoryUrl() >> 'https://github.com/nextflow-io/nextflow' }
        manager = new AssetManager(provider: p)
                .setLocalPath(dir.toFile())
                .setProject('nextflow-io/nextflow')
        script = manager.getScriptFile('foo.nf')

        then: 'script file object is correctly configured'
        script.localPath == dir
        script.commitId == commit.name()
        script.revision == getLocalDefaultBranch()
        script.parent == dir
        script.text == "this is foo content"
        script.repository == 'https://github.com/nextflow-io/nextflow'
        script.projectName == 'nextflow-io/nextflow'
    }

    // ----------------------------------------------------------------
    // Project listing and finding tests
    // ----------------------------------------------------------------

    def "should list and find projects correctly"() {
        given: 'a set of test repositories'
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        when: 'listing all projects'
        def list = AssetManager.list()

        then: 'correct list is returned'
        list.sort() == ['cbcrg/pipe1','cbcrg/pipe2','ncbi/blast']

        expect: 'finding specific projects works'
        AssetManager.find('blast') == 'ncbi/blast'
        AssetManager.find('pipe1') == 'cbcrg/pipe1'
        AssetManager.find('pipe') as Set == ['cbcrg/pipe1', 'cbcrg/pipe2'] as Set
    }

}
