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
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Config
import org.junit.Rule
import spock.lang.Ignore
import spock.lang.PendingFeature
import spock.lang.Requires
import spock.lang.Specification
import test.TemporaryPath
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Test suite for AssetManager class which handles Git repository operations
 * and pipeline asset management.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class AssetManagerTest extends Specification {

    static final String GIT_CONFIG_TEXT = '''
            [remote "origin"]
                url = https://github.com/nextflow-io/nextflow.git
                fetch = +refs/heads/*:refs/remotes/origin/*
            [branch "master"]
                remote = origin
                merge = refs/heads/master
            '''.stripIndent()

    static final String GIT_CONFIG_LONG = '''
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
        '''.stripIndent()

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def setup() {
        AssetManager.root = tempDir.root.toFile()
    }

    // Helper method to grab the default branch if set in ~/.gitconfig
    private String getLocalDefaultBranch() {
        def defaultBranch = 'master'
        def gitconfig = Paths.get(System.getProperty('user.home'),'.gitconfig')
        if(gitconfig.exists()) {
            def config = new Config()
            config.fromText(gitconfig.text)
            defaultBranch = config.getString('init', null, 'defaultBranch') ?: 'master'
        }
        defaultBranch
    }

    def "should list available pipeline assets"() {
        given: 'A set of pipeline directories'
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        when: 'Listing available assets'
        def list = AssetManager.list()

        then: 'Should return sorted list of pipelines'
        list.sort() == ['cbcrg/pipe1','cbcrg/pipe2','ncbi/blast']

        and: 'Should find specific pipelines'
        AssetManager.find('blast') == 'ncbi/blast'
        AssetManager.find('pipe1') == 'cbcrg/pipe1'
        AssetManager.find('pipe') as Set == ['cbcrg/pipe1', 'cbcrg/pipe2'] as Set
    }

    def "should resolve pipeline names correctly"() {
        given: 'A set of pipeline directories and manager instance'
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()
        def manager = new AssetManager()

        when: 'Resolving exact path'
        def result = manager.resolveName('x/y')
        then: 'Should return exact path'
        result == 'x/y'

        when: 'Resolving simple name'
        result = manager.resolveName('blast')
        then: 'Should return full path'
        result == 'ncbi/blast'

        when: 'Resolving path with script'
        result = manager.resolveName('ncbi/blast/script.nf')
        then: 'Should return base path'
        result == 'ncbi/blast'

        when: 'Resolving script with simple name'
        result = manager.resolveName('blast/script.nf')
        then: 'Should return full path'
        result == 'ncbi/blast'

        when: 'Resolving ambiguous name'
        manager.resolveName('pipe')
        then: 'Should throw exception'
        thrown(AbortOperationException)

        when: 'Resolving deep path'
        manager.resolveName('pipe/alpha/beta')
        then: 'Should throw exception'
        thrown(AbortOperationException)

        when: 'Resolving path with parent reference'
        result = manager.resolveName('../blast/script.nf')
        then: 'Should throw exception'
        thrown(AbortOperationException)

        when: 'Resolving path with current directory reference'
        result = manager.resolveName('./blast/script.nf')
        then: 'Should throw exception'
        thrown(AbortOperationException)
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should pull repository successfully"() {
        given: 'Asset manager with GitHub token'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'Downloading repository'
        manager.download()
        then: 'Git repository should exist'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'Downloading repository again'
        manager.download()
        then: 'Should not throw exception'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should pull tag multiple times successfully"() {
        given: 'Asset manager with GitHub token'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'Downloading specific tag'
        manager.download("v1.2")
        then: 'Git repository should exist'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'Downloading same tag again'
        manager.download("v1.2")
        then: 'Should not throw exception'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should pull commit hash multiple times successfully"() {
        given: 'Asset manager with GitHub token'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'Downloading specific commit'
        manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")
        then: 'Git repository should exist'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'Downloading same commit again'
        def result = manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")
        then: 'Should not throw exception and indicate no changes'
        noExceptionThrown()
        result == "Already-up-to-date"
    }

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

    def testCreateProviderFor(){

        when:
        def manager = new AssetManager()
        def repo = manager.createHubProvider('github')
        then:
        repo instanceof GithubRepositoryProvider

        when:
        manager = new AssetManager()
        repo = manager.createHubProvider('bitbucket')
        then:
        repo instanceof BitbucketRepositoryProvider

        when:
        manager = new AssetManager()
        repo = manager.createHubProvider('gitlab')
        then:
        repo instanceof GitlabRepositoryProvider

        when:
        manager = [:] as AssetManager
        manager.createHubProvider('xxx')
        then:
        thrown(AbortOperationException)

    }


    def "should read manifest file correctly"() {
        given: 'A manifest configuration file'
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

        when: 'Building asset manager with manifest'
        def holder = new AssetManager()
        holder.build('foo/bar')

        then: 'Manifest values should be correctly read'
        holder.getMainScriptName() == 'hello.nf'
        holder.manifest.getDefaultBranch() == 'super-stuff'
        holder.manifest.getHomePage() == 'http://foo.com'
        holder.manifest.getDescription() == 'This pipeline do this and that'
        holder.manifest.getAuthor() == 'Hi Dude'
    }

    def "should read default tag from manifest file"() {
        given: 'A manifest configuration with default tag'
        def config = '''
            manifest {
                homePage = 'http://foo.com'
                mainScript = 'hello.nf'
                defaultBranch = '1.0.0'
                description = 'This pipeline do this and that'
                author = 'Hi Dude'
            }
            '''
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = config
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when: 'Building asset manager'
        def holder = new AssetManager()
        holder.build('foo/bar')

        then: 'Default branch should be the tag'
        holder.manifest.getDefaultBranch() == '1.0.0'
    }

    def "should return default main script file"() {
        given: 'A basic configuration file'
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = 'empty: 1'
        dir.resolve('foo/bar/.git').mkdir()
        dir.resolve('foo/bar/.git/config').text = GIT_CONFIG_TEXT

        when: 'Building asset manager'
        def holder = new AssetManager()
        holder.build('foo/bar')

        then: 'Should return default values'
        holder.getMainScriptName() == 'main.nf'
        holder.getHomePage() == 'https://github.com/foo/bar'
        holder.manifest.getDefaultBranch() == null
        holder.manifest.getDescription() == null
    }

    def "should parse git config and return remote URL"() {
        given: 'A Git configuration file'
        def dir = tempDir.root
        dir.resolve('.git').mkdir()
        dir.resolve('.git/config').text = GIT_CONFIG_LONG

        when: 'Getting remote URL'
        def manager = new AssetManager().setLocalPath(dir.toFile())

        then: 'Should return correct remote URL'
        manager.getGitConfigRemoteUrl() == 'git@github.com:nextflow-io/nextflow.git'
    }

    def "should parse git config and return remote domain"() {
        given: 'A Git configuration file'
        def dir = tempDir.root
        dir.resolve('.git').mkdir()
        dir.resolve('.git/config').text = GIT_CONFIG_LONG

        when: 'Getting remote domain'
        def manager = new AssetManager().setLocalPath(dir.toFile())

        then: 'Should return correct remote domain'
        manager.getGitConfigRemoteDomain() == 'github.com'
    }

    def "should create script file object correctly"() {
        given: 'A Git repository with files'
        def dir = tempDir.root
        dir.resolve('main.nf').text = "println 'Hello world'"
        dir.resolve('nextflow.config').text = 'manifest {  }'
        dir.resolve('foo.nf').text = 'this is foo content'

        and: 'Initialize Git repository'
        def init = Git.init()
        def repo = init.setDirectory(dir.toFile()).call()
        repo.add().addFilepattern('.').call()
        def commit = repo.commit().setSign(false).setAll(true).setMessage('First commit').call()
        repo.close()

        and: 'Add remote configuration'
        dir.resolve('.git/config').text = GIT_CONFIG_TEXT

        when: 'Creating script file object for main script'
        def p = Mock(RepositoryProvider) { getRepositoryUrl() >> 'https://github.com/nextflow-io/nextflow' }
        def manager = new AssetManager(provider: p)
                .setLocalPath(dir.toFile())
                .setProject('nextflow-io/nextflow')
        def script = manager.getScriptFile()

        then: 'Script file object should have correct properties'
        script.localPath == dir
        script.commitId == commit.name()
        script.revision == getLocalDefaultBranch()
        script.parent == dir
        script.text == "println 'Hello world'"
        script.repository == 'https://github.com/nextflow-io/nextflow'
        script.projectName == 'nextflow-io/nextflow'

        when: 'Creating script file object for specific script'
        p = Mock(RepositoryProvider) { getRepositoryUrl() >> 'https://github.com/nextflow-io/nextflow' }
        manager = new AssetManager(provider: p)
                .setLocalPath(dir.toFile())
                .setProject('nextflow-io/nextflow')
        script = manager.getScriptFile('foo.nf')

        then: 'Script file object should have correct properties'
        script.localPath == dir
        script.commitId == commit.name()
        script.revision == getLocalDefaultBranch()
        script.parent == dir
        script.text == "this is foo content"
        script.repository == 'https://github.com/nextflow-io/nextflow'
        script.projectName == 'nextflow-io/nextflow'
    }

    def "should return project name from git url"() {

        AssetManager manager
        String result

        when:
        manager = new AssetManager()
        result = manager.resolveNameFromGitUrl('nextflow/pipe')
        then:
        result == null
        manager.hub == null

        when:
        manager = new AssetManager()
        result = manager.resolveNameFromGitUrl('https://gitlab.com/pditommaso/hello.git')
        then:
        result == 'pditommaso/hello'
        manager.hub == 'gitlab'

        when:
        manager = new AssetManager()
        result = manager.resolveNameFromGitUrl('file:/user/repo/projects/hello.git')
        then:
        result == 'local/hello'
        manager.hub == 'file:/user/repo/projects'

        when:
        manager = new AssetManager()
        manager.providerConfigs.add( new ProviderConfig('local-scm', [platform: 'github', server: 'http://foo.bar.com']) )
        result = manager.resolveNameFromGitUrl('https://foo.bar.com/project/xyz.git')
        then:
        result == 'project/xyz'
        manager.hub == 'local-scm'

        when:
        manager = new AssetManager()
        manager.providerConfigs.add( new ProviderConfig('gitea', [platform: 'gitea', server: 'http://my-server.org/sub1']) )
        result = manager.resolveNameFromGitUrl('http://my-server.org/sub1/foo/bar.git')
        then:
        result == 'foo/bar'
        manager.hub == 'gitea'

    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should download branch specified"() {

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

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should fetch main script from branch specified"() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])

        expect:
        manager.checkValidRemoteRepo('dev')
        and:
        manager.getMainScriptName() == 'workflow.nf'

    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should download tag specified"() {

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

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should use a default tag"() {

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

    def "should work with defaultBranch = master"() {
        given:
        def config = '''
            manifest {
                defaultBranch = 'master'
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
        holder.manifest.getDefaultBranch() == 'master'
        holder.manifest.getDefaultRevision() == 'master'
    }

    def "should work with defaultRevision"() {
        given:
        def config = '''
            manifest {
                defaultRevision = '1.0.0'
                defaultBranch = 'master'
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
        holder.manifest.getDefaultRevision() == '1.0.0'
        holder.manifest.getDefaultBranch() == 'master'
    }

    def "should use version as defaultRevision when available"() {
        given:
        def config = '''
            manifest {
                version = '2.0.0'
                defaultBranch = 'master'
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
        holder.manifest.getVersion() == '2.0.0'
        holder.manifest.getDefaultRevision() == '2.0.0'
        holder.manifest.getDefaultBranch() == 'master'
    }

    def "should prioritize defaultRevision over version"() {
        given:
        def config = '''
            manifest {
                version = '2.0.0' // Development version for main branch
                defaultRevision = '1.0.0' // Latest stable version
                defaultBranch = 'master'
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
        holder.manifest.getVersion() == '2.0.0'
        holder.manifest.getDefaultRevision() == '1.0.0'
        holder.manifest.getDefaultBranch() == 'master'
    }

    def "should work with no defaultBranch"() {
        given:
        def config = '''
            manifest {
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
        holder.manifest.getDefaultBranch() == 'master'
        holder.manifest.getDefaultRevision() == 'master'
    }

    def "should default to version tag if manifest version and no defaultBranch"() {
        given:
        def config = '''
            manifest {
                version = '3.0.0'
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
        holder.manifest.getVersion() == '3.0.0'
        holder.manifest.getDefaultRevision() == '3.0.0'
        holder.manifest.getDefaultBranch() == 'master'
    }

    def "should handle commit hash in defaultRevision"() {
        given:
        def config = '''
            manifest {
                defaultRevision = '6b9515aba6c7efc6a9b3f273ce116fc0c224bf68'
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
        holder.manifest.getDefaultRevision() == '6b9515aba6c7efc6a9b3f273ce116fc0c224bf68'
        holder.manifest.getDefaultBranch() == 'master'
    }

    @PendingFeature
    def "should not warn if project uses a tag as a defaultBranch"() {
        given:
        def ENV = [FOO: '/something', NXF_DEBUG: 'true']

        when:
        new CmdRun(revision: 'xyz')

        then:
        def warning = capture
                .toString()
                .readLines()
                .findResults { line -> line.contains('WARN') ? line : null }
                .join('\n')
        and:
        !warning
        noExceptionThrown()
    }

    @PendingFeature
    def "should handle development version with stable defaultRevision"() {
        given:
        def config = '''
        manifest {
            version = '2.3.0dev'
            defaultRevision = '2.2.0'
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
        holder.manifest.getVersion() == '2.3.0dev'
        holder.manifest.getDefaultRevision() == '2.2.0'
        holder.manifest.isDevelopmentVersion() == true
    }

    @PendingFeature
    def "should correctly compare development and release versions"() {
        given:
        def config = '''
        manifest {
            version = '2.3.0dev'
            defaultRevision = '2.2.0'
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
        holder.manifest.isVersionGreaterThan('2.2.0') == true
        holder.manifest.isVersionCompatibleWith('2.2.0') == true
    }

    @PendingFeature
    def "should not warn if project uses a tag as a defaultBranch"() {
        given:
        def ENV = [FOO: '/something', NXF_DEBUG: 'true']

        when:
        new CmdRun(revision: 'xyz')

        then:
        def warning = capture
            .toString()
            .readLines()
            .findResults { line -> line.contains('WARN') ? line : null }
            .join('\n')
        and:
        !warning
        noExceptionThrown()
    }

    // Test version = '2.3.0-RC1' with defaultRevision
    def "should handle release candidate versions"() {
        given:
        def config = '''
    manifest {
        version = '2.3.0-RC1'
        defaultRevision = '2.2.0'
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
        holder.manifest.getVersion() == '2.3.0-RC1'
        holder.manifest.getDefaultRevision() == '2.2.0'
        holder.manifest.isVersionGreaterThan('2.2.0') == true
    }

    // Test version = '2.2.1-hotfix' with defaultRevision = '2.2.0'
    def "should handle hotfix versions"() {
        given:
        def config = '''
        manifest {
        version = '2.2.1-hotfix'
        defaultRevision = '2.2.0'
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
        holder.manifest.getVersion() == '2.2.1-hotfix'
        holder.manifest.getDefaultRevision() == '2.2.0'
        holder.manifest.isVersionGreaterThan('2.2.0') == true
    }

    // Test handling of feature branches while maintaining stable defaultRevision
    def "should support multiple development branches"() {
        given:
        def config = '''
        manifest {
            version = '2.3.0-dev'
            defaultRevision = '2.2.0'
            defaultBranch = 'feature/new-feature'
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
        holder.manifest.getVersion() == '2.3.0-dev'
        holder.manifest.getDefaultRevision() == '2.2.0'
        holder.manifest.getDefaultBranch() == 'feature/new-feature'
    }

    // Test downgrading defaultRevision for emergency rollbacks
    def "should handle version rollback scenarios"() {
        given:
        def config = '''
        manifest {
        version = '2.2.0'
        defaultRevision = '2.3.0'  // Attempting to rollback to older version
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
        holder.manifest.getVersion() == '2.2.0'
        holder.manifest.getDefaultRevision() == '2.3.0'
        holder.manifest.isVersionGreaterThan('2.3.0') == false
    }

    // Test that development version is always ahead of defaultRevision
    def "should validate version and defaultRevision compatibility"() {
        given:
        def config = '''
        manifest {
        version = '2.1.0'  // Version older than defaultRevision
        defaultRevision = '2.2.0'
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
        holder.manifest.getVersion() == '2.1.0'
        holder.manifest.getDefaultRevision() == '2.2.0'
        holder.manifest.isVersionGreaterThan('2.2.0') == false
        holder.manifest.isVersionCompatibleWith('2.2.0') == false
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should pull branch multiple times successfully"() {
        given: 'Asset manager with GitHub token'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'Downloading branch first time'
        manager.download("mybranch")
        then: 'Git repository should exist'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'Downloading branch second time'
        manager.download("mybranch")
        then: 'Should not throw exception'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should handle tag to branch transition"() {
        given: 'Asset manager with GitHub token'
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when: 'Downloading specific tag'
        manager.download("v1.2")
        then: 'Git repository should exist'
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when: 'Downloading without revision (default branch)'
        manager.download()
        then: 'Should not throw exception'
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def "should clone repository successfully"() {
        given: 'Asset manager with GitHub token'
        def dir = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers:[github: [auth: token]]])

        when: 'Cloning repository'
        manager.clone(dir.toFile())

        then: 'Repository files should exist'
        dir.resolve('README.md').exists()
        dir.resolve('.git').isDirectory()
    }

    def "should get script name correctly"() {
        given: 'Test directories and configurations'
        def dir = tempDir.getRoot()
        dir.resolve('sub1').mkdir()
        dir.resolve('sub1/nextflow.config').text = "manifest.mainScript = 'pippo.nf'"
        dir.resolve('sub2').mkdir()

        when: 'Getting script name from directory with config'
        def holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then: 'Should return configured script name'
        holder.getMainScriptName() == 'pippo.nf'

        when: 'Getting script name from directory without config'
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub2').toFile()
        then: 'Should return default script name'
        holder.getMainScriptName() == 'main.nf'

        when: 'Getting script name with resolved project'
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        holder.resolveName('nextflow/hello')
        then: 'Should return configured script name'
        holder.getMainScriptName() == 'pippo.nf'

        when: 'Getting script name with specific script path'
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        holder.resolveName('nextflow/hello/my-script.nf')
        then: 'Should return specified script name'
        holder.getMainScriptName() == 'my-script.nf'

        when: 'Getting script name with deep path'
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then: 'Should resolve name and return script path'
        holder.resolveName('nextflow-io/hello/x/y/z/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'x/y/z/my-script.nf'

        when: 'Getting script name with simple path'
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then: 'Should resolve name and return script name'
        holder.resolveName('nextflow-io/hello/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'my-script.nf'

        when: 'Getting script name with short path'
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then: 'Should resolve name and return script name'
        holder.resolveName('hello/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'my-script.nf'
    }

    def "should create provider for different platforms"() {
        when: 'Creating GitHub provider'
        def manager = new AssetManager()
        def repo = manager.createHubProvider('github')
        then: 'Should return GitHub provider'
        repo instanceof GithubRepositoryProvider

        when: 'Creating Bitbucket provider'
        manager = new AssetManager()
        repo = manager.createHubProvider('bitbucket')
        then: 'Should return Bitbucket provider'
        repo instanceof BitbucketRepositoryProvider

        when: 'Creating GitLab provider'
        manager = new AssetManager()
        repo = manager.createHubProvider('gitlab')
        then: 'Should return GitLab provider'
        repo instanceof GitlabRepositoryProvider

        when: 'Creating unknown provider'
        manager = [:] as AssetManager
        manager.createHubProvider('xxx')
        then: 'Should throw exception'
        thrown(AbortOperationException)
    }
}
