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

import static nextflow.scm.AssetManager.BARE_REPO
import static nextflow.scm.AssetManager.REVISION_MAP
import static nextflow.scm.AssetManager.REVISION_SUBDIR

import spock.lang.IgnoreIf

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

    def testListRevisions() {
        given:
        String revisionMap1 = '''branch1,12345\nbranch2,67890'''
        String revisionMap2 = '''branchA,abcde\nbranchB,fghij'''

        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1/.nextflow/').mkdirs()
        folder.resolve('cbcrg/pipe2/.nextflow/').mkdirs()
        folder.resolve('cbcrg/pipe1/' + REVISION_MAP).text = revisionMap1
        folder.resolve('cbcrg/pipe2/' + REVISION_MAP).text = revisionMap2

        when:
        def manager = new AssetManager('cbcrg/pipe1')
        def list = manager.listRevisions()
        then:
        list == ['branch1','branch2']

        when:
        manager = new AssetManager('cbcrg/pipe2')
        list = manager.listRevisions()
        then:
        list == ['branchA', 'branchB']
    }

    def testListRevisionsAndCommits() {
        given:
        String revisionMap1 = '''branch1,12345\nbranch2,67890'''
        String revisionMap2 = '''branchA,abcde\nbranchB,fghij'''

        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1/.nextflow/').mkdirs()
        folder.resolve('cbcrg/pipe2/.nextflow/').mkdirs()
        folder.resolve('cbcrg/pipe1/' + REVISION_MAP).text = revisionMap1
        folder.resolve('cbcrg/pipe2/' + REVISION_MAP).text = revisionMap2

        when:
        def manager = new AssetManager('cbcrg/pipe1')
        def dict = manager.listRevisionsAndCommits()
        then:
        dict == Map.of('branch1','12345','branch2','67890')

        when:
        manager = new AssetManager('cbcrg/pipe2')
        dict = manager.listRevisionsAndCommits()
        then:
        dict == Map.of('branchA','abcde','branchB','fghij')
    }

    def testListCommits() {
        given:
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1/' + REVISION_SUBDIR + '/12345').mkdirs()
        folder.resolve('cbcrg/pipe1/' + REVISION_SUBDIR + '/67890').mkdirs()
        folder.resolve('cbcrg/pipe2/' + REVISION_SUBDIR + '/abcde').mkdirs()
        folder.resolve('cbcrg/pipe2/' + REVISION_SUBDIR + '/fghij').mkdirs()

        when:
        def manager = new AssetManager('cbcrg/pipe1')
        def list = manager.listCommits()
        then:
        list.sort() == ['12345','67890']

        when:
        manager = new AssetManager('cbcrg/pipe2')
        list = manager.listCommits()
        then:
        list.sort() == ['abcde','fghij']
    }


    def testResolveName() {

        given:
        def folder = tempDir.getRoot()
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        def manager = new AssetManager()

        when:
        def result = manager.resolveName('x/y')
        then:
        result == 'x/y'

        when:
        result = manager.resolveName('blast')
        then:
        result == 'ncbi/blast'

        when:
        result = manager.resolveName('ncbi/blast/script.nf')
        then:
        result == 'ncbi/blast'

        when:
        result = manager.resolveName('blast/script.nf')
        then:
        result == 'ncbi/blast'

        when:
        manager.resolveName('pipe')
        then:
        thrown(AbortOperationException)

        when:
        manager.resolveName('pipe/alpha/beta')
        then:
        thrown(AbortOperationException)

        when:
        result = manager.resolveName('../blast/script.nf')
        then:
        thrown(AbortOperationException)

        when:
        result = manager.resolveName('./blast/script.nf')
        then:
        thrown(AbortOperationException)

    }


    def testUpdateRevisionMap() {

        given:
        def folder = tempDir.getRoot()
        String revision = null
        def manager = new AssetManager().build('nextflow-io/hello')
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)
        String revisionMap1 = '''v1.2,1b420d060d3fad67027154ac48e3bdea06f058da\n'''

        when:
        manager.updateRevisionMap('v1.2','1b420d060d3fad67027154ac48e3bdea06f058da')
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_MAP).exists()
        folder.resolve('nextflow-io/hello/' + REVISION_MAP).text == revisionMap1
    }


    def testRevisionToCommitWithMap() {

        given:
        def folder = tempDir.getRoot()
        String revision = null
        def manager = new AssetManager().build('nextflow-io/hello')
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)
        String revisionMap1 = '''v1.2,1b420d060d3fad67027154ac48e3bdea06f058da\n'''

        when:
        folder.resolve('nextflow-io/hello/.nextflow').mkdirs()
        folder.resolve('nextflow-io/hello/' + REVISION_MAP).text = revisionMap1
        then:
        manager.revisionToCommitWithMap('v1.2') == '1b420d060d3fad67027154ac48e3bdea06f058da'
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testCloneBareRepo() {

        given:
        def folder = tempDir.getRoot()
        String revision = null
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)

        when:
        manager.checkBareRepo()
        then:
        folder.resolve('nextflow-io/hello/' + BARE_REPO).isDirectory()
        folder.resolve('nextflow-io/hello/' + BARE_REPO + '/config').exists()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testRevisionToCommitWithBareRepo() {

        given:
        def folder = tempDir.getRoot()
        String revision = null
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)

        when:
        manager.checkBareRepo()
        then:
        manager.revisionToCommitWithBareRepo('v1.2') == '1b420d060d3fad67027154ac48e3bdea06f058da'
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPull() {

        given:
        def folder = tempDir.getRoot()
        String revision = '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' // easier to fix commit for a generic test
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' + '/.git').isDirectory()

        when:
        def result = manager.download()
        then:
        noExceptionThrown()
        result == "Already-up-to-date"
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullTagTwice() {

        given:
        def folder = tempDir.getRoot()
        String revision = 'v1.2'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)
        // tag v1.2 -> commit 1b420d060d3fad67027154ac48e3bdea06f058da

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '1b420d060d3fad67027154ac48e3bdea06f058da' + '/.git').isDirectory()

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

    // The hashes used here are NOT associated with tags.
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullHashTwice() {

        given:
        def folder = tempDir.getRoot()
        String revision = '6b9515aba6c7efc6a9b3f273ce116fc0c224bf68'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '6b9515aba6c7efc6a9b3f273ce116fc0c224bf68' + '/.git').isDirectory()

        when:
        def result = manager.download()
        then:
        noExceptionThrown()
        result == "Already-up-to-date"
    }


    // Downloading a branch first and then pulling the branch
    // should work fine, unlike with tags.
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPullBranchTwice() {

        given:
        def folder = tempDir.getRoot()
        String revision = 'mybranch'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)
        // as of Jun 2024, branch "mybranch" -> commit "1c3e9e7404127514d69369cd87f8036830f5cf64"

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '1c3e9e7404127514d69369cd87f8036830f5cf64' + '/.git').isDirectory()

        when:
        def result = manager.download()
        then:
        noExceptionThrown()
        result == "Already-up-to-date"
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testClone() {

        given:
        def dir = tempDir.getRoot()
        String revision = null
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers:[github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)

        when:
        manager.clone(dir.toFile())

        then:
        dir.resolve('README.md').exists()
        dir.resolve('.git').isDirectory()

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
        dir.resolve('foo/bar/' + REVISION_SUBDIR + '/' + 'mockup_dir').mkdirs()
        dir.resolve('foo/bar/' + REVISION_SUBDIR + '/' + 'mockup_dir' + '/nextflow.config').text = config
        dir.resolve('foo/bar/' + BARE_REPO).mkdirs()
        dir.resolve('foo/bar/' + BARE_REPO + '/config').text = GIT_CONFIG_TEXT

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')
        holder.setLocalPath(new File(dir.toString() + '/foo/bar/' + REVISION_SUBDIR + '/' + 'mockup_dir'))
        then:
        holder.getMainScriptName() == 'hello.nf'
        holder.manifest.getDefaultBranch() == 'super-stuff'
        holder.manifest.getHomePage() == 'http://foo.com'
        holder.manifest.getDescription() == 'This pipeline do this and that'
        holder.manifest.getAuthor() == 'Hi Dude'

    }

    def 'should return default main script file' () {

        given:
        def dir = tempDir.getRoot()
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = 'empty: 1'

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')
        holder.setLocalPath(new File(dir.toString() + '/foo/bar'))

        then:
        holder.getMainScriptName() == 'main.nf'
        holder.getHomePage() == 'https://github.com/foo/bar'
        holder.manifest.getDefaultBranch() == 'master'
        holder.manifest.getDescription() == null

    }

    def 'should parse git config and return the remote url' () {

        given:
        def dir = tempDir.root
        dir.resolve('foo/bar/' + BARE_REPO).mkdirs()
        dir.resolve('foo/bar/' + BARE_REPO + '/config').text = GIT_CONFIG_LONG

        when:
        def manager = new AssetManager().build('foo/bar')
        then:
        manager.getGitConfigRemoteUrl() == 'git@github.com:nextflow-io/nextflow.git'

    }

    def 'should parse git config and return the remote host' () {

        given:
        def dir = tempDir.root
        dir.resolve('foo/bar/' + BARE_REPO).mkdirs()
        dir.resolve('foo/bar/' + BARE_REPO + '/config').text = GIT_CONFIG_LONG

        when:
        def manager = new AssetManager().build('foo/bar')
        then:
        manager.getGitConfigRemoteDomain() == 'github.com'

    }

    def 'should create a script file object' () {

        given:
        def dir = tempDir.root
        // create the repo dir
        dir.resolve('main.nf').text = "println 'Hello world'"
        dir.resolve('nextflow.config').text = 'manifest {  }'
        dir.resolve('foo.nf').text = 'this is foo content'

        def init = Git.init()
        def repo = init.setDirectory( dir.toFile() ).call()
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

    def 'should return project name from git url' () {

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
    def 'should download branch specified'() {

        given:
        def folder = tempDir.getRoot()
        String revision = 'dev'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/nf-test-branch', revision)
        // as of June 2024, branch "dev" -> commit "6f882561d589365c3950d170df8445e3c0dc8028"

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/' + '6f882561d589365c3950d170df8445e3c0dc8028' + '/.git').isDirectory()
        and:
        folder.resolve('nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/' + '6f882561d589365c3950d170df8445e3c0dc8028' + '/workflow.nf').text == "println 'Hello'\n"
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should fetch main script from branch specified'() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        String revision = 'dev'
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)

        expect:
        manager.checkValidRemoteRepo()
        and:
        manager.getMainScriptName() == 'workflow.nf'
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should download tag specified'() {

        given:
        def folder = tempDir.getRoot()
        String revision = 'v0.1'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])
        manager.setRevisionAndLocalPath('nextflow-io/hello', revision)
        // tag "v0.1" -> commit "6f882561d589365c3950d170df8445e3c0dc8028"

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/' + '6f882561d589365c3950d170df8445e3c0dc8028' + '/.git').isDirectory()
        and:
        folder.resolve('nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/' + '6f882561d589365c3950d170df8445e3c0dc8028' + '/workflow.nf').text == "println 'Hello'\n"
    }

}
