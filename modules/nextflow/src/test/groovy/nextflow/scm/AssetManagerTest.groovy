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

package nextflow.scm

import java.util.concurrent.CountDownLatch
import java.util.concurrent.Executors

import static nextflow.scm.MultiRevisionRepository.BARE_REPO
import static nextflow.scm.MultiRevisionRepository.REPOS_SUBDIR
import static nextflow.scm.MultiRevisionRepository.REVISION_SUBDIR

import spock.lang.IgnoreIf

import nextflow.scm.HubOptions
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Config
import org.junit.Rule
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Timeout
import test.OutputCapture
import test.TemporaryPath
import java.nio.file.Paths
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(60)
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

    @Rule
    OutputCapture capture = new OutputCapture()

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

    /**
     * Create a multi-revision repository for the given project with a single commit,
     * checked out in its clone directory, without accessing the network
     */
    private String createLocalRepo(String project, Map<String,String> files, String remoteUrl = "https://github.com/${project}.git") {
        def source = tempDir.root.resolve('_source/' + project)
        source.mkdirs()
        files.each { name, text -> source.resolve(name).text = text }
        def git = Git.init().setDirectory(source.toFile()).call()
        git.add().addFilepattern('.').call()
        def commitId = git.commit().setSign(false).setMessage('First commit').call().name()
        git.close()

        def projectDir = tempDir.root.resolve(REPOS_SUBDIR + '/' + project)
        def bare = Git.cloneRepository()
            .setBare(true)
            .setURI(source.toFile().absolutePath)
            .setGitDir(projectDir.resolve(BARE_REPO).toFile())
            .call()
        def config = bare.getRepository().getConfig()
        config.setString('remote', 'origin', 'url', remoteUrl)
        config.save()
        bare.close()

        Git.cloneRepository()
            .setURI(source.toFile().absolutePath)
            .setDirectory(projectDir.resolve(REVISION_SUBDIR + '/' + commitId).toFile())
            .call()
            .close()
        return commitId
    }

    def testList() {

        given:
        def folder = tempDir.getRoot()
        folder.resolve(REPOS_SUBDIR +'/cbcrg/pipe1').mkdirs()
        folder.resolve(REPOS_SUBDIR +'/cbcrg/pipe2').mkdirs()
        folder.resolve(REPOS_SUBDIR +'/ncbi/blast').mkdirs()
        folder.resolve(REPOS_SUBDIR +'/new/repo').mkdirs()
        and: 'leftover clones of the legacy layout are ignored'
        folder.resolve('ncbi/blast/.git').mkdirs()
        folder.resolve('legacy/pipe3/.git').mkdirs()

        when:
        def list = AssetManager.list()
        then:
        list.sort() == ['cbcrg/pipe1','cbcrg/pipe2','ncbi/blast', 'new/repo']

        expect:
        AssetManager.find('blast') == 'ncbi/blast'
        AssetManager.find('pipe1') == 'cbcrg/pipe1'
        AssetManager.find('pipe') as Set == ['cbcrg/pipe1', 'cbcrg/pipe2'] as Set
        AssetManager.find('repo') == 'new/repo'

    }


    def testResolveName() {

        given:
        def folder = tempDir.getRoot()
        folder.resolve(REPOS_SUBDIR + '/cbcrg/pipe1').mkdirs()
        folder.resolve(REPOS_SUBDIR + '/cbcrg/pipe2').mkdirs()
        folder.resolve(REPOS_SUBDIR + '/ncbi/blast').mkdirs()

        when:
        def result = new AssetManager().resolveName('x/y')
        then:
        result == 'x/y'

        when:
        result = new AssetManager().resolveName('blast')
        then:
        result == 'ncbi/blast'

        when:
        result = new AssetManager().resolveName('ncbi/blast/script.nf')
        then:
        result == 'ncbi/blast'

        when:
        result = new AssetManager().resolveName('blast/script.nf')
        then:
        result == 'ncbi/blast'

        when:
        new AssetManager().resolveName('pipe')
        then:
        thrown(AbortOperationException)

        when:
        new AssetManager().resolveName('pipe/alpha/beta')
        then:
        thrown(AbortOperationException)

        when:
        new AssetManager().resolveName('../blast/script.nf')
        then:
        thrown(AbortOperationException)

        when:
        new AssetManager().resolveName('./blast/script.nf')
        then:
        thrown(AbortOperationException)

        when:
        // mainScript already set via `-main-script` AND name ends with a script file
        def manager = new AssetManager()
        manager.mainScript = 'my-script.nf'
        manager.resolveName('ncbi/blast/script.nf')
        then:
        def err = thrown(AbortOperationException)
        err.message.contains('Project name must be a directory when main script is provided')
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'test download from tag twice with multi-revision'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("v1.2")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git').isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()
        manager.getLocalPath().toString() == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da').toString()
        when:
        manager.download("v1.2")
        then:
        noExceptionThrown()

    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'test download from hash twice with multi-revision'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/6b9515aba6c7efc6a9b3f273ce116fc0c224bf68/.git').isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()
        manager.getLocalPath().toString() == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/6b9515aba6c7efc6a9b3f273ce116fc0c224bf68').toString()

        when:
        def result = manager.download("6b9515aba6c7efc6a9b3f273ce116fc0c224bf68")
        then:
        noExceptionThrown()
        result == "Already-up-to-date"
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'test download from branch twice with multi-revision'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("mybranch")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git').isDirectory()
        manager.getLocalPath().toString() == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64').toString()
        when:
        manager.download("mybranch")
        then:
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'test download from branch twice with multi-revision private repo'() {

        given:
        def folder = tempDir.getRoot()
        def token  = System.getenv('NXF_AZURE_REPOS_TOKEN')

        def manager = new AssetManager().build('https://dev.azure.com/pditommaso/nf-azure-repo', [providers: [azurerepos: [ auth: token]]])

        when:
        manager.download("master")
        then:
        folder.resolve(REPOS_SUBDIR + '/pditommaso/nf-azure-repo/nf-azure-repo/' + BARE_REPO).isDirectory()
        when:
        manager.download("master")
        then:
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'test download tag then branch with multi-revision'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])

        when:
        manager.download("v1.2")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git').isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()
        manager.getLocalPath().toString() == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da').toString()

        when:
        manager.download("mybranch")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git').isDirectory()
        manager.getLocalPath().toString() == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64').toString()
        noExceptionThrown()
    }


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

    def 'should pass the scm config http client opts to the provider' () {
        given:
        def config = [
            httpClient: [connectTimeout: '30s', requestTimeout: '90s'],
            providers: [github: [user: 'foo', password: 'bar']] ]

        when:
        def manager = new AssetManager().build('nextflow-io/hello', config)
        def provider = manager.createHubProvider('github')
        then:
        provider.getHttpClientOpts().connectTimeout() == java.time.Duration.ofSeconds(30)
        provider.getHttpClientOpts().requestTimeout() == java.time.Duration.ofSeconds(90)
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
                    description = 'This pipeline do this and that'
                    author = 'Hi Dude'
                }
                '''
        createLocalRepo('foo/bar', ['nextflow.config': config])

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')
        then:
        holder.getMainScriptName() == 'hello.nf'
        holder.manifest.getHomePage() == 'http://foo.com'
        holder.manifest.getDescription() == 'This pipeline do this and that'
        holder.manifest.getAuthor() == 'Hi Dude'

    }

    def 'should propagate explicit main script' () {

        given:
        createLocalRepo('foo/bar', ['nextflow.config': 'empty: 1'])

        when:
        // pass an explicit script file via the new build(...) signature
        def holder = new AssetManager()
        holder.build('foo/bar', null, null, null, 'sub/entry.nf')
        then:
        holder.getMainScriptName() == 'sub/entry.nf'

        when:
        // exercise the new (pipelineName, revision, scriptFile, cliOpts) constructor
        holder = new AssetManager('foo/bar', null, 'sub/entry.nf')
        then:
        holder.getMainScriptName() == 'sub/entry.nf'
    }

    def 'should return default main script file' () {

        given:
        createLocalRepo('foo/bar', ['nextflow.config': 'empty: 1'])

        when:
        def holder = new AssetManager()
        holder.build('foo/bar')

        then:
        holder.getMainScriptName() == 'main.nf'
        holder.getHomePage() == 'https://github.com/foo/bar'
        holder.manifest.getDescription() == null

    }

    def 'should parse git config and return the remote url' () {

        given:
        def dir = tempDir.root
        dir.resolve('.git').mkdir()
        dir.resolve('.git/config').text = GIT_CONFIG_LONG

        when:
        def manager = new AssetManager().setLocalPath(dir.toFile())
        then:
        manager.getGitConfigRemoteUrl() == 'git@github.com:nextflow-io/nextflow.git'

    }

    def 'should parse git config and return the remote host' () {

        given:
        def dir = tempDir.root
        dir.resolve('.git').mkdir()
        dir.resolve('.git/config').text = GIT_CONFIG_LONG

        when:
        def manager = new AssetManager().setLocalPath(dir.toFile())
        then:
        manager.getGitConfigRemoteDomain() == 'github.com'

    }

    def 'should warn when manifest default branch differs from the remote default branch' () {

        given:
        createLocalRepo('foo/bar', ['main.nf': "println 'Hello'", 'nextflow.config': "manifest.defaultBranch = 'dev'"])

        when:
        new AssetManager().build('foo/bar').tryCheckout()
        then:
        capture.toString().contains("The `manifest.defaultBranch` setting is ignored -- running the default branch of the remote repository (`${getLocalDefaultBranch()}`). Use `-r dev` to run the `dev` branch")
    }

    def 'should not warn when manifest default branch matches or a revision is specified' () {

        given:
        def commitId = createLocalRepo('foo/bar', ['main.nf': "println 'Hello'", 'nextflow.config': "manifest.defaultBranch = '${defaultBranch}'"])

        when:
        def manager = new AssetManager().build('foo/bar')
        manager.tryCheckout(useRevision ? commitId : null)
        then:
        manager.getManifest().getDefaultBranch() == defaultBranch
        !capture.toString().contains('manifest.defaultBranch')

        where:
        defaultBranch               | useRevision
        getLocalDefaultBranch()     | false
        'other'                     | true
    }

    def 'should detect hub provider from the bare repository config' () {

        given:
        createLocalRepo('foo/bar', ['main.nf': "println 'Hello'"], 'https://gitlab.com/foo/bar.git')

        when:
        def manager = new AssetManager().build('foo/bar')
        then:
        manager.getHub() == 'gitlab'
        manager.getLocalGitConfig() == tempDir.root.resolve(REPOS_SUBDIR + '/foo/bar/' + BARE_REPO + '/config').toFile()
    }

    def 'should reject a hub provider that does not match the bare repository config' () {

        given:
        createLocalRepo('foo/bar', ['main.nf': "println 'Hello'"], 'https://gitlab.com/foo/bar.git')

        when:
        new AssetManager().build('foo/bar', null, new HubOptions('github', null))
        then:
        def e = thrown(AbortOperationException)
        e.message.contains('has already been downloaded from a different provider: `gitlab`')
    }

    def 'should create a script file object from a local directory' () {

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

        // append fake remote data
        dir.resolve('.git/config').text = GIT_CONFIG_TEXT

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

        when:
        // a repository URL ending with a `.nf` script file is not a valid repo URL
        manager = new AssetManager()
        manager.resolveNameFromGitUrl('https://github.com/foo/bar/main.nf')
        then:
        def err = thrown(AbortOperationException)
        err.message.contains('-main-script')
    }

    def 'should detect local scm source' () {
        when: 'provider is a LocalRepositoryProvider'
        def config = new ProviderConfig('local', [path: '/tmp/workflows'])
        def manager = new AssetManager(provider: new LocalRepositoryProvider('socks', config))
        then:
        manager.isLocalScmSource()

        when: 'provider is a remote RepositoryProvider'
        manager = new AssetManager(provider: Mock(RepositoryProvider))
        then:
        !manager.isLocalScmSource()

        when: 'provider is null'
        manager = new AssetManager()
        then:
        !manager.isLocalScmSource()
    }

    def 'should detect local scm for pipelineNames' () {
        when:
        def manager = new AssetManager('my-repo', new HubOptions(null, null))
        then:
        !manager.isLocalScmSource()

        when:
        manager = new AssetManager('my-project/repo',  new HubOptions(null, null))
        then:
        !manager.isLocalScmSource()

        when:
        manager = new AssetManager('https://github.com/nextflow-io/socks.git',  new HubOptions(null, null))
        then:
        !manager.isLocalScmSource()

        when:
        manager = new AssetManager( 'file:/path/my-project/repo.git', new HubOptions(null, null))
        then:
        manager.isLocalScmSource()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should download branch specified multi-revision'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])

        when:
        manager.download("dev")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/6f882561d589365c3950d170df8445e3c0dc8028/.git').isDirectory()
        and:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/6f882561d589365c3950d170df8445e3c0dc8028/workflow.nf').text == "println 'Hello'\n"

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

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

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should download tag specified with revision'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/nf-test-branch', [providers: [github: [auth: token]]])

        when:
        manager.download("v0.1")
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/6f882561d589365c3950d170df8445e3c0dc8028/.git').isDirectory()
        and:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/nf-test-branch/' + REVISION_SUBDIR + '/6f882561d589365c3950d170df8445e3c0dc8028/workflow.nf').text == "println 'Hello'\n"

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should identify default branch when downloading repo multi-revision'() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/socks', [providers: [github: [auth: token]]])

        when:
        // simulate calling `nextflow run nextflow-io/socks` without specifying a revision
        manager.download()
        manager.checkout(null)
        then:
        manager.getCurrentRevision() == 'main'

        when:
        manager.download()
        then:
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list revisions and commits'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def pipelineName ='nextflow-io/hello'
        def revision1 = 'v1.2'
        def revision2 = 'mybranch'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download(revision1)
        manager.download(revision2)

        when:
        def branchesAndTags = manager.getBranchesAndTags(false)
        then:
        def pulled = branchesAndTags.pulled as List
        pulled.size() == 2
        revision1 in pulled && revision2 in pulled
        def branches = branchesAndTags.branches as List<Map>
        def tags = branchesAndTags.tags as List<Map>
        tags.find {it.name == revision1 }.commitId == '0ec2ecd0ac13bc7e32594c0258ebce55e383d241'
        branches.find { it.name == revision2 }.commitId == '1c3e9e7404127514d69369cd87f8036830f5cf64'
    }

    def 'should not be initialized when no local repository exists'() {
        when:
        def manager = new AssetManager().build('test-org/test-repo')

        then:
        manager.isNotInitialized()
        manager.getLocalPath() == tempDir.root.resolve(REPOS_SUBDIR + '/test-org/test-repo/' + REVISION_SUBDIR + '/HEAD').toFile()
        !manager.getLocalPath().exists()
        !manager.isLocal()
        manager.getProjectPath() == tempDir.root.resolve(REPOS_SUBDIR + '/test-org/test-repo').toFile()
    }

    def 'should ignore a leftover legacy clone'() {
        given:
        def folder = tempDir.getRoot()
        def legacyPath = folder.resolve('test-org/test-repo')
        legacyPath.mkdirs()
        legacyPath.resolve('main.nf').text = "println 'Hello'"
        def repo = Git.init().setDirectory(legacyPath.toFile()).call()
        repo.close()
        legacyPath.resolve('.git/config').text = GIT_CONFIG_TEXT

        when:
        def manager = new AssetManager().build('test-org/test-repo')

        then:
        capture.toString().contains("The legacy pipeline cache at ${legacyPath} is no longer used")
        manager.isNotInitialized()
        !manager.getLocalPath().exists()
        !manager.getLocalPath().toPath().startsWith(legacyPath)
        !manager.isLocal()
        !manager.isRunnable()
        manager.getLocalGitConfig() == null
        manager.getBranchList() == []
        AssetManager.list() == []

        when:
        manager.drop()

        then:
        thrown(AbortOperationException)
        legacyPath.resolve('.git').exists()
        legacyPath.resolve('main.nf').exists()
    }

    def 'should use bare repository when it exists'() {
        given:
        def folder = tempDir.getRoot()
        def barePath = folder.resolve(REPOS_SUBDIR + '/test-org/test-repo/' + BARE_REPO)
        barePath.mkdirs()

        // Create a proper bare git repository
        def init = Git.init()
        init.setDirectory(barePath.toFile())
        init.setBare(true)
        def repo = init.call()
        repo.close()

        when:
        def manager = new AssetManager().build('test-org/test-repo')

        then:
        !manager.isNotInitialized()
        manager.getLocalGitConfig() == barePath.resolve('config').toFile()
    }

    def 'should ignore legacy clone when multi-revision repository exists'() {
        given:
        def folder = tempDir.getRoot()
        // Create legacy repository
        def legacyPath = folder.resolve('test-org/test-repo')
        legacyPath.mkdirs()
        def initLegacy = Git.init()
        def repoLegacy = initLegacy.setDirectory(legacyPath.toFile()).call()
        repoLegacy.close()
        legacyPath.resolve('.git/config').text = GIT_CONFIG_TEXT

        // Create bare repository
        def barePath = folder.resolve(REPOS_SUBDIR + '/test-org/test-repo/' + BARE_REPO)
        barePath.mkdirs()
        def initBare = Git.init()
        initBare.setDirectory(barePath.toFile())
        initBare.setBare(true)
        def repoBare = initBare.call()
        repoBare.close()
        barePath.resolve('config').text = GIT_CONFIG_TEXT

        when:
        def manager = new AssetManager().build('test-org/test-repo')

        then:
        !manager.isNotInitialized()
        manager.getLocalGitConfig() == barePath.resolve('config').toFile()
        !manager.getLocalPath().exists()
        AssetManager.list() == ['test-org/test-repo']
    }

    def 'should detect repository status correctly'() {
        given:
        def folder = tempDir.getRoot()

        expect: 'no repository exists'
        new AssetManager().build('test-org/test-repo').isNotInitialized()


        when: 'only legacy exists'
        def legacyPath2 = folder.resolve('test-org/repo2')
        legacyPath2.mkdirs()
        legacyPath2.resolve('.git').mkdir()
        legacyPath2.resolve('.git/config').text = GIT_CONFIG_TEXT
        def manager = new AssetManager().build('test-org/repo2')
        then:
        manager.isNotInitialized()

        when: 'only bare exists'
        def barePath3 = folder.resolve(REPOS_SUBDIR + '/test-org/repo3/' + BARE_REPO)
        barePath3.mkdirs()
        def initBare3 = Git.init()
        initBare3.setDirectory(barePath3.toFile())
        initBare3.setBare(true)
        def repoBare3 = initBare3.call()
        repoBare3.close()
        barePath3.resolve('config').text = GIT_CONFIG_TEXT
        manager = new AssetManager().build('test-org/repo3')
        then:
        !manager.isNotInitialized()
        MultiRevisionRepository.checkProject(folder.toFile(), 'test-org/repo3')

        when: 'both exist'
        def legacyPath4 = folder.resolve('test-org/repo4')
        legacyPath4.mkdirs()
        legacyPath4.resolve('.git').mkdir()
        legacyPath4.resolve('.git/config').text = GIT_CONFIG_TEXT
        def barePath4 = folder.resolve(REPOS_SUBDIR + '/test-org/repo4/' + BARE_REPO)
        barePath4.mkdirs()
        def initBare4 = Git.init()
        initBare4.setDirectory(barePath4.toFile())
        initBare4.setBare(true)
        def repoBare4 = initBare4.call()
        repoBare4.close()
        barePath4.resolve('config').text = GIT_CONFIG_TEXT
        manager = new AssetManager().build('test-org/repo4')
        then:
        !manager.isNotInitialized()
        MultiRevisionRepository.checkProject(folder.toFile(), 'test-org/repo4')
    }

    // ============================================
    // DROP OPERATIONS TESTS
    // ============================================

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should drop a specific revision with multi-revision strategy'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def folder = tempDir.getRoot()
        def pipelineName = 'nextflow-io/hello'
        def revision1 = 'v1.2'
        def revision2 = 'mybranch'

        and: 'download two revisions'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download(revision1)
        manager.download(revision2)

        and: 'verify both revisions exist'
        def commit1Path = folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da')
        def commit2Path = folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64')
        commit1Path.exists()
        commit2Path.exists()

        when: 'drop the first revision'
        manager.drop(revision1)

        then: 'first revision is deleted but second remains'
        !commit1Path.exists()
        commit2Path.exists()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).exists()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should drop all revisions with multi-revision strategy'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def folder = tempDir.getRoot()
        def pipelineName = 'nextflow-io/hello'
        def revision1 = 'v1.2'
        def revision2 = 'mybranch'

        and: 'download two revisions'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download(revision1)
        manager.download(revision2)

        and: 'verify both revisions and bare repo exist'
        def projectPath = folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello')
        def commit1Path = projectPath.resolve(REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da')
        def commit2Path = projectPath.resolve(REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64')
        def barePath = projectPath.resolve(BARE_REPO)
        commit1Path.exists()
        commit2Path.exists()
        barePath.exists()

        when: 'drop all revisions (no revision specified)'
        manager.drop(null)

        then: 'everything is deleted including bare repo'
        !commit1Path.exists()
        !commit2Path.exists()
        !barePath.exists()
        !projectPath.exists()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should not drop revision with uncommitted changes unless forced'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def folder = tempDir.getRoot()
        def pipelineName = 'nextflow-io/hello'
        def revision = 'v1.2'

        and: 'download a revision'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download(revision)

        and: 'make local changes'
        def commitPath = folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da')
        commitPath.resolve('test-file.txt').text = 'uncommitted change'

        when: 'try to drop without force'
        manager.drop(revision, false)

        then: 'operation fails'
        thrown(AbortOperationException)
        commitPath.exists()

        when: 'drop with force flag'
        manager.drop(revision, true)

        then: 'revision is deleted'
        !commitPath.exists()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should handle drop of non-existent revision gracefully'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def pipelineName = 'nextflow-io/hello'

        and: 'create manager downloading a revision'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download("v1.2")

        when: 'try to drop a revision that was never downloaded'
        manager.drop('nonexistent-revision')

        then: 'no exception is thrown'
        noExceptionThrown()
    }

    // ============================================
    // REVISION SWITCHING TESTS
    // ============================================

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should switch between downloaded revisions'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def folder = tempDir.getRoot()
        def pipelineName = 'nextflow-io/hello'
        def revision1 = 'v1.2'
        def revision2 = 'mybranch'

        and: 'download two revisions'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download(revision1)
        manager.download(revision2)

        when: 'checkout first revision'
        manager.checkout(revision1)

        then: 'local path points to first revision'
        manager.getLocalPath().toString().contains('1b420d060d3fad67027154ac48e3bdea06f058da')
        manager.getCurrentRevision() == revision1

        when: 'checkout second revision'
        manager.checkout(revision2)

        then: 'local path points to second revision'
        manager.getLocalPath().toString().contains('1c3e9e7404127514d69369cd87f8036830f5cf64')
        manager.getCurrentRevision() == revision2
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should use setRevision to switch between commits'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def pipelineName = 'nextflow-io/hello'
        def revision1 = 'v1.2'
        def revision2 = 'mybranch'

        and: 'create manager and download revisions'
        def manager = new AssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.download(revision1)
        manager.download(revision2)

        when: 'use setRevision'
        manager.setRevision(revision1)

        then: 'revision is updated'
        manager.getLocalPath().toString().contains('1b420d060d3fad67027154ac48e3bdea06f058da')

        when: 'switch to another revision'
        manager.setRevision(revision2)

        then: 'local path updated'
        manager.getLocalPath().toString().contains('1c3e9e7404127514d69369cd87f8036830f5cf64')
    }


}
