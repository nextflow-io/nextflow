/*
 * Copyright 2013-2025, Seqera Labs
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

import nextflow.config.Manifest
import org.eclipse.jgit.api.Git

import static MultiRevisionRepositoryStrategy.BARE_REPO
import static MultiRevisionRepositoryStrategy.REVISION_SUBDIR
import static nextflow.scm.MultiRevisionRepositoryStrategy.REPOS_SUBDIR

import spock.lang.IgnoreIf

import org.junit.Rule
import spock.lang.Requires
import spock.lang.Specification
import test.TemporaryPath

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class MultiRevisionRepositoryStrategyTest extends Specification {

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def setup() {
        AssetManager.root = tempDir.root.toFile()
    }

    private MultiRevisionRepositoryStrategy createStrategy(String project, String token) {
        final strategy = new MultiRevisionRepositoryStrategy(project)
        if( token )
            strategy.setProvider(new GithubRepositoryProvider(project, new ProviderConfig('github').setAuth(token)))
        return strategy
    }

    def 'should list commits'() {
        given:
        def folder = tempDir.getRoot()

        when:
        def strategy = createStrategy('cbcrg/pipe1', null)
        folder.resolve(REPOS_SUBDIR + '/cbcrg/pipe1/' + REVISION_SUBDIR + '/12345').mkdirs()
        folder.resolve(REPOS_SUBDIR + '/cbcrg/pipe1/' + REVISION_SUBDIR + '/67890').mkdirs()
        def list = strategy.listDownloadedCommits()
        then:
        list.sort() == ['12345', '67890']

        when:
        strategy = createStrategy('cbcrg/pipe2', null)
        folder.resolve(REPOS_SUBDIR + '/cbcrg/pipe2/' + REVISION_SUBDIR + '/abcde').mkdirs()
        folder.resolve(REPOS_SUBDIR + '/cbcrg/pipe2/' + REVISION_SUBDIR + '/fghij').mkdirs()
        list = strategy.listDownloadedCommits()
        then:
        list.sort() == ['abcde', 'fghij']
    }

    @IgnoreIf({ System.getenv('NXF_SMOKE') })
    @Requires({ System.getenv('NXF_GITHUB_ACCESS_TOKEN') })
    def 'should clone bare repo and get revisions'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def strategy = createStrategy('nextflow-io/hello', token)
        def manifest = Mock(Manifest) {
            getDefaultBranch() >> 'master'
            getRecurseSubmodules() >> false
        }

        when:
        strategy.checkBareRepo(manifest)
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO + '/config').exists()

        expect:
        strategy.revisionToCommitWithBareRepo('v1.2') == '1b420d060d3fad67027154ac48e3bdea06f058da'

    }

    @IgnoreIf({ System.getenv('NXF_SMOKE') })
    @Requires({ System.getenv('NXF_GITHUB_ACCESS_TOKEN') })
    def 'should create shared clone from commit, tag and branch'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def strategy = createStrategy('nextflow-io/hello', token)
        def manifest = Mock(Manifest) {
            getDefaultBranch() >> 'master'
            getRecurseSubmodules() >> false
        }

        when:
        strategy.download('7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8', 1, manifest)
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8/.git').isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8/.git/objects/info/alternates').text == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO + '/objects').toAbsolutePath().toString()

        when:
        // tag v1.2 -> commit 1b420d060d3fad67027154ac48e3bdea06f058da
        strategy.download('v1.2', 1, manifest)
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git').isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git/objects/info/alternates').text == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO + '/objects').toAbsolutePath().toString()

        when:
        strategy.download('mybranch', 1, manifest)
        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git').isDirectory()
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git/objects/info/alternates').text == folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO + '/objects').toAbsolutePath().toString()
    }

    def 'should create correct RefSpec for branches tags and commits'() {
        given:
        def folder = tempDir.getRoot()
        // Create a bare repo manually for testing
        def bareRepoPath = folder.resolve(REPOS_SUBDIR + '/test/project/' + BARE_REPO).toFile()
        bareRepoPath.mkdirs()
        def git = Git.init().setDirectory(bareRepoPath).setBare(true).call()
        // Configure origin to point to itself so lsRemote works (returns empty refs)
        def config = git.getRepository().getConfig()
        config.setString("remote", "origin", "url", bareRepoPath.absolutePath)
        config.save()
        git.close()

        def strategy = new MultiRevisionRepositoryStrategy('test/project')

        when:
        // Access the private method via Groovy's metaclass
        def refSpec = strategy.invokeMethod('refSpecForName', 'some-revision')

        then:
        // When revision is not found locally or remotely, it's treated as a commit
        refSpec.toString() == 'some-revision:refs/tags/some-revision'
    }

    @IgnoreIf({ System.getenv('NXF_SMOKE') })
    @Requires({ System.getenv('NXF_GITHUB_ACCESS_TOKEN') })
    def 'should fetch new remote branch not in local bare repo'() {
        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        // Use a test repo with known branches
        def strategy = createStrategy('nextflow-io/hello', token)
        def manifest = Mock(Manifest) {
            getDefaultBranch() >> 'master'
            getRecurseSubmodules() >> false
        }

        when:
        // First, create bare repo with default branch only
        strategy.checkBareRepo(manifest)

        then:
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + BARE_REPO).isDirectory()

        when:
        // Now try to download 'mybranch' which exists remotely but wasn't fetched initially
        // Clear the local branch ref to simulate it not existing locally
        def bareGit = Git.open(strategy.getBareRepo())
        def localBranchRef = bareGit.getRepository().findRef("refs/heads/mybranch")
        // If mybranch was fetched during clone, delete it to simulate fresh fetch scenario
        if (localBranchRef) {
            bareGit.branchDelete().setBranchNames("mybranch").setForce(true).call()
        }
        bareGit.close()

        // This should succeed by querying remote refs
        strategy.download('mybranch', 1, manifest)

        then:
        noExceptionThrown()
        // mybranch commit is 1c3e9e7404127514d69369cd87f8036830f5cf64
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git').isDirectory()

        when:
        // Now try to download tag 'v1.2' which exists remotely but wasn't fetched initially
        // Clear the local tag ref to simulate it not existing locally
        bareGit = Git.open(strategy.getBareRepo())
        def localTagRef = bareGit.getRepository().findRef("refs/tags/v1.2")
        // If v1.2 was fetched during clone, delete it to simulate fresh fetch scenario
        if (localTagRef) {
            bareGit.tagDelete().setTags("v1.2").call()
        }
        bareGit.close()

        // This should succeed by querying remote refs
        strategy.download('v1.2', 1, manifest)

        then:
        noExceptionThrown()
        // v1.2 ref is 1b420d060d3fad67027154ac48e3bdea06f058
        folder.resolve(REPOS_SUBDIR + '/nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git').isDirectory()
    }

}
