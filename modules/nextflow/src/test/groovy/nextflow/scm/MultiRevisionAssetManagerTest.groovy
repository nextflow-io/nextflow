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

import static MultiRevisionAssetManager.BARE_REPO
import static MultiRevisionAssetManager.REVISION_SUBDIR

import spock.lang.IgnoreIf

import org.junit.Rule
import spock.lang.Requires
import spock.lang.Specification
import test.TemporaryPath

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class MultiRevisionAssetManagerTest extends Specification {

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def setup() {
        MultiRevisionAssetManager.root = tempDir.root.toFile()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list revisions and commits'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def pipelineName ='nextflow-io/hello'
        def revision1 = 'v1.2'
        def revision2 = 'mybranch'
        def manager = new MultiRevisionAssetManager().build(pipelineName, [providers: [github: [auth: token]]])
        manager.setRevision(revision1)
        manager.createSharedClone()
        manager.setRevision(revision2)
        manager.createSharedClone()

        when:
        def list = manager.listRevisions()
        then:
        list.size() == 2
        list.contains(revision1)
        list.contains(revision2)

        when:
        def branchesAndTags = manager.getBranchesAndTags(false)
        then:
        branchesAndTags.size() == 2
        branchesAndTags[revision1] == '1b420d060d3fad67027154ac48e3bdea06f058da'
        branchesAndTags[revision2] == '1c3e9e7404127514d69369cd87f8036830f5cf64'
    }

    def 'should list commits' () {
        given:
        def folder = tempDir.getRoot()

        when:
        def manager = new MultiRevisionAssetManager('cbcrg/pipe1')
        folder.resolve('cbcrg/pipe1/' + REVISION_SUBDIR + '/12345').mkdirs()
        folder.resolve('cbcrg/pipe1/' + REVISION_SUBDIR + '/67890').mkdirs()
        def list = manager.listCommits()
        then:
        list.sort() == ['12345','67890']

        when:
        manager = new MultiRevisionAssetManager('cbcrg/pipe2')
        folder.resolve('cbcrg/pipe2/' + REVISION_SUBDIR + '/abcde').mkdirs()
        folder.resolve('cbcrg/pipe2/' + REVISION_SUBDIR + '/fghij').mkdirs()
        list = manager.listCommits()
        then:
        list.sort() == ['abcde','fghij']
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should clone bare repo'() {

        given:
        def folder = tempDir.getRoot()
        String revision = null
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new MultiRevisionAssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevision(revision)

        when:
        manager.checkBareRepo()
        then:
        folder.resolve('nextflow-io/hello/' + BARE_REPO).isDirectory()
        folder.resolve('nextflow-io/hello/' + BARE_REPO + '/config').exists()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should get revision form commit using bare repo' () {

        given:
        def folder = tempDir.getRoot()
        String revision = null
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new MultiRevisionAssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevision(revision)

        when:
        manager.checkBareRepo()
        then:
        manager.revisionToCommitWithBareRepo('v1.2') == '1b420d060d3fad67027154ac48e3bdea06f058da'
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should create shared clone from commit'() {

        given:
        def folder = tempDir.getRoot()
        String revision = '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' // easier to fix commit for a generic test
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new MultiRevisionAssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevision(revision)

        when:
        manager.createSharedClone()
        then:
        def git = folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' + '/.git')
        git.isDirectory()
        git.resolve('objects/info/alternates').text == folder.resolve('nextflow-io/hello/'+ BARE_REPO + '/objects').toAbsolutePath().toString()

        when:
        def result = manager.createSharedClone()
        then:
        noExceptionThrown()
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should create shared clone from tag'() {

        given:
        def folder = tempDir.getRoot()
        String revision = 'v1.2'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new MultiRevisionAssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevision( revision)
        // tag v1.2 -> commit 1b420d060d3fad67027154ac48e3bdea06f058da

        when:
        manager.createSharedClone()
        then:
        def git = folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '1b420d060d3fad67027154ac48e3bdea06f058da' + '/.git')
        git.isDirectory()
        git.resolve('objects/info/alternates').text == folder.resolve('nextflow-io/hello/'+ BARE_REPO + '/objects').toAbsolutePath().toString()

        when:
        manager.createSharedClone()
        then:
        noExceptionThrown()
    }


    // Downloading a branch first and then pulling the branch
    // should work fine, unlike with tags.
    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should create shared clone from branch'() {

        given:
        def folder = tempDir.getRoot()
        String revision = 'mybranch'
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new MultiRevisionAssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.setRevision(revision)
        // as of Jun 2024, branch "mybranch" -> commit "1c3e9e7404127514d69369cd87f8036830f5cf64"

        when:
        manager.createSharedClone()
        then:
        def git = folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '1c3e9e7404127514d69369cd87f8036830f5cf64' + '/.git').isDirectory()

        when:
        def result = manager.createSharedClone()
        then:
        noExceptionThrown()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should identify default branch when pulling the repo'() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new MultiRevisionAssetManager().build('nextflow-io/socks', [providers: [github: [auth: token]]])

        when:
        // simulate calling `nextflow run nextflow-io/socks` without specifying a revision
        manager.createSharedClone()
        then:
        manager.getCurrentRevision() == 'main'

    }

}
