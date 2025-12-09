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

import nextflow.config.Manifest

import static MultiRevisionRepositoryStrategy.BARE_REPO
import static MultiRevisionRepositoryStrategy.REVISION_SUBDIR

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
class MultiRevisionRepositoryStrategyTest extends Specification {

    @Rule
    TemporaryPath tempDir = new TemporaryPath()


    private MultiRevisionRepositoryStrategy createStrategy (String project, String token){
        def provider = token ? new GithubRepositoryProvider(project, new ProviderConfig('github').setAuth(token)): null
        def context = new RepositoryContext(project, tempDir.root.toFile(), tempDir.root.resolve(project).toFile(), provider)
        return new MultiRevisionRepositoryStrategy(context)
    }
    def 'should list commits' () {
        given:
        def folder = tempDir.getRoot()

        when:
        def strategy = createStrategy('cbcrg/pipe1', null)
        folder.resolve('cbcrg/pipe1/' + REVISION_SUBDIR + '/12345').mkdirs()
        folder.resolve('cbcrg/pipe1/' + REVISION_SUBDIR + '/67890').mkdirs()
        def list = strategy.listDownloadedCommits()
        then:
        list.sort() == ['12345','67890']

        when:
        def context2 = new RepositoryContext('cbcrg/pipe2', folder.toFile(), folder.resolve('cbcrg/pipe2').toFile(), null)
        strategy = new MultiRevisionRepositoryStrategy(context2)
        folder.resolve('cbcrg/pipe2/' + REVISION_SUBDIR + '/abcde').mkdirs()
        folder.resolve('cbcrg/pipe2/' + REVISION_SUBDIR + '/fghij').mkdirs()
        list = strategy.listDownloadedCommits()
        then:
        list.sort() == ['abcde','fghij']
    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should clone bare repo and get revisions'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def strategy = createStrategy('nextflow-io/hello', token)
        def manifest = Mock(Manifest){
            getDefaultBranch() >> 'master'
            getRecurseSubmodules() >> false
        }

        when:
        strategy.checkBareRepo(manifest)
        then:
        folder.resolve('nextflow-io/hello/' + BARE_REPO).isDirectory()
        folder.resolve('nextflow-io/hello/' + BARE_REPO + '/config').exists()

        expect:
        strategy.revisionToCommitWithBareRepo('v1.2') == '1b420d060d3fad67027154ac48e3bdea06f058da'

    }



    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should create shared clone from commit, tag and branch'() {

        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def strategy = createStrategy('nextflow-io/hello', token)
        def manifest = Mock(Manifest){
            getDefaultBranch() >> 'master'
            getRecurseSubmodules() >> false
        }

        when:
        strategy.download('7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8', 1, manifest)
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8/.git').isDirectory()
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8/.git/objects/info/alternates').text == folder.resolve('nextflow-io/hello/'+ BARE_REPO + '/objects').toAbsolutePath().toString()



        when:
        // tag v1.2 -> commit 1b420d060d3fad67027154ac48e3bdea06f058da
        strategy.download('v1.2', 1, manifest)
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git').isDirectory()
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/1b420d060d3fad67027154ac48e3bdea06f058da/.git/objects/info/alternates').text == folder.resolve('nextflow-io/hello/'+ BARE_REPO + '/objects').toAbsolutePath().toString()

        when:
        strategy.download('mybranch', 1, manifest)
        then:
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git').isDirectory()
        folder.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/1c3e9e7404127514d69369cd87f8036830f5cf64/.git/objects/info/alternates').text == folder.resolve('nextflow-io/hello/'+ BARE_REPO + '/objects').toAbsolutePath().toString()
    }

}
