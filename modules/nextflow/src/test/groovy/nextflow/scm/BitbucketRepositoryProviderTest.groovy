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
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Timeout

@Timeout(10)
@IgnoreIf({System.getenv('NXF_SMOKE')})
class BitbucketRepositoryProviderTest extends Specification {

    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def testBitbucketCloneURL() {

        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)

        when:
        def url = new BitbucketRepositoryProvider('pditommaso/tutorial',config).getCloneUrl()
        then:
        url == "https://${config.user}@bitbucket.org/pditommaso/tutorial.git".toString()
    }


    def testGetHomePage() {
        expect:
        new BitbucketRepositoryProvider('pditommaso/tutorial').getRepositoryUrl() == "https://bitbucket.org/pditommaso/tutorial"
    }


    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def testReadContent() {

        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)

        when:
        def repo = new BitbucketRepositoryProvider('pditommaso/tutorial', config)
        def result = repo.readText('main.nf')

        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

    }

    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def 'should list tags' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)
        and:
        def repo = new BitbucketRepositoryProvider('pditommaso/tutorial', config)

        when:
        def result = repo.getTags()
        then:
        result.contains( new RepositoryProvider.TagInfo('v1.1', '755ba829cbc4f28dcb3c16b9dcc1c49c7ee47ff5') )
        result.contains( new RepositoryProvider.TagInfo('v1.0', '755ba829cbc4f28dcb3c16b9dcc1c49c7ee47ff5') )
    }

    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def 'should list branches' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)
        and:
        def repo = new BitbucketRepositoryProvider('pditommaso/tutorial', config)

        when:
        def result = repo.getBranches()
        then:
        result.contains( new RepositoryProvider.BranchInfo('test-branch', '755ba829cbc4f28dcb3c16b9dcc1c49c7ee47ff5') )
    }

    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def 'should return content URL' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)

        expect:
        new BitbucketRepositoryProvider('pditommaso/tutorial', config)
                .setRevision('test-branch')
                .getContentUrl('main.nf') == 'https://bitbucket.org/api/2.0/repositories/pditommaso/tutorial/src/test-branch/main.nf'

        and:
        new BitbucketRepositoryProvider('pditommaso/tutorial', config)
            .setRevision('feature/with-slash')
            .getContentUrl('main.nf') == 'https://bitbucket.org/api/2.0/repositories/pditommaso/tutorial/src/a6b825b22d46758cdeb496ae6cf26aef839ace52/main.nf'

        and:
        new BitbucketRepositoryProvider('pditommaso/tutorial', config)
            .setRevision('test/tag/v2')
            .getContentUrl('main.nf') == 'https://bitbucket.org/api/2.0/repositories/pditommaso/tutorial/src/8f849beceb2ea479ef836809ca33d3daeeed25f9/main.nf'

    }

    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def 'should read content on main and on separate branch' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)
        and:
        def repo = new BitbucketRepositoryProvider('pditommaso/tutorial', config)

        expect:
        repo.readText('main.nf').contains('world')
        !repo.readText('main.nf').contains('WORLD')

        and:
        repo.setRevision('test-branch')
        repo.readText('main.nf').contains('world')
        !repo.readText('main.nf').contains('WORLD')

        and:
        repo.setRevision('feature/with-slash')
        !repo.readText('main.nf').contains('world')
        repo.readText('main.nf').contains('WORLD')

        and:
        repo.setRevision('v1.1')
        repo.readText('main.nf').contains('world')
        !repo.readText('main.nf').contains('WORLD')

        and:
        repo.setRevision('test/tag/v2')
        !repo.readText('main.nf').contains('world')
        repo.readText('main.nf').contains('mundo')

    }
}
