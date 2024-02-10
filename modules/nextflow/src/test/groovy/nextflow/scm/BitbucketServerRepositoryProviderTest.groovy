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

import spock.lang.Ignore
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Piotr Faba <piotr.faba@ardigen.com>
 */
class BitbucketServerRepositoryProviderTest extends Specification {

    static final String CONFIG = '''
        providers {

            bbserver {
                server = 'https://bitbucket.server.com'
                endpoint = 'https://bitbucket.server.com'
                platform = 'bitbucketserver'
                user = 'myname'
                password = 'mypassword'
            }

        }
        '''

    @Unroll
    def 'should return endpoint url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider(NAME, obj).getEndpointUrl() == ENDPOINT

        where:
        NAME                                    | ENDPOINT
        'pditommaso/hello'                      | 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello'
        'scm/pditommaso/hello'                  | 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello'
        'projects/DA/repos/crispr-pipeline'     | 'https://bitbucket.server.com/rest/api/1.0/projects/DA/repos/crispr-pipeline'
    }

    @Unroll
    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider(NAME, obj).getRepositoryUrl() == REPO_URL

        where:
        NAME                                | REPO_URL
//        'pditommaso/hello'                  | 'https://bitbucket.server.com/pditommaso/hello'
        '/pditommaso/hello'                 | 'https://bitbucket.server.com/pditommaso/hello'
//        and:
//        'scm/pditommaso/hello'              | 'https://bitbucket.server.com/scm/pditommaso/hello'
//        '/scm/pditommaso/hello'             | 'https://bitbucket.server.com/scm/pditommaso/hello'
//        and:
//        'projects/DA/repos/crispr-pipeline' | 'https://bitbucket.server.com/projects/DA/repos/crispr-pipeline'
//        '/projects/DA/repos/crispr-pipeline'| 'https://bitbucket.server.com/projects/DA/repos/crispr-pipeline'
    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('main.nf') == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello/raw/main.nf'

        and:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj)
                .setRevision('foo')
                .getContentUrl('main.nf') == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello/raw/main.nf?at=foo'

    }

    @Ignore
    @Requires( { System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN') } )
    def 'should list branches' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN')
        def config = new ProviderConfig('bbs', [server:'http://slurm.seqera.io:7990', platform:'bitbucketsever']).setAuth(token)
        and:
        def repo = new BitbucketServerRepositoryProvider('scm/hello/hello', config)

        when:
        def result = repo.getBranches()
        then:
        result.contains( new RepositoryProvider.BranchInfo('master', 'c62df3d9c2464adcaa0fb6c978c8e32e2672b191') )
    }

    @Ignore
    @Requires( { System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN') } )
    def 'should list tags' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN')
        def config = new ProviderConfig('bbs', [server:'http://slurm.seqera.io:7990', platform:'bitbucketsever']).setAuth(token)
        and:
        def repo = new BitbucketServerRepositoryProvider('scm/hello/hello', config)

        when:
        def result = repo.getTags()
        then:
        result.contains( new RepositoryProvider.TagInfo('v1.0', 'c62df3d9c2464adcaa0fb6c978c8e32e2672b191') )
    }
}
