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

/**
 *
 * @author Akira Sekiguchi <pachiras.yokohama@gmail.com>
 */
class GiteaRepositoryProviderTest extends Specification {

    static final String CONFIG = '''
        providers {

            mygitea {
                server = 'https://git.seqera.io'
                endpoint = 'https://git.seqera.io/api/v1'
                platform = 'gitea'
                user = 'myname'
                password = 'mypassword'
            }

        }
        '''

    def 'should return repo url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('gitea', config.providers.mygitea as ConfigObject)

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getEndpointUrl() == 'https://git.seqera.io/api/v1/repos/pditommaso/hello'
    }

    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('gitea', config.providers.mygitea as ConfigObject)

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getRepositoryUrl() == 'https://git.seqera.io/pditommaso/hello'

    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('gitea', config.providers.mygitea as ConfigObject)

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('main.nf') == 'https://git.seqera.io/api/v1/repos/pditommaso/hello/raw/main.nf'
        and:
        new GiteaRepositoryProvider('pditommaso/hello', obj)
                .setRevision('12345')
                .getContentUrl('main.nf') == 'https://git.seqera.io/api/v1/repos/pditommaso/hello/raw/main.nf?ref=12345'

    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_GITEA_ACCESS_TOKEN')})
    def 'should read file content'() {

        given:
        def token = System.getenv('NXF_GITEA_ACCESS_TOKEN')
        def config = new ProviderConfig('gitea', new ConfigSlurper().parse(CONFIG).providers.mygitea as ConfigObject).setAuth(token)

        when:
        def repo = new GiteaRepositoryProvider('test-org/nextflow-ci-repo', config)
        def result = repo.readText('README.md')
        then:
        result.contains('nextflow-ci-repo')

//        when:
//        repo = new GiteaRepositoryProvider('test-org/nextflow-ci-repo', config)
//                        .setRevision('foo')
//        result = repo.readText('README.md')
//        then:
//        result.contains("foo branch")

    }

}
