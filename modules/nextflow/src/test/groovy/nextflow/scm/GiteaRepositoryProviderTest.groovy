/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
                server = 'https://git.mydomain.com'
                endpoint = 'https://git.mydomain.com/api/v1'
                platform = 'gitea'
                user = 'myname'
                password = 'mypassword'
                token = '65eaa9c8ef52460d22a93307fe0aee76289dc675'
            }

        }
        '''

    def 'should return repo url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('gitea', config.providers.mygitea as ConfigObject)

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getEndpointUrl() == 'https://git.mydomain.com/api/v1/repos/pditommaso/hello'
    }

    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('gitea', config.providers.mygitea as ConfigObject)

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getRepositoryUrl() == 'https://git.mydomain.com/pditommaso/hello'

    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('gitea', config.providers.mygitea as ConfigObject)

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getContentUrl('main.nf') == 'https://git.mydomain.com/api/v1/repos/pditommaso/hello/raw/main.nf'

    }

}
