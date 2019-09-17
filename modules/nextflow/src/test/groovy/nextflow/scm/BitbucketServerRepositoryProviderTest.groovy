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

    def 'should return repo url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj).getEndpointUrl() == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello'
    }

    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)
        def testString = new BitbucketServerRepositoryProvider('pditommaso/hello', obj).getRepositoryUrl()

        expect:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj).getRepositoryUrl() == 'https://bitbucket.server.com/scm/pditommaso/hello'

    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)
        def testString = new GiteaRepositoryProvider('pditommaso/hello', obj).getContentUrl('main.nf')

        expect:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj).getContentUrl('main.nf') == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello/raw/main.nf'

    }

}
