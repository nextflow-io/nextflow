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
 * @author Tobias Neumann <tobias.neumann.at@gmail.com>
 */
class AzureDevOpsRepositoryProviderTest extends Specification {

    static final String CONFIG = '''
        providers {

            azurerepos {
                server = 'https://dev.azure.com'
                endpoint = 'https://dev.azure.com'
                platform = 'azurerepos'
                user = 'myname'
                password = 'mypassword'
                token = '65eaa9c8ef52460d22a93307fe0aee76289dc675'
            }

        }
        '''

    def 'should return repo url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureDevOpsRepositoryProvider('t-neumann/hello', obj).getEndpointUrl() == 'https://dev.azure.com/t-neumann/hello/_apis/git/repositories/hello'
    }

    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureDevOpsRepositoryProvider('t-neumann/hello', obj).getRepositoryUrl() == 'https://dev.azure.com/t-neumann/hello'

    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureDevOpsRepositoryProvider('t-neumann/hello', obj).getContentUrl('main.nf') == 'https://dev.azure.com/t-neumann/hello/_apis/git/repositories/hello/items?download=false&includeContent=true&includeContentMetadata=false&api-version=6.0&\$format=json&path=main.nf'

    }

}
