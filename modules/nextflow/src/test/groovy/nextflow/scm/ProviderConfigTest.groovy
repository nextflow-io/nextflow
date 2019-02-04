/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProviderConfigTest extends Specification {

    static final String CONFIG = '''
        providers {
              github {
                user = '12732'
                password = '35454'
                server = 'https://github.com'
              }

              custom {
                server = 'http://local.host'
                platform = 'gitlab'
              }

              local {
                path = '/home/data/projects'
              }
        }
        '''

    def 'should create a ProviderConfig object' () {

        when:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj1 = new ProviderConfig('github', config.providers.github as ConfigObject)
        def obj2 = new ProviderConfig('custom', config.providers.custom as ConfigObject)
        def obj3 = new ProviderConfig('local', config.providers.local as ConfigObject )

        then:
        obj1.name == 'github'
        obj1.server == 'https://github.com'
        obj1.auth == '12732:35454'
        obj1.domain == 'github.com'
        obj1.platform == 'github'
        obj1.endpoint == 'https://api.github.com'

        obj2.name == 'custom'
        obj2.server == 'http://local.host'
        obj2.platform == 'gitlab'
        obj2.endpoint == 'http://local.host'

        obj3.platform == 'file'
        obj3.path == '/home/data/projects'
    }


    def 'should return defaults attributes' () {

        when:
        def config = new ProviderConfig('github')
        then:
        config.name == 'github'
        config.server == 'https://github.com'
        config.endpoint == 'https://api.github.com'
        config.platform == 'github'
        config.domain == 'github.com'

        when:
        config = new ProviderConfig('gitlab')
        then:
        config.name == 'gitlab'
        config.server == 'https://gitlab.com'
        config.endpoint == 'https://gitlab.com'
        config.platform == 'gitlab'
        config.domain == 'gitlab.com'

        when:
        config = new ProviderConfig('bitbucket')
        then:
        config.name == 'bitbucket'
        config.server == 'https://bitbucket.org'
        config.endpoint == 'https://bitbucket.org'
        config.platform == 'bitbucket'
        config.domain == 'bitbucket.org'
    }

    def 'should return provider attributes' () {

        when:
        def config = new ProviderConfig('custom',[platform: 'github', server:'http://local.host'])
        then:
        config.name == 'custom'
        config.platform == 'github'
        config.domain == 'local.host'
        config.server == 'http://local.host'
        config.endpoint == 'http://local.host'

        when:
        config = new ProviderConfig('my-git',[platform: 'gitea', server:'http://my-domain.org/gitea/'])
        then:
        config.name == 'my-git'
        config.platform == 'gitea'
        config.domain == 'my-domain.org'
        config.server == 'http://my-domain.org/gitea'
        config.endpoint == 'http://my-domain.org/gitea'

        when:
        config = new ProviderConfig('github')
        config.setPassword('abc')
        then:
        config.password == 'abc'
        config.auth == 'abc'
        config.authObfuscated == '-:***'
        config.user == null
        config.password == 'abc'

        when:
        config = new ProviderConfig('github')
        config.setUser('yo').setPassword('123')
        then:
        config.auth == 'yo:123'
        config.authObfuscated == 'yo:***'
        config.user == 'yo'
        config.password == '123'

        when:
        config = new ProviderConfig('github')
        config.setUser('hello')
        then:
        config.auth == 'hello:'
        config.authObfuscated == 'hello:-'
        config.user == 'hello'
        config.password == null

        when:
        config = new ProviderConfig('github', [auth: 'yo:123'])
        then:
        config.auth == 'yo:123'
        config.user == 'yo'
        config.password == '123'

        when:
        config = new ProviderConfig('github', [auth: 'xyz'])
        then:
        config.auth == null
        config.user == null
        config.password == null
        config.token == 'xyz'
    }

    def 'should ending slash and add protocol prefix' () {
        when:
        def config = new ProviderConfig('github',[server:'host.com///'])
        then:
        config.server == 'https://host.com'

    }

    def 'should create all config providers'() {

        when:
        def result = ProviderConfig.createFromText(CONFIG)
        then:
        result.size() == 6

        result.find { it.name == 'github' }.server == 'https://github.com'
        result.find { it.name == 'github' }.auth == '12732:35454'

        result.find { it.name == 'gitlab' }.server == 'https://gitlab.com'

        result.find { it.name == 'bitbucket' }.server == 'https://bitbucket.org'

        result.find { it.name == 'custom' }.server == 'http://local.host'
        result.find { it.name == 'custom' }.platform == 'gitlab'

        result.find { it.name == 'local' }.path == '/home/data/projects'
        result.find { it.name == 'local' }.platform == 'file'

        !result.find { it.name == 'xxxx' }
    }


}
