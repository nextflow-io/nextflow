/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
                auth = '12732:35454'
                host = 'https://github.com'
              }

              custom {
                host = 'http://local.host'
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
        obj1.host == 'https://github.com'
        obj1.auth == '12732:35454'
        obj1.domain == 'github.com'
        obj1.platform == 'github'
        obj1.endpoint == 'https://api.github.com'

        obj2.name == 'custom'
        obj2.host == 'http://local.host'
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
        config.host == 'https://github.com'
        config.endpoint == 'https://api.github.com'
        config.platform == 'github'
        config.domain == 'github.com'

        when:
        config = new ProviderConfig('gitlab')
        then:
        config.name == 'gitlab'
        config.host == 'https://gitlab.com'
        config.endpoint == 'https://gitlab.com'
        config.platform == 'gitlab'
        config.domain == 'gitlab.com'

        when:
        config = new ProviderConfig('bitbucket')
        then:
        config.name == 'bitbucket'
        config.host == 'https://bitbucket.org'
        config.endpoint == 'https://bitbucket.org'
        config.platform == 'bitbucket'
        config.domain == 'bitbucket.org'
    }

    def 'should return provider attributes' () {

        when:
        def config = new ProviderConfig('custom',[platform: 'github', auth: 'xyz', host:'http://local.host'])
        then:
        config.name == 'custom'
        config.platform == 'github'
        config.domain == 'local.host'
        config.host == 'http://local.host'
        config.endpoint == 'http://local.host'
        config.auth == 'xyz'

        when:
        config.setAuth('abc')
        then:
        config.auth == 'abc'
        config.authObfuscated == '***'
        config.user == null
        config.password == 'abc'

        when:
        config.setAuth('yo','123')
        then:
        config.auth == 'yo:123'
        config.authObfuscated == 'yo:***'
        config.user == 'yo'
        config.password == '123'

        when:
        config.setAuth(null,'123')
        then:
        config.auth == ':123'
        config.authObfuscated == '-:***'
        config.user == null
        config.password == '123'

        when:
        config.setAuth('hello',null)
        then:
        config.auth == 'hello:'
        config.authObfuscated == 'hello:-'
        config.user == 'hello'
        config.password == null
    }

    def 'should create all config providers'() {

        given:
        when:
        def result = ProviderConfig.createFromText(CONFIG)
        then:
        result.size() == 5

        result.find { it.name == 'github' }.host == 'https://github.com'
        result.find { it.name == 'github' }.auth == '12732:35454'

        result.find { it.name == 'gitlab' }.host == 'https://gitlab.com'

        result.find { it.name == 'bitbucket' }.host == 'https://bitbucket.org'

        result.find { it.name == 'custom' }.host == 'http://local.host'
        result.find { it.name == 'custom' }.platform == 'gitlab'

        result.find { it.name == 'local' }.path == '/home/data/projects'
        result.find { it.name == 'local' }.platform == 'file'

        !result.find { it.name == 'xxxx' }
    }


}
