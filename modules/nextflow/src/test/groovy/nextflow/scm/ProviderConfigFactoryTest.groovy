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

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProviderConfigFactoryTest extends Specification {

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


    def 'should create config providers from text'() {

        when:
        def result = ProviderConfigFactory.createFromText(CONFIG)
        then:
        result.size() == 7

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
