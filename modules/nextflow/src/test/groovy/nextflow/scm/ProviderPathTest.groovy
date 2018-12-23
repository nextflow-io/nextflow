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

import nextflow.config.ConfigBuilder
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProviderPathTest extends Specification {


    def 'should resolve with a github remote provider path' () {

        given:
        def provider = Mock(RepositoryProvider)
        def path = new ProviderPath(provider, 'nextflow.config')

        def MAIN_CONFIG = '''
            params.foo = 1 
            
            includeConfig 'conf/nested.config'
            '''

        def NESTED_CONFIG = '''
            params.bar = 2 
            process {
                cpus = 3 
                memory = '4 GB'
            }
            '''
        
        when:
        def cfg = new ConfigBuilder().buildGivenFiles(path)
        then:
        provider.readText('nextflow.config') >> MAIN_CONFIG
        provider.readText('conf/nested.config') >> NESTED_CONFIG
        cfg.params.foo == 1
        cfg.params.bar == 2
        cfg.process.cpus == 3
        cfg.process.memory == '4 GB'

    }

    def 'should return uri path' () {

        given:
        def provider = Mock(RepositoryProvider)
        def path = new ProviderPath(provider, 'conf/nextflow.config')
        def URI = 'http://github.com/foo/bar/content/conf/nextflow.config'
        when:
        def str = path.toUriString()
        then:
        provider.getContentUrl('conf/nextflow.config') >> URI
        str == URI
    }

}
