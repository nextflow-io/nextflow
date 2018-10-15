/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
