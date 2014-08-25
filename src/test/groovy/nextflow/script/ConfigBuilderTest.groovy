/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.script

import nextflow.cli.CliOptions
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigBuilderTest extends Specification {


    def buildConfigObject () {

        setup:
        def env = [PATH:'/local/bin', HOME:'/home/my']

        when:
        def config = ConfigBuilder.buildConfig0(env,null)

        then:
        ('PATH' in config.env )
        ('HOME' in config.env )
        !('XXX' in config.env )

        config.env.PATH == '/local/bin'
        config.env.HOME == '/home/my'

    }

    def buildConfigObject2 () {

        setup:
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        '''

        def text2 = '''
        task { field2 = 'Hello' }
        env { beta = 'b2'; delta = 'd2'; HOME="$HOME:/local/path"; XXX="$PATH:/xxx"; YYY = "$XXX:/yyy"; WWW = "${WWW?:''}:value"   }
        '''

        when:
        def config1 = ConfigBuilder.buildConfig0(env, [text1])
        def config2 = ConfigBuilder.buildConfig0(env, [text1, text2])

        // note: configuration object can be modified like any map
        config2.env ['ZZZ'] = '99'

        then:
        config1.task.field1 == 1
        config1.task.field2 == 'hola'
        config1.env.HOME == '/home/my:/some/path'
        config1.env.PATH == '/local/bin'
        config1.env.'dot.key.name' == 'any text'

        config2.task.field1 == 1
        config2.task.field2 == 'Hello'
        config2.env.alpha == 'a1'
        config2.env.beta == 'b2'
        config2.env.delta == 'd2'
        config2.env.HOME == '/home/my:/local/path'
        config2.env.XXX == '/local/bin:/xxx'
        config2.env.PATH == '/local/bin'
        config2.env.YYY == '/local/bin:/xxx:/yyy'
        config2.env.ZZZ == '99'
        config2.env.WWW == ':value'

    }

    def buildConfigObject3 () {

        setup:
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        params { demo = 1   }
        '''


        when:
        def config1 = ConfigBuilder.buildConfig0(env, [text1])

        then:
        config1.task.field1 == 1
        config1.task.field2 == 'hola'
        config1.env.HOME == '/home/my:/some/path'
        config1.env.PATH == '/local/bin'
        config1.env.'dot.key.name' == 'any text'

        config1.params.demo == 1


    }


    def testValidateConfigFiles () {
        when:
        def files = new ConfigBuilder().validateConfigFiles(['file1','file2'])
        then:
        thrown(AbortOperationException)

        when:
        def f1 = File.createTempFile('file1','x').absoluteFile
        def f2 = File.createTempFile('file1','x').absoluteFile
        files = new ConfigBuilder().validateConfigFiles([f1.toString(), f2.toString()])
        then:
        files == [f1, f2]

        cleanup:
        f1?.delete()
        f2?.delete()

    }


    def testConfigToMap  () {

        setup:
        def config = new ConfigSlurper().parse( 'task {field1=1; field2="two"}; env { x = 99 } ' )

        when:
        def map = ConfigBuilder.configToMap(config)
        map.env.PATH = '/local/bin'

        then:
        map.task.field1 == 1
        map.task.field2 == "two"
        map.env.x == 99
        map.env.y == null
        map.env.PATH == '/local/bin'

    }

    def testCommandDaemonOptions() {

        when:
        def opt = new CliOptions(daemonOptions: [group:'pippo', join:'192.168.0.1', 'x.y.z': 123])
        def result = new ConfigBuilder().setOptions(opt).buildConfig([])
        then:
        result.daemon == [group:'pippo', join:'192.168.0.1', x:[y:[z:123]]]

    }

    def testCommandExecutorOptions() {

        when:
        def opt = new CliOptions()
        def run = new CmdRun(executorOptions: [ alpha: 1, 'beta.x': 'hola', 'beta.y': 'ciao' ])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([])
        then:
        result.executor.alpha == 1
        result.executor.beta.x == 'hola'
        result.executor.beta.y == 'ciao'

    }

}
