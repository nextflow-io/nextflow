/*
 * Copyright (c) 2012, the authors.
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
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import org.apache.commons.io.FileUtils
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CliRunnerTest extends Specification {

 //TODO check why it stops testing with gradle

//    def 'test task' () {
//
//        setup:
//        def runner = new CliRunner([task:[processor:'nope']])
//
//        /*
//         * Test a task with a very simple body.
//         * For testing purposes the processor just return the script itself as result
//         */
//        when:
//        def script =
//            """
//            task('task1')  {
//                "echo Hello world"
//            }
//            """
//
//        def result = runner.execute(script)
//
//        then:
//        result == "echo Hello world"
//
//        cleanup:
//        runner.workDirectory?.deleteDir()
//
//    }


    def 'test cli params' () {
        setup:
        def runner = new CliRunner()

        when:
        def script =
            '''
            params['p1'] = 1
            params['p2'] = 2
            params['p3'] = 3

            return params
            '''
        runner.setParam( ['p1':'10','p2':'Hola'] )
        def result = runner.execute( script )

        then:
        result.p1 == 10
        result.p2 == 'Hola'
        result.p3 == 3


    }

    def 'test task with assignment' () {
        setup:
        def runner = new CliRunner([task:[processor:'nope']])

        when:
        def script =
            '''
            task('task2')  {

                input x: 1
                input y: [3,4]

                def z = 'Hello world!'

                """echo $z"""
            }

            '''
        def result = runner.execute(script)

        then:
        result == 'echo Hello world!'
        runner.getScript().getTaskProcessor().getName() == 'task2'
        runner.getScript().getTaskProcessor().getInput('x').getVal() == 1
        runner.getScript().getTaskProcessor().getInput('y') instanceof DataflowQueue
        runner.getScript().getTaskProcessor().getOutput('-') instanceof DataflowWriteChannel

        cleanup:
        runner.workDirectory?.deleteDir()

    }

    def 'test task with syntax error' () {
        setup:
        def runner = new CliRunner([task:[processor:'nope']])

        /*
         * this declaration returns a syntax error because the task code block
         * does not terminate with a script to execute
         */
        when:
        def script =
            """
            task('task1')  {
                input x: 'hola'
            }
            """
        runner.execute(script)

        then:
        thrown(MultipleCompilationErrorsException)

        cleanup:
        runner.workDirectory?.deleteDir()

    }

    def 'test task variables' () {


        setup:
        def runner = new CliRunner([task:[processor:'nope']])

        def script = '''
            X = 1
            Y = 2
            task {
                input Y
                def Z = 3

                """$X-$Y-$Z"""
            }

            '''


        expect:
        runner.execute(script) == '1-2-3'

        cleanup:
        runner.workDirectory?.deleteDir()

    }

    def 'test version' () {

        when:
        def opt = CliRunner.parseMainArgs('-v')

        then:
        assert opt.version

    }

    def 'test help' () {

        when:
        def opt = CliRunner.parseMainArgs('-h')

        then:
        assert opt.help

    }


    def 'test usage' () {

        when:
        CliRunner.parseMainArgs([] as String)
        CliRunner.jcommander.usage()

        then:
        noExceptionThrown()

    }

    def 'validateWorkDir' () {

        /*
         * the user specified a path that does not exist
         * it will created
         */
        when:
        def path = CliRunner.validateWorkDirectory( 'x1_testValidateFolder', 'scriptName' )

        then:
        path == new File('x1_testValidateFolder').canonicalFile
        path.exists()

        cleanup:
        if( path ) FileUtils.deleteDirectory(path)

    }

    def 'validateWorkDir 2' () {

        setup:
        def folder = new File('x2_testValidateFolder')
        folder.mkdirs()
        def file = File.createTempFile('test','test',folder)

        /*
         * The folder specified by the user contains a file
         * -> it is not a valid folder, an exception is thrown
         */
        when:
        def path = CliRunner.validateWorkDirectory( folder.toString(), 'scriptName' )

        then:
        thrown(CliRunner.CliArgumentException.class)

        cleanup:
        if( folder ) FileUtils.deleteDirectory(folder)

    }

    def 'validateWorkDir 3' () {

        /*
         * No folder is specified, create a temp folder in the current dir, using the script name
         */
        when:
        def path = CliRunner.validateWorkDirectory( null, 'myScript' )

        then:
        path == new File('run-myScript').canonicalFile
        path.exists()

        cleanup:
        path?.delete()

    }

    def 'buildConfigObject' () {

        setup:
        def env = [PATH:'/local/bin', HOME:'/home/my']

        when:
        def config = CliRunner.buildConfig0(env,null)

        then:
        ('PATH' in config.env )
        ('HOME' in config.env )
        !('XXX' in config.env )

        config.env.PATH == '/local/bin'
        config.env.HOME == '/home/my'

    }

    def 'buildConfigObject 2 ' () {

        setup:
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        '''

        def text2 = '''
        task { field2 = 'Hello' }
        env { beta = 'b2'; delta = 'd2'; HOME="$HOME:/local/path"; XXX="$PATH:/xxx"; YYY = "$XXX:/yyy"   }
        '''

        when:
        def config1 = CliRunner.buildConfig0(env, [text1])
        def config2 = CliRunner.buildConfig0(env, [text1, text2])

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

    }


    def 'test validateConfigFiles '() {
        when:
        def files = CliRunner.validateConfigFiles(['file1','file2'])

        then:
        thrown(CliRunner.CliArgumentException)

        when:
        def f1 = File.createTempFile('file1','x').absoluteFile
        def f2 = File.createTempFile('file1','x').absoluteFile
        files = CliRunner.validateConfigFiles([f1.toString(), f2.toString()])

        then:
        files == [f1, f2]

        cleanup:
        f1?.delete()
        f2?.delete()

    }


    def 'test config to map ' () {

        setup:
        def config = new ConfigSlurper().parse( 'task {field1=1; field2="two"}; env { x = 99 } ' )

        when:
        def map = CliRunner.configToMap(config)
        map.env.PATH = '/local/bin'

        then:
        map.task.field1 == 1
        map.task.field2 == "two"
        map.env.x == 99
        map.env.y == null
        map.env.PATH == '/local/bin'

    }


    def 'test parseValue' () {

        expect:
        CliRunner.parseValue(str) == value

        where:
        str         | value
        'hola'      | 'hola'
        '1'         | 1
        "${Long.MAX_VALUE}" | Long.MAX_VALUE
        'True'      | true
        'False'     | false
        "10.2"      | 10.2

    }

    def 'test setParams' () {
        when:
        def cli = new CliRunner()
        cli.setParam('x', '1')
        cli.setParam('y', 'true')
        cli.setParam('z', 'Hola')

        then:
        cli.params.x == 1
        cli.params.y == true
        cli.params.z == 'Hola'




    }



}
