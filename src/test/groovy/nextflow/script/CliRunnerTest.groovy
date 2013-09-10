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
import nextflow.exception.MissingLibraryException
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



    def 'test task with assignment' () {
        setup:
        def runner = new CliRunner( executor: 'nope' )

        when:
        def script =
            '''
            process task2  {
                input:
                val x using 1
                val y using ([3])
                output:
                stdout result

                """echo $x - $y"""
            }

            '''
        def result = runner.execute(script)

        then:
        result == 'echo 1 - 3'
        runner.getScript().getTaskProcessor().getName() == 'task2'
        runner.getScript().getTaskProcessor().taskConfig.name == 'task2'
        runner.getScript().getTaskProcessor().taskConfig.inputs[0].channel.getVal() == 1
        runner.getScript().getTaskProcessor().taskConfig.inputs[1].channel instanceof DataflowQueue
        runner.getScript().getTaskProcessor().taskConfig.outputs[0].channel instanceof DataflowWriteChannel
    }


    def 'test task echo' () {

        setup:
        def runner = new CliRunner( executor: 'nope' )

        when:
        def script =
            '''
            process test  {
                input:
                val x using 1
                output:
                stdout result

                "echo $x"
            }
            '''
        def result = runner.execute(script)

        then:
        result == 'echo 1'
        runner.script.taskProcessor.taskConfig.name == 'test'

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
            process task1  {
                input x: 'hola'
            }
            """
        runner.execute(script)

        then:
        thrown(MultipleCompilationErrorsException)

    }

    def 'test task variables' () {


        setup:
        def runner = new CliRunner( executor: 'nope' )

        def script = '''
            X = 1
            Y = 2
            process test {
                input Y
                def Z = 3

                "$X-$Y-$Z"
            }

            '''


        expect:
        runner.execute(script) == '1-2-3'

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

    def 'buildConfigObject 3 ' () {

        setup:
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        params { demo = 1   }
        '''


        when:
        def config1 = CliRunner.buildConfig0(env, [text1])

        then:
        config1.task.field1 == 1
        config1.task.field2 == 'hola'
        config1.env.HOME == '/home/my:/some/path'
        config1.env.PATH == '/local/bin'
        config1.env.'dot.key.name' == 'any text'

        config1.params.demo == 1


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


    def 'test normalize cmdline' () {

        expect:
        CliRunner.normalizeArgs('a','-bb','-ccc','dddd') == ['a','-bb','-ccc','dddd']
        CliRunner.normalizeArgs('a','-bb','-ccc','-resume', 'last') == ['a','-bb','-ccc','-resume','last']
        CliRunner.normalizeArgs('a','-bb','-ccc','-resume') == ['a','-bb','-ccc','-resume','last']
        CliRunner.normalizeArgs('a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz') == ['a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz']

        CliRunner.normalizeArgs('x','-test') == ['x','-test','%all']
        CliRunner.normalizeArgs('x','-test','alpha') == ['x','-test','alpha']
        CliRunner.normalizeArgs('x','-test','-other') == ['x','-test','%all','-other']

        CliRunner.normalizeArgs('--alpha=1') == ['--alpha=1']
        CliRunner.normalizeArgs('--alpha','1') == ['--alpha=1']
        CliRunner.normalizeArgs('-x', '1', 'script.nf', '--long', 'v1', '--more', 'v2', '--flag') == ['-x','1','script.nf','--long=v1','--more=v2','--flag=true']

        CliRunner.normalizeArgs('-x', '1', '-process.alpha','2', '3') == ['-x', '1', '-process.alpha=2', '3']
        CliRunner.normalizeArgs('-x', '1', '-process.echo') == ['-x', '1', '-process.echo=true']


        CliRunner.normalizeArgs('-x', '1', '-that.alpha','2', '3') == ['-x', '1', '-that.alpha','2', '3']
    }



    def testAddLibPath() {

        setup:
        def path1 = File.createTempDir()
        def path2 = File.createTempDir()
        def jar1 = new File(path2, 'lib1.jar'); jar1.createNewFile()
        def jar2 = new File(path2, 'lib2.jar'); jar2.createNewFile()

        when:
        def runner = new CliRunner()
        runner.addLibPaths( path1 )
        then:
        runner.libraries == [ path1 ]


        when:
        runner = new CliRunner()
        runner.addLibPaths( path2 )
        then:
        runner.libraries == [ path2, jar1, jar2 ]


        cleanup:
        path1.deleteDir()
        path2.deleteDir()

    }


    def testSetLibPath() {


        setup:
        def path1 = File.createTempDir()
        def path2 = File.createTempDir()
        def jar1 = new File(path2, 'lib1.jar'); jar1.createNewFile()
        def jar2 = new File(path2, 'lib2.jar'); jar2.createNewFile()

        when:
        def runner = new CliRunner()
        runner.setLibPath(jar1.toString())
        then:
        runner.libraries == [ jar1 ]

        when:
        runner = new CliRunner()
        runner.libPath = "${jar1}:${jar2}"
        then:
        runner.libraries == [ jar1, jar2 ]


        when:
        runner = new CliRunner()
        runner.setLibPath('some/lib.jar')

        then:
        thrown(MissingLibraryException)

        cleanup:
        path1.deleteDir()
        path2.deleteDir()



    }


}
