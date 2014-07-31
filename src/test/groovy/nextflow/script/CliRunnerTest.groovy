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
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.exception.MissingLibraryException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import spock.lang.Specification
import test.TestParser
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CliRunnerTest extends Specification {

    def testProcess () {

        setup:
        def runner = new CliRunner([process:[executor:'nope']])

        /*
         * Test a task with a very simple body.
         * For testing purposes the processor just return the script itself as result
         */
        when:
        def script =
            """
            process sayHello  {
                "echo Hello world"
            }
            """

        runner.execute(script)

        // when no outputs are specified, the 'stdout' is the default output
        then:
        runner.result instanceof DataflowVariable
        runner.result.val == "echo Hello world"

    }


    def testProcessorConfig() {

        /*
         * test that the *instanceType* attribute is visible in the taskConfig object
         */
        when:
        def runner = new CliRunner( process: [executor:'nope', instanceType:'alpha'] )
        def script =
            '''
            process simpleTask  {
                input:
                val x from 1
                output:
                stdout into result

                """echo $x"""
            }

            '''
        runner.execute(script)
        then:
        runner.getScript().getTaskProcessor().taskConfig.name == 'simpleTask'
        runner.getScript().getTaskProcessor().taskConfig.instanceType == 'alpha'


        /*
         * test that the *instanceType* property defined by the task (beta)
         * override the one define in the main config (alpha)
         */
        when:
        def runner2 = new CliRunner( process: [executor:'nope', instanceType:'alpha'] )
        def script2 =
            '''
            process otherTask  {
                instanceType 'beta'
                input:
                val x from 1
                output:
                stdout result

                """echo $x"""
            }

            '''
        runner2.execute(script2)

        then:
        runner2.getScript().getTaskProcessor().taskConfig.name == 'otherTask'
        runner2.getScript().getTaskProcessor().taskConfig.instanceType == 'beta'

    }


    def testProcessWithArgs () {
        setup:
        def runner = new CliRunner( executor: 'nope' )

        when:
        def script =
            '''
            process task2  {
                input:
                val x from 1
                val y from ([3])
                output:
                stdout result

                """echo $x - $y"""
            }

            '''
        runner.execute(script)

        then:
        runner.getResult().val == 'echo 1 - 3'
        runner.getScript().getTaskProcessor().getName() == 'task2'
        runner.getScript().getTaskProcessor().taskConfig.name == 'task2'
        runner.getScript().getTaskProcessor().taskConfig.inputs[0].inChannel.getVal() == 1
        runner.getScript().getTaskProcessor().taskConfig.inputs[1].inChannel instanceof DataflowQueue
        runner.getScript().getTaskProcessor().taskConfig.outputs[0].outChannel instanceof DataflowWriteChannel
    }


    def testProcessEcho () {

        setup:
        def runner = new CliRunner( executor: 'nope' )

        when:
        def script =
            '''
            process test  {
                input:
                val x from 1
                output:
                stdout result

                "echo $x"
            }
            '''
        runner.execute(script)

        then:
        runner.getResult().val == 'echo 1'
        runner.script.taskProcessor.taskConfig.name == 'test'

    }



    def testProcessVariables () {


        setup:
        def runner = new CliRunner( executor: 'nope' )

        def script = '''
            X = 1
            Y = 2
            process test {
                input:
                val Y

                "$X-$Y-3"
            }

            '''


        expect:
        runner.execute(script).val == '1-2-3'

    }

    def testProcessVariables2 () {

        setup:
        def runner = new CliRunner( executor: 'nope' )

        def script = '''
            X = 1
            Y = 2
            process test {
                input:
                val Y

                exec:
                def Z = 3
                "$X-$Y-$Z"
            }

            '''

        expect:
        runner.execute(script).val == '1-2-3'

    }


    def testProcessOutFile () {


        setup:
        def runner = new CliRunner( executor: 'nope' )

        def script = '''
            X = file('filename')
            process test {
                input:
                file X

                "cat $X"
            }

            '''

        expect:
        runner.execute(script).val == 'cat filename'

    }

    def testVersion () {

        when:
        def opt = CliRunner.parseMainArgs('-v')

        then:
        assert opt.version

    }

    def testHelp () {

        when:
        def opt = CliRunner.parseMainArgs('-h')

        then:
        assert opt.help

    }


    def testUsage () {

        when:
        CliRunner.parseMainArgs([] as String)
        CliRunner.jcommander.usage()

        then:
        noExceptionThrown()

    }



    def buildConfigObject () {

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
        def config1 = CliRunner.buildConfig0(env, [text1])

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
        def files = CliRunner.validateConfigFiles(['file1','file2'])

        then:
        thrown(CliArgumentException)

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


    def testConfigToMap  () {

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
        runner.libraries.sort() == [ path2, jar1, jar2 ].sort()


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

    def testProcessNameOptions ( ) {

        setup:
        // -- this represent the configuration file
        def config = '''
            executor = 'nope'

            process.delta = '333'
            process.$hola.beta = '222'
            process.$hola.gamma = '555'

            process.$ciao.beta = '999'

            '''

        def script = '''
            process hola {
              alpha 1
              beta 2

              input:
              val x

              return ''
            }
            '''

        def session = new Session( new ConfigSlurper().parse(config))

        when:
        TaskProcessor process = new TestParser(session).parseScript(script).run()

        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.alpha == 1
        process.taskConfig.beta == '222'  // !! this value is overridden by the one in the config file
        process.taskConfig.delta == '333'
        process.taskConfig.gamma == '555'

    }

    def testProcessNameOptions2( ) {

        setup:
        // -- this represent the configuration file
        def config = '''
            executor = 'nope'

            process {
                delta = '333'

                $hola {
                    beta = '222'
                    gamma = '555'
                }

                $ciao {
                    beta = '999'
                }
            }


            '''

        def script = '''
            process hola {
              alpha 1
              beta 2

              input:
              val x

              return ''
            }
            '''

        def session = new Session( new ConfigSlurper().parse(config))

        when:
        TaskProcessor process = new TestParser(session).parseScript(script).run()

        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.alpha == 1
        process.taskConfig.beta == '222'  // !! this value is overridden by the one in the config file
        process.taskConfig.delta == '333'
        process.taskConfig.gamma == '555'

    }

    def testModulesConfig() {

        setup:
        // -- this represent the configuration file
        def config = '''
            executor = 'nope'
            process.module = 'a/1'
            '''

        def script = '''
            process hola {
              module 'b/2'
              module 'c/3'

              'echo 1'
            }
            '''

        def session = new Session(new ConfigSlurper().parse(config))

        when:
        TaskProcessor process = new TestParser(session).parseScript(script).run()

        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.getModule() == ['b/2','c/3']
    }

    def testModulesConfig2() {

        setup:
        /*
         * the module defined in the config file 'b/2' has priority and overrides the 'a/1' and 'c/3'
         */
        def config = '''
            executor = 'nope'
            process.module = 'a/1'
            process.$hola.module = 'b/2:z/9'
            '''

        def script = '''
            process hola {
              module 'c/3'
              module 'd/4'

              'echo 1'
            }
            '''

        def session = new Session(new ConfigSlurper().parse(config))

        when:
        TaskProcessor process = new TestParser(session).parseScript(script).run()

        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.getModule() == ['b/2','z/9']
    }

    def testModulesConfig3() {

        setup:
        /*
         * the module defined in the config file 'b/2' has priority and overrides the 'a/1' and 'c/3'
         */
        def config = '''
            executor = 'nope'
            process.module = 'a/1'
            '''

        def script = '''
            process hola {
              'echo 1'
            }
            '''

        def session = new Session(new ConfigSlurper().parse(config))

        when:
        TaskProcessor process = new TestParser(session).parseScript(script).run()

        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.getModule() == ['a/1']
    }


    def testCommandLineOptions() {

        when:
        def opt = CliRunner.parseMainArgs('-daemon.x=1', '-daemon.y.z=2')


        then:
        opt.daemonOptions.x == '1'
        opt.daemonOptions.'y.z'== '2'

    }

    def testCommandDaemonOptions() {

        when:
        def opt = new CliOptions(daemonOptions: [group:'pippo', join:'192.168.0.1', 'x.y.z': 123])
        def result = CliRunner.buildConfig([], opt)
        then:
        result.daemon == [group:'pippo', join:'192.168.0.1', x:[y:[z:123]]]

    }

    def testCommandExecutorOptions() {

        when:
        def opt = new CliOptions(executorOptions: [ alpha: 1, 'beta.x': 'hola', 'beta.y': 'ciao' ])
        def result = CliRunner.buildConfig([], opt)
        then:
        result.executor.alpha == 1
        result.executor.beta.x == 'hola'
        result.executor.beta.y == 'ciao'


    }

//
//    def testParamsOverride() {
//
//        given:
//        def config = '''
//            params.cpus = 10
//
//            derived.value = "x ${params.cpus}"
//            '''
//
//        when:
//        def obj = new ConfigSlurper().parse(config)
//        then:
//        obj.params.cpus == 10
//        obj.derived.value == 'x 10'
//
//
//        when:
//        def map = [params:[cpus:20] ]
//        def parser = new ConfigSlurper()
//        parser.setBinding(map)
//        obj = parser.parse(config)
//        then:
//        obj.params.cpus == 20
//        obj.derived.value == 'x 20'
//    }


}
