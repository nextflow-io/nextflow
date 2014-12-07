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
import nextflow.exception.ProcessScriptException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import test.TestParser
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptRunnerTest extends Specification {

    def 'test process' () {

        setup:
        def runner = new ScriptRunner([process:[executor:'nope']])

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

        runner.setScript(script).execute()

        // when no outputs are specified, the 'stdout' is the default output
        then:
        runner.result instanceof DataflowVariable
        runner.result.val == "echo Hello world"

    }


    def 'test processor config'() {

        /*
         * test that the *instanceType* attribute is visible in the taskConfig object
         */
        when:
        def runner = new ScriptRunner( process: [executor:'nope', instanceType:'alpha'] )
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
        runner.setScript(script).execute()
        then:
        runner.getScriptObj().getTaskProcessor().taskConfig.name == 'simpleTask'
        runner.getScriptObj().getTaskProcessor().taskConfig.instanceType == 'alpha'


        /*
         * test that the *instanceType* property defined by the task (beta)
         * override the one define in the main config (alpha)
         */
        when:
        def runner2 = new ScriptRunner( process: [executor:'nope', instanceType:'alpha'] )
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
        runner2.setScript(script2).execute()

        then:
        runner2.getScriptObj().getTaskProcessor().taskConfig.name == 'otherTask'
        runner2.getScriptObj().getTaskProcessor().taskConfig.instanceType == 'beta'

    }


    def 'test process with args' () {
        setup:
        def runner = new ScriptRunner( executor: 'nope' )

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
        runner.setScript(script).execute()

        then:
        runner.getResult().val == 'echo 1 - 3'
        runner.getScriptObj().getTaskProcessor().getName() == 'task2'
        runner.getScriptObj().getTaskProcessor().taskConfig.name == 'task2'
        runner.getScriptObj().getTaskProcessor().taskConfig.inputs[0].inChannel.getVal() == 1
        runner.getScriptObj().getTaskProcessor().taskConfig.inputs[1].inChannel instanceof DataflowQueue
        runner.getScriptObj().getTaskProcessor().taskConfig.outputs[0].outChannel instanceof DataflowWriteChannel
    }


    def 'test process echo' () {

        setup:
        def runner = new ScriptRunner( executor: 'nope' )

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
        runner.setScript(script).execute()

        then:
        runner.getResult().val == 'echo 1'
        runner.scriptObj.taskProcessor.taskConfig.name == 'test'

    }



    def 'test process variables' () {

        setup:
        def runner = new ScriptRunner( executor: 'nope' )

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
        runner.setScript(script).execute().val == '1-2-3'

    }

    def 'test process variables 2' () {

        setup:
        def runner = new ScriptRunner( executor: 'nope' )

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
        runner.setScript(script).execute().val == '1-2-3'

    }


    def 'test process missing variable' () {

        setup:
        def session = new Session( executor: 'nope' ) {
            Throwable error
            @Override
            void abort(Throwable cause) {
                this.error = cause
                forceTermination()
            }
        }

        def runner = new ScriptRunner(session)

        def script = '''
            process test {
                script:
                "$HELLO"
            }

            '''

        when:
        runner.setScript(script).execute()
        then:
        session.error instanceof ProcessScriptException
        session.error.cause instanceof MissingPropertyException
        session.error.cause.message =~ /Unknown variable 'HELLO' -- .*/

    }


    def 'test process fallback variable' () {

        setup:
        def runner = new ScriptRunner( executor: 'nope', env: [HELLO: 'Hello world!'] )

        def script = '''
            process test {
                exec:
                "$HELLO"
            }

            '''

        expect:
        runner.setScript(script).execute().val == 'Hello world!'

    }



    def 'test process output file' () {


        setup:
        def runner = new ScriptRunner( executor: 'nope' )

        def script = '''
            X = file('filename')
            process test {
                input:
                file X

                "cat $X"
            }

            '''

        expect:
        runner.setScript(script).execute().val == 'cat filename'

    }


    def 'test process name options' ( ) {

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

    def 'test process name options 2'( ) {

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

    def 'test module config'() {

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

    def 'test module config 2'() {

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

    def 'test module config 3'() {

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


    def 'test resource'() {

        setup:
        // -- this represent the configuration file
        def configText = '''
            executor = 'nope'
            process {
              queue = 'short'
              cpus  = 2
              time  = '6 hour'
              penv  = 'mpi'
              memory = '10G'
            }
            '''

        def script = '''
            process hola {
              """
              queue: ${task.queue}
              cpus: ${task.cpus}
              time: ${task.time}
              penv: ${task.penv}
              nodes: ${task.nodes}
              memory: ${task.memory}
              """
            }
            '''

        def config = new ConfigSlurper().parse(configText)


        when:
        def session = new Session(config)
        TaskProcessor process = new TestParser(session).parseScript(script).run()
        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.queue == 'short'
        process.taskConfig.cpus == 2
        process.taskConfig.penv == 'mpi'
        process.taskConfig.memory == new MemoryUnit('10G')
        process.taskConfig.time == new Duration('6h')

        when:
        def result = new ScriptRunner(config)
                    .setScript(script)
                    .execute()
                    .getVal()
                    .toString()
                    .stripIndent()
                    .trim()
                    .readLines()

        then:
        result[0] == 'queue: short'
        result[1] == 'cpus: 2'
        result[2] == 'time: 6h'
        result[3] == 'penv: mpi'
        result[5] == 'memory: 10 GB'

    }


    def 'test resource with default'() {

        setup:
        def script = '''
            process hola {
              """
              cpus: ${task.cpus}
              """
            }
            '''

        when:
        def session = new Session()
        TaskProcessor process = new TestParser(session).parseScript(script).run()
        then:
        process.taskConfig instanceof TaskConfig
        process.taskConfig.cpus == null

        when:
        def result = new ScriptRunner(process: [executor:'nope'])
                .setScript(script)
                .execute()
                .getVal()
                .toString()
                .stripIndent()
                .trim()
                .readLines()

        then:
        result[0] == 'cpus: 1'

    }

}
