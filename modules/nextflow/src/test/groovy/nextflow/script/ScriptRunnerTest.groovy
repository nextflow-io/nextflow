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

package nextflow.script

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.config.ConfigParser
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Timeout
import test.TestParser
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ScriptRunnerTest extends Specification {

    def 'test process' () {

        setup:
        def runner = new TestScriptRunner([process:[executor:'nope']])

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
        runner.result instanceof DataflowQueue
        runner.result.val == "echo Hello world"

    }


    def 'test processor config'() {

        /*
         * test that the *machineType* attribute is visible in the taskConfig object
         */
        when:
        def runner = new TestScriptRunner( process: [executor:'nope', machineType:'alpha'] )
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
        runner.getScriptObj().getTaskProcessor().name == 'simpleTask'
        runner.getScriptObj().getTaskProcessor().config.machineType == 'alpha'
    }


    def 'test process with args' () {
        setup:
        def runner = new TestScriptRunner( executor: 'nope' )

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
        runner.getScriptObj().getTaskProcessor().config.getInputs()[0].inChannel.getVal() == 1
        runner.getScriptObj().getTaskProcessor().config.getInputs()[1].inChannel instanceof DataflowQueue
        runner.getScriptObj().getTaskProcessor().config.getOutputs()[0].outChannel instanceof DataflowWriteChannel
    }


    def 'test process echo' () {

        setup:
        def runner = new TestScriptRunner( executor: 'nope' )

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
        runner.scriptObj.taskProcessor.name == 'test'

    }


    def 'test process variables' () {

        setup:
        def runner = new TestScriptRunner( executor: 'nope' )

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
        def runner = new TestScriptRunner( executor: 'nope' )

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
            @Override
            void abort(Throwable cause) {
                forceTermination()
            }
        }

        def runner = new TestScriptRunner(session)

        def script = '''
            process test {
                script:
                "$HELLO"
            }

            '''

        when:
        runner.setScript(script).execute()
        then:
        session.fault.error instanceof ProcessUnrecoverableException
        session.fault.error.cause instanceof MissingPropertyException
        session.fault.error.cause.message =~ /Unknown variable 'HELLO' -- .*/
        session.fault.report =~ /Unknown variable 'HELLO' -- .*/

    }


    def 'test process fallback variable' () {

        setup:
        def runner = new TestScriptRunner( executor: 'nope', env: [HELLO: 'Hello world!'] )

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
        def runner = new TestScriptRunner( executor: 'nope' )

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

        given:
        // -- this represent the configuration file
        def config = '''
            executor = 'nope'
            process {
                memory = '333'
                withName: hola { cpus = '222'; time = '555' }
                withName: ciao { cpus = '999' }
            }
            '''

        def script = '''
            process hola {
              penv 1
              cpus 2

              return ''
            }
            '''

        def session = new Session( new ConfigParser().parse(config))

        when:
        def process = new TestParser(session).parseAndGetProcess(script)

        then:
        process.config instanceof ProcessConfig
        process.config.penv == 1
        process.config.cpus == '222'  // !! this value is overridden by the one in the config file
        process.config.memory == '333'
        process.config.time == '555'

    }

    def 'test process name options 2'( ) {

        given:
        // -- this represent the configuration file
        def config = '''
            executor = 'nope'

            process {
                memory = '333'

                withName: hola {
                    cpus = '222'
                    time = '555'
                }

                withName: ciao {
                    cpus = '999'
                }
            }
            '''

        def script = '''
            process hola {
              penv 1
              cpus 2

              return ''
            }
            '''

        def session = new Session( new ConfigParser().parse(config))

        when:
        def process = new TestParser(session).parseAndGetProcess(script)

        then:
        process.config instanceof ProcessConfig
        process.config.penv == 1
        process.config.cpus == '222'  // !! this value is overridden by the one in the config file
        process.config.memory == '333'
        process.config.time == '555'

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

        def session = new Session(new ConfigParser().parse(config))

        when:
        def process = new TestParser(session).parseAndGetProcess(script)

        then:
        process.config instanceof ProcessConfig
        process.config.module == ['b/2','c/3']
        process.config.createTaskConfig().module == ['b/2','c/3']
        process.config.createTaskConfig().getModule() == ['b/2','c/3']
    }

    def 'test module config 2'() {

        setup:
        /*
         * the module defined in the config file 'b/2' has priority and overrides the 'a/1' and 'c/3'
         */
        def config = '''
            executor = 'nope'
            process {
                module = 'a/1'
                withName: hola { module = 'b/2:z/9' }
            }
            '''

        def script = '''
            process hola {
              module 'c/3'
              module 'd/4'

              'echo 1'
            }
            '''

        def session = new Session(new ConfigParser().parse(config))

        when:
        def process = new TestParser(session).parseAndGetProcess(script)

        then:
        process.config instanceof ProcessConfig
        process.config.module == ['b/2:z/9']
        process.config.createTaskConfig().module == ['b/2','z/9']
        process.config.createTaskConfig().getModule() == ['b/2','z/9']
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

        def session = new Session(new ConfigParser().parse(config))

        when:
        def process = new TestParser(session).parseAndGetProcess(script)

        then:
        process.config instanceof ProcessConfig
        process.config.module == ['a/1']
        process.config.createTaskConfig().module ==  ['a/1']
        process.config.createTaskConfig().getModule() ==  ['a/1']

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

        def config = new ConfigParser().parse(configText)


        when:
        def session = new Session(config)
        def process = new TestParser(session).parseAndGetProcess(script)
        then:
        process.config instanceof ProcessConfig
        process.config.queue == 'short'
        process.config.cpus == 2
        process.config.penv == 'mpi'
        process.config.memory == '10G'
        process.config.time == '6 hour'

        when:
        def result = new TestScriptRunner(config)
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
        def process = new TestParser(session).parseAndGetProcess(script)
        then:
        process.config instanceof ProcessConfig
        process.config.cpus == null

        when:
        def result = new TestScriptRunner(process: [executor:'nope'])
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

    def 'should parse mem and duration units' () {
        given:
        def script = '''  
            def result = [:] 
            result.mem1 = 1.GB
            result.mem2 = 1_000_000.toMemory()
            result.mem3 = MemoryUnit.of(2_000)
            result.time1 = 2.hours
            result.time2 = 60_000.toDuration()
            result.time3 = Duration.of(120_000)
            result.flag = 10000 < 1.GB 
            result // return result object
           '''

        when:
        def runner = new TestScriptRunner([executor:'nope'])
        def result = (Map)runner.setScript(script).execute()
        then:
        result.mem1 instanceof MemoryUnit
        result.mem1 == MemoryUnit.of('1 GB')
        result.mem2 == MemoryUnit.of(1_000_000)
        result.mem3 == MemoryUnit.of(2_000)
        result.time1 instanceof Duration
        result.time1 == Duration.of('2 hours')
        result.time2 == Duration.of(60_000)
        result.time3 == Duration.of(120_000)
        result.flag == true
    }

    def 'should define directive with a negative value' () {

        when:
        def script = '''
            X = 10
            process taskHello {
                maxRetries -1
                maxErrors -X
                ''
            }
            '''
        def runner = new TestScriptRunner([executor:'nope'])
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor

        then:
        processor.config.maxRetries == -1
        processor.config.maxErrors == -10

    }


    def 'should not thrown duplicate channel exception' () {

        when:
        def script = '''

            process foo {
              output: file '*.pdf'
              'touch x.pdf'
            }

            process bar {
              output: file '*.pdf'
              'touch x.pdf'
            }
                        '''
        def runner = new TestScriptRunner([executor:'nope'])
        runner.setScript(script).execute()

        then:
        noExceptionThrown()

    }


}
