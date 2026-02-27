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

package nextflow.script

import groovyx.gpars.dataflow.DataflowVariable
import nextflow.exception.ScriptCompilationException
import nextflow.extension.Bolts
import nextflow.processor.TaskProcessor
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class ScriptRunnerTest extends Dsl2Spec {

    def 'test process' () {

        given:
        def config = [process:[executor:'nope']]
        /*
         * Test a task with a very simple body.
         * For testing purposes the processor just return the script itself as result
         */
        and:
        def script =
            """
            process sayHello  {
              output:
                stdout
              script:
                "echo Hello world"
            }

            workflow {
              sayHello()
            }
            """

        when:
        def result = runScript(script, config)

        then:
        result instanceof DataflowVariable
        result.val == "echo Hello world"

    }


    def 'test processor config'() {

        /*
         * test that the *machineType* attribute is visible in the taskConfig object
         */
        given:
        def config = [process: [executor:'nope', machineType:'alpha']]
        and:
        def script =
            '''
            process simpleTask  {
                input:
                val x
                output:
                stdout

                script:
                """echo $x"""
            }

            workflow {
                simpleTask(1)
            }
            '''
        when:
        runScript(script, config)
        def processor = TaskProcessor.currentProcessor()
        then:
        processor.name == 'simpleTask'
        processor.config.machineType == 'alpha'
    }


    def 'test process with args' () {
        given:
        def script =
            '''
            process simpleTask  {
                input:
                val x
                val y
                output:
                stdout

                script:
                """echo $x - $y"""
            }

            workflow {
                simpleTask(1, channel.of(3))
            }
            '''

        when:
        def result = runScript(script)

        then:
        result.val == 'echo 1 - 3'

    }


    def 'test process echo' () {

        given:
        def script =
            '''
            process simpleTask  {
                input:
                val x
                output:
                stdout

                script:
                "echo $x"
            }

            workflow {
                simpleTask(1)
            }
            '''

        when:
        def result = runScript(script)

        then:
        result.val == 'echo 1'
        TaskProcessor.currentProcessor().name == 'simpleTask'

    }


    def 'test process variables' () {

        given:
        def script = '''
            process simpleTask {
                input:
                val X
                val Y
                output:
                stdout

                script:
                "$X-$Y-3"
            }

            workflow {
                simpleTask(1, 2)
            }
            '''

        when:
        def result = runScript(script)
        then:
        result.val == '1-2-3'

    }

    def 'test process variables 2' () {

        given:
        def script = '''
            process simpleTask {
                input:
                val X
                val Y
                output:
                stdout

                script:
                def Z = 3
                "$X-$Y-$Z"
            }

            workflow {
                simpleTask(1, 2)
            }
            '''

        when:
        def result = runScript(script)
        then:
        result.val == '1-2-3'

    }

    def 'test process missing variable' () {

        given:
        def script = '''
            process simpleTask {
                script:
                "echo $HELLO"
            }

            workflow { simpleTask() }
            '''

        when:
        def config = [process:[executor: 'nope']]
        runScript(script, config)

        then:
        def e = thrown(ScriptCompilationException)
        and:
        e.cause.message.contains '`HELLO` is not defined'

    }


    def 'test process output file' () {
        given:
        def script = '''
            process simpleTask {
                input:
                file X
                output:
                stdout

                script:
                "cat $X"
            }

            workflow {
                X = file('filename')
                simpleTask(X)
            }
            '''
        and:
        def config = [process: [executor: 'nope']]

        when:
        def result = runScript(script, config)
        then:
        result.val == 'cat filename'

    }


    def 'test process name options' ( ) {

        given:
        def config = loadConfig('''
            process {
                executor = 'nope'
                memory = '333'
                withName: hola { cpus = '222'; time = '555' }
                withName: ciao { cpus = '999' }
            }
            ''')

        def script = '''
            process hola {
              penv 1
              cpus 2

              script:
              'echo hola'
            }

            workflow { hola() }
            '''

        when:
        runScript(script, config)
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.penv == 1
        process.config.cpus == '222'  // !! this value is overridden by the one in the config file
        process.config.memory == '333'
        process.config.time == '555'

    }

    def 'test process name options 2'( ) {

        given:
        def config = loadConfig('''
            process {
                executor = 'nope'
                memory = '333'

                withName: hola {
                    cpus = '222'
                    time = '555'
                }

                withName: ciao {
                    cpus = '999'
                }
            }
            ''')

        def script = '''
            process hola {
              penv 1
              cpus 2

              script:
              'echo hola'
            }

            workflow { hola() }
            '''

        when:
        runScript(script, config)
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.penv == 1
        process.config.cpus == '222'  // !! this value is overridden by the one in the config file
        process.config.memory == '333'
        process.config.time == '555'

    }

    def 'test module config'() {

        given:
        def config = loadConfig('''
            process.executor = 'nope'
            process.module = 'a/1'
            ''')

        def script = '''
            process hola {
              module 'b/2'
              module 'c/3'

              script:
              'echo 1'
            }

            workflow { hola() }
            '''

        when:
        runScript(script, config)
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.module == ['b/2','c/3']
        process.config.createTaskConfig().module == ['b/2','c/3']
        process.config.createTaskConfig().getModule() == ['b/2','c/3']
    }

    def 'test module config 2'() {
        given:
        /*
         * the module defined in the config file 'b/2' has priority and overrides the 'a/1' and 'c/3'
         */
        def config = loadConfig('''
            process {
                executor = 'nope'
                module = 'a/1'
                withName: hola { module = 'b/2:z/9' }
            }
            ''')

        def script = '''
            process hola {
              module 'c/3'
              module 'd/4'

              script:
              'echo 1'
            }

            workflow { hola() }
            '''

        when:
        runScript(script, config)
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.module == ['b/2:z/9']
        process.config.createTaskConfig().module == ['b/2','z/9']
        process.config.createTaskConfig().getModule() == ['b/2','z/9']
    }

    def 'test module config 3'() {
        given:
        /*
         * the module defined in the config file 'b/2' has priority and overrides the 'a/1' and 'c/3'
         */
        def config = loadConfig('''
            process.executor = 'nope'
            process.module = 'a/1'
            ''')

        def script = '''
            process hola {
              'echo 1'
            }

            workflow { hola() }
            '''

        when:
        runScript(script, config)
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.module == ['a/1']
        process.config.createTaskConfig().module ==  ['a/1']
        process.config.createTaskConfig().getModule() ==  ['a/1']

    }


    def 'test resource'() {
        given:
        def config = loadConfig('''
            process {
              executor = 'nope'
              queue = 'short'
              cpus  = 2
              time  = '6 hour'
              penv  = 'mpi'
              memory = '10G'
            }
            ''')

        def script = '''
            process hola {
              output: stdout
              script:
              """
              queue: ${task.queue}
              cpus: ${task.cpus}
              time: ${task.time}
              penv: ${task.penv}
              nodes: ${task.nodes}
              memory: ${task.memory}
              """
            }

            workflow {
              hola()
            }
            '''

        when:
        def result = runScript(script, config)
                .getVal()
                .toString()
                .stripIndent()
                .trim()
                .readLines()
        and:
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.queue == 'short'
        process.config.cpus == 2
        process.config.penv == 'mpi'
        process.config.memory == '10G'
        process.config.time == '6 hour'

        and:
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
              output: stdout
              script:
              """
              cpus: ${task.cpus}
              """
            }

            workflow {
              hola()
            }
            '''

        and:
        def config = [process: [executor:'nope']]

        when:
        def result = runScript(script, config)
                .getVal()
                .toString()
                .stripIndent()
                .trim()
                .readLines()
        and:
        def process = TaskProcessor.currentProcessor()

        then:
        process.config instanceof ProcessConfig
        process.config.cpus == null

        and:
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
        def result = runScript(script)
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

        given:
        def script = '''
            process taskHello {
                maxRetries -1
                maxErrors -X
                input:
                val X
                script:
                'echo hello'
            }

            workflow { taskHello(10) }
            '''

        when:
        def result = runScript(script)
        def processor = TaskProcessor.currentProcessor()

        then:
        processor.config.maxRetries == -1
        Bolts.resolveLazy([X: 10], processor.config.maxErrors) == -10

    }


    def 'test stub command'() {

        given:
        /*
         * the module defined in the config file 'b/2' has priority and overrides the 'a/1' and 'c/3'
         */
        def config = loadConfig('''
            process.executor = 'nope'
            stubRun = true
            ''')

        def script = '''
            process hola {
                output:
                stdout

                script:
                /echo bar/

                stub:
                /echo foo/
            }

            workflow {
                hola()
            }
            '''

        when:
        def result = runScript(script, config)

        then:
        result instanceof DataflowVariable
        result.val == "echo foo"

    }
}
