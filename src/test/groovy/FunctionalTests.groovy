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
import nextflow.processor.TaskProcessor
import nextflow.script.ScriptRunner
import nextflow.util.MemoryUnit
import spock.lang.Shared
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FunctionalTests extends Specification {

    @Shared
    File scriptFile

    // run before the first feature method
    def setupSpec() {
        scriptFile = new File('test.nf')
        scriptFile.deleteOnExit()
    }

    // run after the last feature method
    def cleanupSpec() {
    }

    /**
     * test passing values through environment variables
     */
    def 'test environment' () {

        setup:
        def environment = [ XXX: 'value1']
        def script = '''

            def x = config.env['XXX']
            def y = config.env['YYY'] ?: -1

            [ x, y ]

            '''

        when:
        def runner = new ScriptRunner( [env: environment] )

        then:
        runner.setScript(script).execute() == ['value1', -1]

    }


    /*
     * test passing values through command line argument (unnamed parameters)
     */
    def 'test args'()  {

        when:
        def script = """
            def len = args.size()
            def x = args[0]
            def y = args[1]

            return [ x, y, len ]
            """
        def runner = new ScriptRunner()
        def result = runner.setScript(script).execute(['hello', 'hola'] )

        then:
        result[0] == 'hello'
        result[1] == 'hola'
        result[2] == 2

    }


    def 'test configure processor'() {

        setup:
        def configStr = '''
             process {
                dummyField = 99
                executor = 'nope'
                echo = true
                shell = 'zsh'
                maxForks = 10
                environment = [a:1, b:2,c:3]
                validExitStatus = [0,11,22,33]
            }
            '''
        def cfg = new ConfigSlurper().parse(configStr)


        when:
        def script = '''

            process taskHello {
                echo true
                maxForks 11
                ''
            }

            '''

        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.getName() == 'taskHello'
        processor.config.echo == true
        processor.config.shell == 'zsh'
        processor.config.maxForks == 11
        processor.config.dummyField == 99
        processor.config.environment.entrySet() == [a:1,b:2,c:3].entrySet()
        processor.config.validExitStatus == [0,11,22,33]

    }

    def 'test merge ext properties' () {
        given:
        def configStr = '''
            process {
                ext {
                    alpha = "aaa"
                    delta = "ddd"
                }
                $foo {
                    ext {
                        beta = "BBB"
                        delta = "DDD"
                    }
                }
            }
        '''
        def cfg = new ConfigSlurper().parse(configStr)

        when:
        def script = '''

            process foo {
                exec:
                println true
            }

            '''

        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor
        println processor.config.ext
        then:
        processor instanceof TaskProcessor
        processor.config.ext.alpha == 'aaa'
        processor.config.ext.beta == 'BBB'
        processor.config.ext.delta == 'DDD'
    }


    def 'test configure processor with dynamic resources'() {

        setup:
        def configStr = '''
             process {
                cpus = { 2 * task.attempt }
                memory = { 1.GB * task.attempt  }
                time = { 1.h * task.attempt }

                $taskHello.errorStrategy = 'finish'
            }
            '''
        def cfg = new ConfigSlurper().parse(configStr)


        when:
        def script = '''

            process taskHello {
                cpus 2
                memory 3.GB
                time '1 h'
                errorStrategy 'ignore'

                script:
                'echo hello'
            }

            '''

        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.config.cpus == 2
        processor.config.memory == MemoryUnit.of('3 GB')
        processor.config.time == '1 h'
        processor.config.errorStrategy == 'finish'

    }

}
