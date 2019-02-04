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


import spock.lang.Shared
import spock.lang.Specification

import nextflow.config.ConfigParser
import nextflow.processor.TaskProcessor
import nextflow.script.ScriptRunner
import nextflow.util.MemoryUnit

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

    def 'should define default ext property' () {

        given:
        def CONFIG = '''
            process.ext.foo = 'hello'
        '''
        def cfg = new ConfigSlurper().parse(CONFIG)

        when:
        def script = '''

            process foo {
                script:
                /
                echo ciao
                /
            }

            '''

        def runner = new ScriptRunner(cfg).setScript(script)
        runner.execute()
        def processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.ext.foo == 'hello'

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

    def 'should set config with labels' () {

        given:
        /*
         * A config with a `memory` definition for all process
         * and two labels `small` and `big`
         */
        String CONFIG = '''
            process {
                executor = 'nope'
                memory = 2.GB
                
                withLabel: small {
                    cpus = 2 
                    queue = 'the-small-one'
                }
                
                withLabel: big {
                    cpus = 8 
                    memory = 4.GB
                    queue = 'big-partition'                
                }
                
                $legacy.cpus = 3 
                $legacy.queue = 'legacy-queue'
            }
            '''


        when:
            /*
             * no label is specified it should only use default directives
             */
            String script = '''   

                process foo {
                    script:
                    'echo hello'
                }
                '''

            def cfg = new ConfigParser().parse(CONFIG)
            def runner = new ScriptRunner(cfg)
            runner.setScript(script).execute()
            def processor = runner.scriptObj.taskProcessor

        then:
            processor instanceof TaskProcessor
            processor.config.memory == MemoryUnit.of('2 GB')
            processor.config.cpus == null
            processor.config.queue == null

        when:
            /*
             * the `small` label is applied
             */
            script = '''   

                process foo {
                    label 'small'
                    script:
                    'echo hello'
                }
                '''

            cfg = new ConfigParser().parse(CONFIG)
            runner = new ScriptRunner(cfg)
            runner.setScript(script).execute()
            processor = runner.scriptObj.taskProcessor

        then:
            processor instanceof TaskProcessor
            processor.config.memory == MemoryUnit.of('2 GB')
            processor.config.cpus == 2
            processor.config.queue == 'the-small-one'


        when:
            /*
             * the directives for the `big` label are applied
             */
            script = '''
                process foo {
                    label 'big'
                    script:
                    'echo hello'
                }
                '''

            cfg = new ConfigParser().parse(CONFIG)
            runner = new ScriptRunner(cfg)
            runner.setScript(script).execute()
            processor = runner.scriptObj.taskProcessor

        then:
            processor instanceof TaskProcessor
            processor.config.cpus == 8
            processor.config.memory == MemoryUnit.of('4 GB')
            processor.config.queue == 'big-partition'


        when:
        /*
        * the directives for the `big` label are applied
        * moreover the directive for `bar` overrides the previous ones
        */
        script = '''
                process legacy {
                    cpus 1 
                    queue 'one'
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.config.memory == MemoryUnit.of('2 GB')
        processor.config.cpus == 3
        processor.config.queue == 'legacy-queue'
    }

    def 'should set setting for process with name' () {


        given:
        /*
         * A config with a `memory` definition for all process
         * and two labels `small` and `big`
         */
        String CONFIG = '''
            process {
                executor = 'nope' 
                
                withLabel: small {
                    cpus = 2 
                    queue = 'the-small-one'
                }
                
                withName: bar {
                    cpus = 8 
                    memory = 4.GB
                    queue = 'big-partition'                
                }
            }
            '''


        when:
        /*
         * no label is specified it should only use default directives
         */
        String script = '''   

                process foo {
                    label 'small'
                    script:
                    'echo hello'
                }
                '''

        def cfg = new ConfigParser().parse(CONFIG)
        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.config.cpus == 2
        processor.config.queue == 'the-small-one'


        when:
        /*
         * no label is specified it should only use default directives
         */
        script = '''   

                process bar {
                    label 'small'
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.config.cpus == 8
        processor.config.queue == 'big-partition'
    }

    def 'should set module directive' () {
        when:
        String config = '''
            process {
                executor = 'nope'
                withLabel: 'my_env' {
                    module = 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'
                }
            }
            '''
        String script = '''   
                process foo {
                    module 'mod-a/1.1:mod-b/2.2'
                    script:
                    'echo hello'
                }
                '''

        def cfg = new ConfigParser().parse(config)
        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.config.module == ['mod-a/1.1:mod-b/2.2']
        processor.config.createTaskConfig().getModule() == ['mod-a/1.1','mod-b/2.2']


        when:
        config = '''
            process {
                executor = 'nope'
                withLabel: 'my_env' {
                    module = 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'
                }
            }
            '''
        script = '''   
                process foo {
                    label 'my_env'
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(config)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor

        then:
        processor instanceof TaskProcessor
        processor.config.module == ['ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1']
        processor.config.createTaskConfig().getModule() == ['ncbi-blast/2.2.27','t_coffee/10.0','clustalw/2.1']

    }

    def 'should set publishDir directive' () {
        when:
        String CONFIG = '''
            process {
                executor = 'nope'
                publishDir = '/some/dir'
            }
            '''
        String script = '''   
                process foo {
                    script:
                    'echo hello'
                }
                '''

        def cfg = new ConfigParser().parse(CONFIG)
        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/some/dir']


        when:
        CONFIG = '''
            process {
                executor = 'nope'
                publishDir = [ '/some/dir', [path:'/other/dir', mode: 'copy'] ]
            }
            '''
        script = '''   
                process foo {
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/some/dir']
        processor.config.publishDir[1] == [path:'/other/dir', mode: 'copy']


        when:
        CONFIG = '''
            process {
                executor = 'nope'
            }
            '''
        script = '''   
                process foo {
                    publishDir '/data1'
                    publishDir '/data2', mode: 'symlink'
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/data1']
        processor.config.publishDir[1] == [path:'/data2', mode: 'symlink']


        when:
        CONFIG = '''
            process {
                executor = 'nope'
                publishDir = '/dir/cfg'
            }
            '''
        script = '''   
                process foo {
                    publishDir '/dir/alpha'
                    publishDir '/dir/bravo'  
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/dir/alpha']
        processor.config.publishDir[1] == [path:'/dir/bravo']
        processor.config.publishDir.size() == 2



        when:
        CONFIG = '''
            process {
                executor = 'nope'
                publishDir = '/dir/cfg'
                withName: foo { publishDir = '/dir/omega' }
            }
            '''
        script = '''   
                process foo {
                    publishDir '/dir/alpha'
                    publishDir '/dir/bravo'  
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/dir/omega']
        processor.config.publishDir.size() == 1
    }



    def 'should set directive label' () {

        when:
        def CONFIG = '''
            process {
                executor = 'nope'
                label = 'alpha'
            }
            '''
        def script = '''   
                process foo {
                    script:
                    'echo hello'
                }
                '''

        def cfg = new ConfigParser().parse(CONFIG)
        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        def processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.label == [ 'alpha' ]


        when:
        CONFIG = '''
            process {
                executor = 'nope'
                label = 'alpha'
            }
            '''
        script = '''   
                process foo {
                    label 'bravo'
                    label 'gamma'  
                    script:
                    'echo hello'
                }
                '''

        cfg = new ConfigParser().parse(CONFIG)
        runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        processor = runner.scriptObj.taskProcessor
        then:
        processor instanceof TaskProcessor
        processor.config.label.size() == 2
        processor.config.label == [ 'bravo', 'gamma' ]
    }

    def 'should create process with repeater'() {

        when:
        def CONFIG = '''
            process {
                executor = 'nope'
            }
            '''

        def script = '''   
                process foo {
                    input:
                    each x from (1,2,3)
                    script:
                    'echo hello'
                }
                '''

        def cfg = new ConfigParser().parse(CONFIG)
        def runner = new ScriptRunner(cfg)
        runner.setScript(script).execute()
        then:
        noExceptionThrown()
    }
}
