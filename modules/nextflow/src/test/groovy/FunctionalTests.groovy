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

import nextflow.config.ConfigParser
import nextflow.exception.AbortRunException
import nextflow.processor.TaskProcessor
import nextflow.util.MemoryUnit
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class FunctionalTests extends Dsl2Spec {

    /*
     * test passing values through command line argument (unnamed parameters)
     */
    def 'test args'()  {

        given:
        def script = """
            def len = args.size()
            def x = args[0]
            def y = args[1]

            return [ x, y, len ]
            """

        when:
        def result = new MockScriptRunner().setScript(script).execute(['hello', 'hola'])

        then:
        result[0] == 'hello'
        result[1] == 'hola'
        result[2] == 2

    }


    def 'test configure processor'() {

        given:
        def config = '''
             process {
                dummyField = 99
                executor = 'nope'
                debug = true
                shell = 'zsh'
                maxForks = 10
                environment = [a:1, b:2,c:3]
            }
            '''

        and:
        def script = '''

            process taskHello {
                debug true
                maxForks 11
                'echo hello'
            }

            workflow { taskHello() }
            '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()

        then:
        processor instanceof TaskProcessor
        processor.getName() == 'taskHello'
        processor.config.debug == true
        processor.config.shell == 'zsh'
        processor.config.maxForks == 11
        processor.config.dummyField == 99
        processor.config.environment.entrySet() == [a:1,b:2,c:3].entrySet()

    }

    def 'should define default ext property' () {

        given:
        def config = '''
            process.ext.foo = 'hello'
        '''
        and:
        def script = '''

            process foo {
                script:
                /
                echo ciao
                /
            }

            workflow { foo() }
            '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.ext.foo == 'hello'

    }


    def 'test merge ext properties' () {
        given:
        def config = '''
            process {
                ext {
                    alpha = "aaa"
                    delta = "ddd"
                }
                withName: foo {
                    ext {
                        beta = "BBB"
                        delta = "DDD"
                    }
                }
            }
        '''
        and:
        def script = '''

            process foo {
                exec:
                println true
            }

            workflow { foo() }
            '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.ext.alpha == 'aaa'
        processor.config.ext.beta == 'BBB'
        processor.config.ext.delta == 'DDD'
    }


    def 'test configure processor with dynamic resources'() {

        setup:
        def config = '''
             process {
                cpus = { 2 * task.attempt }
                memory = { 1.GB * task.attempt  }
                time = { 1.h * task.attempt }
                withName: taskHello{ errorStrategy = 'finish' }
            }
            '''

        and:
        def script = '''

            process taskHello {
                cpus 2
                memory 3.GB
                time '1 h'
                errorStrategy 'ignore'

                script:
                'echo hello'
            }

            workflow { taskHello() }
            '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()

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
        String config = '''
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
                
                withName: legacy {
                    cpus = 3 
                    queue = 'legacy-queue'
                }
            }
            '''

        and:
            /*
             * no label is specified it should only use default directives
             */
            String script = '''   

                process foo {
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        when:
            new MockScriptRunner(new ConfigParser().parse(config))
                    .setScript(script)
                    .execute()
            def processor = TaskProcessor.currentProcessor()

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
                
                workflow { foo() }
                '''

        and:
            new MockScriptRunner(new ConfigParser().parse(config))
                    .setScript(script)
                    .execute()
            processor = TaskProcessor.currentProcessor()

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
                
                workflow { foo() }
                '''
        and:
            new MockScriptRunner(new ConfigParser().parse(config))
                    .setScript(script)
                    .execute()
            processor = TaskProcessor.currentProcessor()

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
                
                workflow { legacy() }
                '''

        and:
            new MockScriptRunner(new ConfigParser().parse(config))
                    .setScript(script)
                    .execute()
            processor = TaskProcessor.currentProcessor()

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
        String config = '''
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

        and:
        /*
         * no label is specified it should only use default directives
         */
        String script = '''   

                process foo {
                    label 'small'
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()

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
                
                workflow { bar() }
                '''
        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()

        then:
        processor instanceof TaskProcessor
        processor.config.cpus == 8
        processor.config.queue == 'big-partition'
    }

    def 'should set module directive' () {
        given:
        String config = '''
            process {
                executor = 'nope'
                withLabel: 'my_env' {
                    module = 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'
                }
            }
            '''
        and:
        String script = '''   
                process foo {
                    module 'mod-a/1.1:mod-b/2.2'
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()
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
        and:
        script = '''   
                process foo {
                    label 'my_env'
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()

        then:
        processor instanceof TaskProcessor
        processor.config.module == ['ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1']
        processor.config.createTaskConfig().getModule() == ['ncbi-blast/2.2.27','t_coffee/10.0','clustalw/2.1']

    }

    def 'should set publishDir directive' () {
        given:
        String config = '''
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
                
                workflow { foo() }
                '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/some/dir']


        when:
        config = '''
            process {
                executor = 'nope'
                publishDir = [ '/some/dir', [path:'/other/dir', mode: 'copy'] ]
            }
            '''
        and:
        script = '''   
                process foo {
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''
        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/some/dir']
        processor.config.publishDir[1] == [path:'/other/dir', mode: 'copy']


        when:
        config = '''
            process {
                executor = 'nope'
            }
            '''
        and:
        script = '''   
                process foo {
                    publishDir '/data1'
                    publishDir '/data2', mode: 'symlink'
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/data1']
        processor.config.publishDir[1] == [path:'/data2', mode: 'symlink']


        when:
        config = '''
            process {
                executor = 'nope'
                publishDir = '/dir/cfg'
            }
            '''
        and:
        script = '''   
                process foo {
                    publishDir '/dir/alpha'
                    publishDir '/dir/bravo'  
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/dir/alpha']
        processor.config.publishDir[1] == [path:'/dir/bravo']
        processor.config.publishDir.size() == 2


        when:
        config = '''
            process {
                executor = 'nope'
                publishDir = '/dir/cfg'
                withName: foo { publishDir = '/dir/omega' }
            }
            '''
        and:
        script = '''   
                process foo {
                    publishDir '/dir/alpha'
                    publishDir '/dir/bravo'  
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.publishDir[0] == [path:'/dir/omega']
        processor.config.publishDir.size() == 1
    }


    def 'should set directive label' () {

        given:
        def config = '''
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
                
                workflow { foo() }
                '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        def processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.label == [ 'alpha' ]


        when:
        config = '''
            process {
                executor = 'nope'
                label = 'alpha'
            }
            '''
        and:
        script = '''   
                process foo {
                    label 'bravo'
                    label 'gamma'  
                    script:
                    'echo hello'
                }
                
                workflow { foo() }
                '''

        and:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        processor = TaskProcessor.currentProcessor()
        then:
        processor instanceof TaskProcessor
        processor.config.label.size() == 2
        processor.config.label == [ 'bravo', 'gamma' ]
    }

    def 'should create process with repeater'() {

        given:
        def config = '''
            process {
                executor = 'nope'
            }
            '''
        and:
        def script = '''   
                process foo {
                    input:
                    each x
                    script:
                    'echo hello'
                }
                
                workflow { foo([1,2,3]) }
                '''

        when:
        new MockScriptRunner(new ConfigParser().parse(config))
                .setScript(script)
                .execute()
        then:
        noExceptionThrown()
    }

    def 'should show the line of the error when throw an exception'() {

        given:
        def script = '''/*1*/
/*2*/   def thisMethodExpectsOnlyOneString(String a){
/*3*/      a
/*4*/   }
/*5*/                   
/*6*/   process foo {
/*7*/       input:
/*8*/           each x
/*9*/       script:
/*10*/          "${thisMethodExpectsOnlyOneString(1)}"
/*11*/      }
/*12*/
/*13*/   workflow { foo(1) }
        '''

        when:
        def config = [process:[executor: 'nope']]
        def runner = new MockScriptRunner(config)
        runner.setScript(script).execute()
        then:
        def abort = thrown(AbortRunException)
        and:
        runner.session.fault.report ==~ /(?s).*-- Check script '(.*?)' at line: 10.*/
    }
}
