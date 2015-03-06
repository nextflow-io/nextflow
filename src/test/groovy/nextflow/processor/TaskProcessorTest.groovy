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

package nextflow.processor

import static test.TestParser.parse

import java.nio.file.Files
import java.nio.file.Paths

import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.executor.NopeExecutor
import nextflow.file.FileHolder
import nextflow.script.BaseScript
import nextflow.script.FileInParam
import nextflow.script.ScriptRunner
import nextflow.script.TaskBody
import nextflow.script.TokenValRef
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.util.CacheHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {

    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(String name, Session session, BaseScript script, ProcessConfig taskConfig) {
            super(name, new NopeExecutor(), session, new DummyScript(), taskConfig, new TaskBody({}, '..', true))
        }

        @Override protected void createOperator() { }
    }

    static class DummyScript extends BaseScript {
        @Override Object run() { return null }
    }


    def filterHidden() {

        setup:
        def processor = [:] as TaskProcessor
        def list = [ Paths.get('file.txt'), Paths.get('.hidden'), Paths.get('file.name') ]

        when:
        def result = processor.filterByRemovingHiddenFiles(list)
        then:
        result == [ Paths.get('file.txt'), Paths.get('file.name')  ]

    }

    def filterStagedInputs() {

        setup:
        def processor = [:] as TaskProcessor
        def binding = new Binding()
        def holder = []

        def inputs = [:]
        def key1 = new FileInParam(binding, holder).bind('file1')
        def key2 = new FileInParam(binding, holder).bind('file_')
        def key3 = new ValueInParam(binding, holder).bind( new TokenVar('xxx') )

        def val1 = [ FileHolder.get('xxx', 'file.txt') ]
        def val2 =  [ FileHolder.get('yyy', 'file.2'), FileHolder.get('zzz', '.hidden') ]
        def val3 =  'just a value'
        inputs[key1] = val1
        inputs[key2] = val2
        inputs[key3] = val3

        def task = [:] as TaskRun
        task.inputs = inputs

        when:
        // three files have been produced
        def files = [ Paths.get('file.1'), Paths.get('file.2'), Paths.get('file.3') ]
        def result = processor.filterByRemovingStagedInputs(task, files)
        then:
        // the *file.2* is removed since it belongs to the inputs list
        result == [ Paths.get('file.1'), Paths.get('file.3')  ]

    }


    def testEnvironment() {

        setup:
        def home = Files.createTempDirectory('test')
        def binFolder = home.resolve('bin')
        binFolder.mkdirs()

        when:
        def wrapper = new DummyScript()
        def session = new Session([env: [X:"1", Y:"2"]])
        session.setBaseDir(home)
        def processor = new DummyProcessor('task1', session, wrapper, new ProcessConfig(wrapper))
        def builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )

        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == "${binFolder.toString()}:\$PATH"

        when:
        session = new Session([env: [X:"1", Y:"2", PATH:'/some']])
        session.setBaseDir(home)
        processor = new DummyProcessor('task1', session, wrapper, new ProcessConfig(wrapper))
        builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )

        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == "${binFolder.toString()}:/some"


        cleanup:
        home.deleteDir()

    }

    def testFetchInterpreter() {

        when:
        def script =
            '''
            #!/bin/perl
            do this
            do that
            '''
        def i = TaskProcessor.fetchInterpreter(script.stripIndent().trim())
        then:
        i == '/bin/perl'

        when:
        i = TaskProcessor.fetchInterpreter('do this')
        then:
        i == null
    }

    def testSingleItemOrCollection() {

        setup:
        def processor = [:] as TaskProcessor
        def path1 = Paths.get('file1')
        def path2 = Paths.get('file2')
        def path3 = Paths.get('file3')

        when:
        def list = [ FileHolder.get(path1, 'x_file_1') ]
        def result = processor.singleItemOrList(list)
        then:
        result.toString() == 'x_file_1'

        when:
        list = [ FileHolder.get(path1, 'x_file_1'), FileHolder.get(path2, 'x_file_2'), FileHolder.get(path3, 'x_file_3') ]
        result = processor.singleItemOrList(list)
        then:
        result*.toString() == [ 'x_file_1',  'x_file_2',  'x_file_3']

    }


    def testWildcardExpand() {

        setup:
        def processor = [:] as TaskProcessor

        /*
         * The name do not contain any wildcards *BUT* when multiple files are provide
         * an index number is added to the specified name
         */
        when:
        def list1 = processor.expandWildcards('file_name', [FileHolder.get('x')])
        def list2 = processor.expandWildcards('file_name', [FileHolder.get('x'), FileHolder.get('y')] )
        then:
        list1 *. stageName  == ['file_name']
        list2 *. stageName  == ['file_name1', 'file_name2']


        /*
         * The star wildcard: when a single item is provided, it is simply ignored
         * When a collection of files is provided, the name is expanded to the index number
         */
        when:
        list1 = processor.expandWildcards('file*.fa', [FileHolder.get('x')])
        list2 = processor.expandWildcards('file_*.fa', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 *. stageName == ['file.fa']
        list2 *. stageName == ['file_1.fa', 'file_2.fa', 'file_3.fa']

        /*
         * The question mark wildcards *always* expand to an index number
         */
        when:
        def p0 = [FileHolder.get('0')]
        def p1_p4 = (1..4).collect { FileHolder.get(it.toString()) }
        def p1_p12 = (1..12).collect { FileHolder.get(it.toString()) }
        list1 = processor.expandWildcards('file?.fa', p0 )
        list2 = processor.expandWildcards('file_???.fa', p1_p4 )
        def list3 = processor.expandWildcards('file_?.fa', p1_p12 )
        then:
        list1 *. stageName == ['file1.fa']
        list2 *. stageName == ['file_001.fa', 'file_002.fa', 'file_003.fa', 'file_004.fa']
        list3 *. stageName == ['file_1.fa', 'file_2.fa', 'file_3.fa', 'file_4.fa', 'file_5.fa', 'file_6.fa', 'file_7.fa', 'file_8.fa', 'file_9.fa', 'file_10.fa', 'file_11.fa', 'file_12.fa']

        when:
        list1 = processor.expandWildcards('*', [FileHolder.get('a')])
        list2 = processor.expandWildcards('*', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 *. stageName == ['a']
        list2 *. stageName == ['x','y','z']


    }


    def testStagePath() {

        when:
        def p1 = Paths.get('/home/data/source.file')
        def path = FileHolder.get(p1, 'target.file')

        then:
        path.storePath == p1
        path.stageName == 'target.file'

    }


    def testGetHashLog() {

        when:
        def h = CacheHelper.hasher('x').hash()
        def task = new TaskRun(hash:h)
        then:
        task.getHashLog() == '76/9f897d'

    }


    def testState() {

        when:
        def state = new Agent<StateObj>(new StateObj())
        int i = 0
        state.addListener { a, b -> i++ }

        state.update { StateObj it ->  it.incSubmitted()  }
        state.update { StateObj it ->  it.incCompleted() }
        state.update  { StateObj it ->  it.poison()  }
        state.await()
        then:
        state.val.finished
        i == 3

        when:
        state = new Agent<StateObj>(new StateObj())
        state.update { StateObj it ->  it.incSubmitted()  }
        state.update  { StateObj it ->  it.poison()  }
        state.await()
        then:
        !state.val.finished

    }

    def testState2() {

        when:
        def agent = new Agent<List>([])
        int i = 0
        agent.addListener { a, b -> println ">>: $a -- $b"; i++ }

        agent << { it.add(1); (this as Agent).updateValue(it) }
        agent << { it.add(2); (this as Agent).updateValue(it) }
        agent.await()
        then:
        agent.val == [1,2]

    }


    def testVarRefs( ) {

        setup:
        def text = '''

        String x = 1

        @Field
        String = 'Ciao'

        z = 'str'

        process hola {

          /
          println $x + $y + $z
          /
        }
        '''
        when:
        TaskProcessor process = parse(text).run()
        then:
        process.taskBody.valRefs == [
                new TokenValRef('x', 13, 20),
                new TokenValRef('y', 13, 25),
                new TokenValRef('z', 13, 30) ] as Set

        process.taskBody.getValNames() == ['x','y','z'] as Set
    }

    def testProcessVariables() {

        when:
        def runner = new ScriptRunner( process: [executor:'nope'] )
        def script =
                '''
                class Foo { def foo() { return [x:1] };  }

                alpha = 1
                params.beta = 2
                params.zeta = new Foo()
                delta = new Foo()
                x = 'alpha'

                process simpleTask  {
                    input:
                    val x from 1

                    """
                    echo ${alpha}
                    echo ${params.beta}
                    echo ${params?.gamma?.omega}
                    echo ${params.zeta.foo()}
                    echo ${params.zeta.foo().x}
                    echo ${delta.foo().x}
                    echo ${params."$x"}
                    """
                }

                '''
        runner.setScript(script).execute()
        then:
        runner.getScriptObj().getTaskProcessor().getTaskBody().getValNames() == ['alpha', 'params.beta', 'params.gamma.omega', 'params.zeta', 'delta', 'params', 'x'] as Set

    }

    def testTaskGlobalVars() {

        setup:
        def processor = [:] as TaskProcessor
        processor.name = 'Hello'

        when:
        def binding = new Binding(x:1, y:2, params: [alpha: 'one'], 'workDir': Paths.get('/work/dir'), baseDir: Paths.get('/base/dir'))
        def vars = ['q','x','y','params.alpha','params.beta.delta', 'workDir', 'baseDir'] as Set
        def local = [q: 'any value']
        def result = processor.getTaskGlobalVars(binding, vars, local)

        then:
        // note: since 'q' is include in the task local scope, is not returned in the var list
        result == [x:1, y:2, 'params.alpha': 'one', 'params.beta.delta': null , baseDir: '/base/dir', workDir: '/work/dir']

    }



}
