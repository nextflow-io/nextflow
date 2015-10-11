/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths
import java.util.concurrent.ExecutorService

import groovyx.gpars.agent.Agent
import nextflow.Global
import nextflow.ISession
import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.executor.NopeExecutor
import nextflow.file.FileHolder
import nextflow.script.BaseScript
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.TaskBody
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.util.CacheHelper
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {

    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(String name, Session session, BaseScript script, ProcessConfig taskConfig) {
            super(name, new NopeExecutor(), session, script, taskConfig, new TaskBody({}, '..'))
        }

        @Override protected void createOperator() { }
    }


    def 'should filter hidden files'() {

        setup:
        def processor = [:] as TaskProcessor
        def list = [ Paths.get('file.txt'), Paths.get('.hidden'), Paths.get('file.name') ]

        when:
        def result = processor.filterByRemovingHiddenFiles(list)
        then:
        result == [ Paths.get('file.txt'), Paths.get('file.name')  ]

    }

    def 'should filter staged inputs'() {

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


    def 'should evaluate environment variables'() {

        setup:
        def home = Files.createTempDirectory('test')
        def binFolder = home.resolve('bin')
        binFolder.mkdirs()

        when:
        def session = new Session([env: [X:"1", Y:"2"]])
        session.setBaseDir(home)
        def processor = new DummyProcessor('task1', session, Mock(BaseScript), Mock(ProcessConfig))
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
        processor = new DummyProcessor('task1', session,  Mock(BaseScript), Mock(ProcessConfig))
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

    def 'should fetch interpreter from shebang line'() {

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

    def 'should return single item or collection'() {

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


    def 'should expand wildcards'() {

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


    def 'shold stage path'() {

        when:
        def p1 = Paths.get('/home/data/source.file')
        def path = FileHolder.get(p1, 'target.file')

        then:
        path.storePath == p1
        path.stageName == 'target.file'

    }


    def 'should return task hash log'() {

        when:
        def h = CacheHelper.hasher('x').hash()
        def task = new TaskRun(hash:h)
        then:
        task.getHashLog() == '76/9f897d'

    }


    def 'should update agent state'() {

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

    def 'should update agent state 2'() {

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


    def 'should return tasks global variables map'() {

        given:
        def processor = [:] as TaskProcessor
        processor.name = 'Hello'

        def binding = new Binding(x:1, y:2, params: [alpha: 'one'], 'workDir': Paths.get('/work/dir'), baseDir: Paths.get('/base/dir'))
        def vars = ['q', 'x', 'y', 'params.alpha', 'params.beta.delta', 'workDir', 'baseDir'] as Set

        when:
        def result = processor.getTaskGlobalVars(vars, binding, [:])
        then:
        // note: since 'q' is include in the task local scope, is not returned in the var list
        result == [x:1, y:2, 'params.alpha': 'one', 'params.beta.delta': null , baseDir: '/base/dir', workDir: '/work/dir']

        when:
        result = processor.getTaskGlobalVars(vars, binding, [x:'foo',y:'bar'])
        then:
        result == ['params.alpha': 'one', 'params.beta.delta': null , baseDir: '/base/dir', workDir: '/work/dir']

    }


    def "should return a file holder" () {

        given:
        FileHolder holder
        def tempFolder = Files.createTempDirectory('test')
        def localFile = Files.createTempFile(tempFolder, 'test','test')
        Global.session = Mock(ISession)
        Global.session.workDir >> tempFolder
        def processor = [:] as TaskProcessor

        /*
         * when the input file is on the local file system
         * simple return a reference to it in the holder object
         */
        when:
        holder = processor.normalizeInputToFile(localFile,null)
        then:
        holder.sourceObj == localFile
        holder.storePath == localFile
        holder.stageName == localFile.getFileName().toString()

        /*
         * when the input is a on a foreign file system
         * it is need to copy it to a local file
         */
        when:
        def remoteFile = TestHelper.createInMemTempFile('remote_file.txt')
        remoteFile.text = 'alpha beta gamma delta'
        holder = processor.normalizeInputToFile(remoteFile,'input.1')
        then:
        holder.sourceObj == remoteFile
        holder.storePath.fileSystem == FileSystems.default
        holder.storePath.text == 'alpha beta gamma delta'
        holder.stageName == remoteFile.getName()

        /*
         * any generic input that is not a file is converted to a string
         * and save to the local file system
         */
        when:
        holder = processor.normalizeInputToFile("text data string",'simple_file_name.txt')
        then:
        holder.sourceObj == "text data string"
        holder.storePath.fileSystem == FileSystems.default
        holder.storePath.text == "text data string"
        holder.stageName == 'simple_file_name.txt'

        cleanup:
        tempFolder?.deleteDir()
    }


    def 'should return `ignore` strategy' () {

        given:
        def task
        def proc = [:] as TaskProcessor
        def error = Mock(ProcessException)

        when:
        task = new TaskRun()
        task.config = new TaskConfig()
        then:
        proc.checkErrorStrategy(task, error, 1,1) == null

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy: 'ignore')
        then:
        proc.checkErrorStrategy(task, error, 10, 10) == ErrorStrategy.IGNORE

    }

    def 'should return `retry` strategy' () {

        given:

        def task
        def error = Mock(ProcessException)
        def session = Mock(Session)
        session.getExecService() >> Mock(ExecutorService)

        def proc = [:] as TaskProcessor
        proc.session = session

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy:'retry', maxErrors: max_errors, maxRetries: max_retries )
        then:
        proc.checkErrorStrategy(task, error, task_err_count , proc_err_count) == strategy

        where:
        max_retries | max_errors    |   task_err_count  |  proc_err_count   | strategy
                1   |        3      |               0   |               0   | ErrorStrategy.RETRY
                1   |        3      |               1   |               0   | null
                1   |        3      |               0   |               1   | ErrorStrategy.RETRY
                1   |        3      |               0   |               2   | ErrorStrategy.RETRY
                1   |        3      |               0   |               3   | null
                3   |       -1      |               0   |               0   | ErrorStrategy.RETRY
                3   |       -1      |               1   |               1   | ErrorStrategy.RETRY
                3   |       -1      |               2   |               2   | ErrorStrategy.RETRY
                3   |       -1      |               2   |               9   | ErrorStrategy.RETRY
                3   |       -1      |               3   |               9   | null

    }



    def 'should return the list of output files'() {

        given:
        def param
        def result

        def folder = Files.createTempDirectory('test')
        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        folder.resolve('.hidden.fa').text = 'hidden'
        folder.resolve('dir1').mkdir()
        folder.resolve('dir1').resolve('file3.txt').text = 'file 3'
        folder.resolve('dir1')
        folder.resolve('dir1').resolve('dir2').mkdirs()
        folder.resolve('dir1').resolve('dir2').resolve('file4.fa').text = 'file '
        Files.createSymbolicLink( folder.resolve('dir_link'), folder.resolve('dir1') )

        def processor = [:] as TaskProcessor

        when:
        result = processor.fetchResultFiles(Mock(FileOutParam), '*.fa', folder )
        then:
        result.collect { it.name }  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('file')
        result = processor.fetchResultFiles(param, '*.fa', folder)
        then:
        result.collect { it.name }  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('dir')
        result = processor.fetchResultFiles(param, '*.fa', folder)
        then:
        result == []

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = processor.fetchResultFiles(param, '**.fa', folder)
        then:
        result.collect { it.name }.sort()  == ['file2.fa','file4.fa','file4.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.followLinks(false)
        result = processor.fetchResultFiles(param, '**.fa', folder)
        then:
        result.collect { it.name }.sort()  == ['file2.fa','file4.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.maxDepth(1)
        result = processor.fetchResultFiles(param, '**.fa', folder)
        then:
        result.collect { it.name }.sort()  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = processor.fetchResultFiles(param, '*', folder)
        then:
        result.collect { it.name }.sort()  == ['dir1', 'dir_link', 'file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('dir')
        result = processor.fetchResultFiles(param, '*', folder)
        then:
        result.collect { it.name }.sort()  == ['dir1', 'dir_link']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('file')
        result = processor.fetchResultFiles(param, '*', folder)
        then:
        result.collect { it.name }.sort()  == ['file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('file')
        param.hidden(true)
        result = processor.fetchResultFiles(param, '*', folder)
        then:
        result.collect { it.name }.sort()  == ['.hidden.fa', 'file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = processor.fetchResultFiles(param,'.*', folder)
        then:
        result.collect { it.name }.sort()  == ['.hidden.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = processor.fetchResultFiles(param,'file{1,2}.{txt,fa}', folder)
        then:
        result.collect { it.name }.sort() == ['file1.txt', 'file2.fa']

        cleanup:
        folder?.deleteDir()

    }

    def 'should create the map of path visit options'() {

        given:
        def param
        def processor = [:] as TaskProcessor

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        then:
        processor.visitOptions(param,'file.txt') == [type:'any', followLinks: true, maxDepth: null, hidden: false, relative: false]
        processor.visitOptions(param,'path/**') == [type:'file', followLinks: true, maxDepth: null, hidden: false, relative: false]
        processor.visitOptions(param,'.hidden_file') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('dir')
        then:
        processor.visitOptions(param,'dir-name') == [type:'dir', followLinks: true, maxDepth: null, hidden: false, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.hidden(true)
        then:
        processor.visitOptions(param,'dir-name') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.followLinks(false)
        then:
        processor.visitOptions(param,'dir-name') == [type:'any', followLinks: false, maxDepth: null, hidden: false, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.maxDepth(5)
        then:
        processor.visitOptions(param,'dir-name') == [type:'any', followLinks: true, maxDepth: 5, hidden: false, relative: false]
    }



}
