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

package nextflow.processor

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ExecutorService

import groovyx.gpars.agent.Agent
import nextflow.Global
import nextflow.ISession
import nextflow.Session
import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.exception.ProcessEvalException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.executor.NopeExecutor
import nextflow.file.FileHolder
import nextflow.file.FilePorter
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ScriptType
import nextflow.script.bundle.ResourcesBundle
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.util.ArrayBag
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {

    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(String name, Session session, BaseScript script, ProcessConfig taskConfig) {
            super(name, new NopeExecutor(session: session), session, script, taskConfig, new BodyDef({}, '..'))
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

        given:
        def task = Spy(TaskRun)
        def processor = [:] as TaskProcessor

        def WORK_DIR = Paths.get('/work/dir')
        def FILE1 = WORK_DIR.resolve('alpha.txt')
        def FILE2 = WORK_DIR.resolve('beta.txt')
        def FILE3 = WORK_DIR.resolve('out/beta.txt')
        def FILE4 = WORK_DIR.resolve('gamma.fasta')

        def collectedFiles = [ FILE1, FILE2, FILE3, FILE4 ]

        when:
        def result = processor.filterByRemovingStagedInputs(task, collectedFiles, WORK_DIR)
        then:
        1 * task.getStagedInputs() >> [ 'beta.txt' ]
        result == [ FILE1, FILE3, FILE4 ]

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
        builder.environment().PATH == "\$PATH:${binFolder.toString()}"

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
        builder.environment().PATH == "/some:${binFolder.toString()}"

        cleanup:
        home.deleteDir()

    }

    @Unroll
    def 'should add module bin paths to task env' () {
        given:
        def session = Mock(Session) { getConfig() >> [:] }
        def executor = Mock(Executor) { getBinDir() >> Path.of('/project/bin')}
        and:
        TaskProcessor processor = Spy(TaskProcessor, constructorArgs: [[session:session, executor:executor]])
        and:
        when:
        def result = processor.getProcessEnvironment()
        then:
        session.enableModuleBinaries() >> MODULE_BIN
        processor.getModuleBundle() >> Mock(ResourcesBundle)  { getBinDirs() >> [Path.of('/foo'), Path.of('/bar')] }
        processor.isLocalWorkDir() >> LOCAL
        and:
        result == EXPECTED

        where:
        LOCAL   | MODULE_BIN    | EXPECTED
        false   | false         | [:]
        true    | false         | [PATH:'$PATH:/project/bin']
        true    | true          | [PATH:'$PATH:/foo:/bar:/project/bin']
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
        def result = processor.singleItemOrList(list, true, ScriptType.SCRIPTLET)
        then:
        result.toString() == 'x_file_1'

        when:
        list = [ FileHolder.get(path1, 'x_file_1') ]
        result = processor.singleItemOrList(list, false, ScriptType.SCRIPTLET)
        then:
        result*.toString() == ['x_file_1']

        when:
        list = [ FileHolder.get(path1, 'x_file_1'), FileHolder.get(path2, 'x_file_2'), FileHolder.get(path3, 'x_file_3') ]
        result = processor.singleItemOrList(list, false, ScriptType.SCRIPTLET)
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
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
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
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list3 instanceof ArrayBag
        list1 *. stageName == ['file1.fa']
        list2 *. stageName == ['file_001.fa', 'file_002.fa', 'file_003.fa', 'file_004.fa']
        list3 *. stageName == ['file_1.fa', 'file_2.fa', 'file_3.fa', 'file_4.fa', 'file_5.fa', 'file_6.fa', 'file_7.fa', 'file_8.fa', 'file_9.fa', 'file_10.fa', 'file_11.fa', 'file_12.fa']

        when:
        list1 = processor.expandWildcards('*', [FileHolder.get('a')])
        list2 = processor.expandWildcards('*', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['a']
        list2 *. stageName == ['x','y','z']

        when:
        list1 = processor.expandWildcards('dir1/*', [FileHolder.get('a')])
        list2 = processor.expandWildcards('dir2/*', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['dir1/a']
        list2 *. stageName == ['dir2/x','dir2/y','dir2/z']

        when:
        list1 = processor.expandWildcards('/dir/file*.fa', [FileHolder.get('x')])
        list2 = processor.expandWildcards('dir/file_*.fa', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['dir/file.fa']
        list2 *. stageName == ['dir/file_1.fa', 'dir/file_2.fa', 'dir/file_3.fa']

        when:
        list1 = processor.expandWildcards('dir/*', [FileHolder.get('file.fa')])
        list2 = processor.expandWildcards('dir/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/file.fa']
        list2 *. stageName == ['dir/titi.fa', 'dir/toto.fa']

        when:
        list1 = processor.expandWildcards('dir/*/*', [FileHolder.get('file.fa')])
        list2 = processor.expandWildcards('dir/*/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/1/file.fa']
        list2 *. stageName == ['dir/1/titi.fa', 'dir/2/toto.fa']

        when:
        list1 = processor.expandWildcards('dir/foo*/*', [FileHolder.get('file.fa')])
        list2 = processor.expandWildcards('dir/foo*/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/foo1/file.fa']
        list2 *. stageName == ['dir/foo1/titi.fa', 'dir/foo2/toto.fa']

        when:
        list1 = processor.expandWildcards('dir/??/*', [FileHolder.get('file.fa')])
        list2 = processor.expandWildcards('dir/??/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/01/file.fa']
        list2 *. stageName == ['dir/01/titi.fa', 'dir/02/toto.fa']

        when:
        list1 = processor.expandWildcards('dir/bar??/*', [FileHolder.get('file.fa')])
        list2 = processor.expandWildcards('dir/bar??/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/bar01/file.fa']
        list2 *. stageName == ['dir/bar01/titi.fa', 'dir/bar02/toto.fa']
    }

    @Unroll
    def 'should expand wildcards rule' () {

        given:
        def processor = [:] as TaskProcessor

        expect:
        processor.expandWildcards0(pattern, 'stage-name.txt', index, size ) == expected

        where:
        pattern             | index | size  | expected
        // just wildcard
        '*'                 | 1     | 1     | 'stage-name.txt'
        '*'                 | 1     | 10    | 'stage-name.txt'
        // wildcard on the file name and single item in the collection
        'foo.txt'           | 1     | 1     | 'foo.txt'
        'foo*.fa'           | 1     | 1     | 'foo.fa'
        'foo?.fa'           | 1     | 1     | 'foo1.fa'
        'foo??.fa'          | 1     | 1     | 'foo01.fa'
        // wildcard on the file name and many items in the collection
        'foo*.fa'           | 1     | 3     | 'foo1.fa'
        'foo*.fa'           | 3     | 3     | 'foo3.fa'
        'foo?.fa'           | 1     | 3     | 'foo1.fa'
        'foo?.fa'           | 3     | 3     | 'foo3.fa'
        'foo??.fa'          | 1     | 3     | 'foo01.fa'
        'foo??.fa'          | 3     | 3     | 'foo03.fa'
        // wildcard on parent path
        'dir/*/foo.txt'     | 1     | 1     | 'dir/1/foo.txt'
        'dir/foo*/bar.txt'  | 1     | 1     | 'dir/foo1/bar.txt'
        'dir/foo?/bar.txt'  | 2     | 2     | 'dir/foo2/bar.txt'
        'dir/foo??/bar.txt' | 2     | 2     | 'dir/foo02/bar.txt'
        // wildcard on parent path and name
        'dir/*/'            | 1     | 1     | 'dir/1/stage-name.txt'
        'dir/*/*'           | 1     | 1     | 'dir/1/stage-name.txt'
        'dir/*/*'           | 1     | 10    | 'dir/1/stage-name.txt'
        'dir/*/foo*.txt'    | 1     | 1     | 'dir/1/foo.txt'
        'dir/*/foo*.txt'    | 1     | 2     | 'dir/1/foo1.txt'
        'dir/*/foo?.txt'    | 2     | 2     | 'dir/2/foo2.txt'
        'dir/???/foo?.txt'  | 5     | 10    | 'dir/005/foo5.txt'
    }

    @Unroll
    def 'should replace question marks' () {
        given:
        def processor = [:] as TaskProcessor

        expect:
        processor.replaceQuestionMarkWildcards(pattern, index) == expected

        where:
        pattern         | index | expected
        'foo.txt'       | 1     | 'foo.txt'
        'foo?.txt'      | 1     | 'foo1.txt'
        'foo???.txt'    | 2     | 'foo002.txt'
        'foo?_???.txt'  | 3     | 'foo3_003.txt'
        'foo??.txt'     | 9999  | 'foo9999.txt'

    }

    def 'should stage path'() {

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
        holder.storePath == localFile.toRealPath()
        holder.stageName == localFile.getFileName().toString()

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
        proc.checkErrorStrategy(task, error, 1, 1, 0) == ErrorStrategy.TERMINATE

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy: 'ignore')
        then:
        proc.checkErrorStrategy(task, error, 10, 10, 0) == ErrorStrategy.IGNORE

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy: 'finish')
        then:
        proc.checkErrorStrategy(task, error, 1, 1, 0) == ErrorStrategy.FINISH

    }

    def 'should return TERMINATE or FINISH error strategy`' () {
        given:
        def task
        def proc = [:] as TaskProcessor
        def error = Mock(ProcessUnrecoverableException)

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy: 'retry')
        then:
        proc.checkErrorStrategy(task, error, 1, 1, 0) == ErrorStrategy.TERMINATE

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy: 'ignore')
        then:
        proc.checkErrorStrategy(task, error, 1, 1, 0) == ErrorStrategy.TERMINATE

        when:
        task = new TaskRun()
        task.config = new TaskConfig(errorStrategy: 'finish')
        then:
        proc.checkErrorStrategy(task, error, 1, 1, 0) == ErrorStrategy.FINISH

    }

    @Unroll
    def 'should return `retry` strategy' () {

        given:

        def task
        def error = Mock(ProcessException)
        def session = Mock(Session)
        session.getExecService() >> Mock(ExecutorService)

        def proc = [:] as TaskProcessor
        proc.session = session

        when:
        task = new TaskRun(context: new TaskContext(holder: [:]))
        task.config = new TaskConfig(errorStrategy: 'retry', maxErrors: MAX_ERRORS, maxRetries: MAX_RETRIES )
        then:
        proc.checkErrorStrategy(task, error, TASK_ERR_COUNT , PROC_ERR_COUNT, SUBMIT_RETRIES) == EXPECTED

        where:
        MAX_RETRIES | MAX_ERRORS    |   TASK_ERR_COUNT  |  PROC_ERR_COUNT   | SUBMIT_RETRIES    | EXPECTED
                1   |        3      |               0   |               0   | 0                 | ErrorStrategy.RETRY
                1   |        3      |               1   |               0   | 0                 | ErrorStrategy.RETRY
                1   |        3      |               2   |               0   | 0                 | ErrorStrategy.TERMINATE
                1   |        3      |               0   |               1   | 0                 | ErrorStrategy.RETRY
                1   |        3      |               0   |               2   | 0                 | ErrorStrategy.RETRY
                1   |        3      |               0   |               3   | 0                 | ErrorStrategy.TERMINATE
                3   |       -1      |               0   |               0   | 0                 | ErrorStrategy.RETRY
                3   |       -1      |               1   |               1   | 0                 | ErrorStrategy.RETRY
                3   |       -1      |               2   |               2   | 0                 | ErrorStrategy.RETRY
                3   |       -1      |               3   |               9   | 0                 | ErrorStrategy.RETRY
                3   |       -1      |               4   |               9   | 0                 | ErrorStrategy.TERMINATE
         and:
         // terminates when the submit retries is greater than the max retries
                1   |       -1      |               0   |               0   | 1                 | ErrorStrategy.RETRY
                1   |       -1      |               0   |               0   | 2                 | ErrorStrategy.TERMINATE
                3   |       -1      |               0   |               0   | 2                 | ErrorStrategy.RETRY
                3   |       -1      |               0   |               0   | 2                 | ErrorStrategy.RETRY
                3   |       -1      |               0   |               0   | 4                 | ErrorStrategy.TERMINATE

    }


    private List<String> fetchResultFiles(TaskProcessor processor, FileOutParam param, String namePattern, Path folder ) {
        processor
                .fetchResultFiles(param, namePattern, folder)
                .collect { Path it -> folder.relativize(it).toString() }
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
        result = fetchResultFiles(processor, Mock(FileOutParam), '*.fa', folder )
        then:
        result == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setType('file')
        result = fetchResultFiles(processor, param, '*.fa', folder)
        then:
        result  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setType('dir')
        result = fetchResultFiles(processor, param, '*.fa', folder)
        then:
        result == []

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = fetchResultFiles(processor, param, '**.fa', folder)
        then:
        result == ['dir1/dir2/file4.fa', 'dir_link/dir2/file4.fa', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setFollowLinks(false)
        result = fetchResultFiles(processor, param, '**.fa', folder)
        then:
        result == ['dir1/dir2/file4.fa', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setMaxDepth(1)
        result = fetchResultFiles(processor, param, '**.fa', folder)
        then:
        result == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = fetchResultFiles(processor, param, '*', folder)
        then:
        result == ['dir1', 'dir_link', 'file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setType('dir')
        result = fetchResultFiles(processor, param, '*', folder)
        then:
        result == ['dir1', 'dir_link']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setType('file')
        result = fetchResultFiles(processor, param, '*', folder)
        then:
        result == ['file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setType('file')
        param.setHidden(true)
        result = fetchResultFiles(processor, param, '*', folder)
        then:
        result == ['.hidden.fa', 'file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = fetchResultFiles(processor, param,'.*', folder)
        then:
        result == ['.hidden.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = fetchResultFiles(processor, param,'file{1,2}.{txt,fa}', folder)
        then:
        result == ['file1.txt', 'file2.fa']

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
        param.setType('dir')
        then:
        processor.visitOptions(param,'dir-name') == [type:'dir', followLinks: true, maxDepth: null, hidden: false, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setHidden(true)
        then:
        processor.visitOptions(param,'dir-name') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setFollowLinks(false)
        then:
        processor.visitOptions(param,'dir-name') == [type:'any', followLinks: false, maxDepth: null, hidden: false, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.setMaxDepth(5)
        then:
        processor.visitOptions(param,'dir-name') == [type:'any', followLinks: true, maxDepth: 5, hidden: false, relative: false]
    }

    def 'should get bin files in the script command' () {

        given:
        def session = Mock(Session)
        session.getBinEntries() >> ['foo.sh': Paths.get('/some/path/foo.sh'), 'bar.sh': Paths.get('/some/path/bar.sh')]
        def processor = [:] as TaskProcessor
        processor.session = session

        when:
        def result = processor.getTaskBinEntries('var=x foo.sh')
        then:
        result.size()==1
        result.contains(Paths.get('/some/path/foo.sh'))

        when:
        result = processor.getTaskBinEntries('echo $(foo.sh); bar.sh')
        then:
        result.size()==2
        result.contains(Paths.get('/some/path/foo.sh'))
        result.contains(Paths.get('/some/path/bar.sh'))

    }

    def 'should make task unique id' () {

        given:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            getBinEntries() >> ['foo.sh': Paths.get('/some/path/foo.sh'), 'bar.sh': Paths.get('/some/path/bar.sh')]
        }
        and:
        def task = Mock(TaskRun) {
            getSource() >> 'hello world'
            isContainerEnabled() >> false
            getContainer() >> null
            getConfig() >> Mock(TaskConfig)
        }
        and:
        def processor = Spy(TaskProcessor)
        processor.@session = session
        processor.@config = Mock(ProcessConfig)

        when:
        def uuid1 = processor.createTaskHashKey(task)
        def uuid2 = processor.createTaskHashKey(task)
        then:
        // global var should *not* change task hash
        processor.getTaskGlobalVars(task) >>> [
                [foo:'a', bar:'b'],
                [bar:'b', foo:'a']
        ]
        and:
        uuid1 == uuid2

    }

    def 'should export env vars' () {

        given:
        def env

        when:
        env = TaskProcessor.bashEnvironmentScript([FOO:'hola',BAR:'ciao mondo'])
        then:
        env ==  '''
                export FOO="hola"
                export BAR="ciao mondo"
                '''
                .stripIndent().leftTrim()

        when:
        env = TaskProcessor.bashEnvironmentScript([PATH: 'foo:$PATH', HOLA: 'one|two'], true)
        then:
        env == '''\
            export PATH="foo:\\$PATH"
            export HOLA="one|two"
            '''.stripIndent()
        env.charAt(env.size()-1) == '\n' as char

        when:
        env = TaskProcessor.bashEnvironmentScript([FOO:null, BAR:''])
        then:
        env == "export FOO=''\nexport BAR=''\n"

    }

    def 'should normalise to path' () {
        given:
        def proc = new TaskProcessor()

        expect:
        proc.normalizeToPath('/foo/bar') == '/foo/bar' as Path
        and:
        proc.normalizeToPath('file:///foo/bar') == '/foo/bar' as Path
        and:
        proc.normalizeToPath(Paths.get('foo.txt')) == Paths.get('foo.txt')

        when:
        proc.normalizeToPath('abc')
        then:
        thrown(ProcessUnrecoverableException)

        when:
        proc.normalizeToPath(null)
        then:
        thrown(ProcessUnrecoverableException)
    }

    def 'should normalize files' () {
        given:
        def batch = Mock(FilePorter.Batch)
        def executor = Mock(Executor)
        def PATH = Paths.get('/some/path')
        def proc = new TaskProcessor(); proc.executor = executor

        when:
        def result = proc.normalizeInputToFiles(PATH.toString(), 0, true, batch)
        then:
        1 * executor.isForeignFile(PATH) >> false
        0 * batch.addToForeign(PATH) >> null
        result.size() == 1
        result[0] == new FileHolder(PATH)

    }

    def 'should get task directive vars' () {
        given:
        def processor = Spy(TaskProcessor)
        processor.@config = Mock(ProcessConfig)
        and:
        def task = Mock(TaskRun)
        and:
        def config = new TaskConfig()
        config.cpus = 4
        config.ext.alpha = 'AAAA'
        config.ext.delta = { foo }
        config.ext.omega = "${-> bar}"
        and:
        config.setContext( foo: 'DDDD', bar: 'OOOO' )

        when:
        def result = processor.getTaskExtensionDirectiveVars(task)
        then:
        1 * task.getVariableNames() >> {[ 'task.cpus', 'task.ext.alpha', 'task.ext.delta', 'task.ext.omega' ] as Set}
        1 * task.getConfig() >> config
        then:
        result == [
                'task.ext.alpha': 'AAAA',
                'task.ext.delta': 'DDDD',
                'task.ext.omega': 'OOOO',
        ]
    }

    def 'should bind fair outputs' () {
        given:
        def processor = Spy(TaskProcessor)
        processor.@config = Mock(ProcessConfig)
        processor.@isFair0 = true
        and:
        def emission3 = new HashMap()
        def task3 = Mock(TaskRun) { getIndex()>>3 }
        and:
        def emission2 = new HashMap()
        def task2 = Mock(TaskRun) { getIndex()>>2 }
        and:
        def emission1 = new HashMap()
        def task1 = Mock(TaskRun) { getIndex()>>1 }
        and:
        def emission5 = new HashMap()
        def task5 = Mock(TaskRun) { getIndex()>>5 }
        and:
        def emission4 = new HashMap()
        def task4 = Mock(TaskRun) { getIndex()>>4 }

        when:
        processor.fairBindOutputs0(emission3, task3)
        then:
        processor.@fairBuffers[2] == emission3
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(emission2, task2)
        then:
        processor.@fairBuffers[1] == emission2
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(emission5, task5)
        then:
        processor.@fairBuffers[4] == emission5
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(emission1, task1)
        then:
        1 * processor.bindOutputs0(emission1)
        then:
        1 * processor.bindOutputs0(emission2)
        then:
        1 * processor.bindOutputs0(emission3)
        and:
        processor.@fairBuffers.size() == 2 
        processor.@fairBuffers[0] == null
        processor.@fairBuffers[1] == emission5

        when:
        processor.fairBindOutputs0(emission4, task4)
        then:
        1 * processor.bindOutputs0(emission4)
        then:
        1 * processor.bindOutputs0(emission5)
        then:
        processor.@fairBuffers.size()==0
    }

    def 'should parse env map' () {
        given:
        def workDir = TestHelper.createInMemTempDir()
        def envFile = workDir.resolve(TaskRun.CMD_ENV)
        envFile.text =  '''
                        ALPHA=one
                        /ALPHA/
                        DELTA=x=y
                        /DELTA/
                        OMEGA=
                        /OMEGA/
                        LONG=one
                        two
                        three
                        /LONG/=exit:0
                        '''.stripIndent()
        and:
        def processor = Spy(TaskProcessor)

        when:
        def result = processor.collectOutEnvMap(workDir, Map.of())
        then:
        result == [ALPHA:'one', DELTA: "x=y", OMEGA: '', LONG: 'one\ntwo\nthree']
    }

    def 'should parse env map with command error' () {
        given:
        def workDir = TestHelper.createInMemTempDir()
        def envFile = workDir.resolve(TaskRun.CMD_ENV)
        envFile.text =  '''
                        ALPHA=one
                        /ALPHA/
                        cmd_out_1=Hola
                        /cmd_out_1/=exit:0
                        cmd_out_2=This is an error message
                        for unknown reason
                        /cmd_out_2/=exit:100
                        '''.stripIndent()
        and:
        def processor = Spy(TaskProcessor)

        when:
        processor.collectOutEnvMap(workDir, [cmd_out_1: 'foo --this', cmd_out_2: 'bar --that'])
        then:
        def e = thrown(ProcessEvalException)
        e.message == 'Unable to evaluate output'
        e.command == 'bar --that'
        e.output == 'This is an error message\nfor unknown reason'
        e.status == 100
    }
    def 'should create a task preview' () {
        given:
        def config = new ProcessConfig([cpus: 10, memory: '100 GB'])
        def EXEC = Mock(Executor) { getName()>>'exec-name'}
        def BODY = Mock(BodyDef) { getType()>>ScriptType.SCRIPTLET }
        def processor = new TaskProcessor(config: config, name: 'proc-name', executor: EXEC, taskBody: BODY)

        when:
        def result = processor.createTaskPreview()
        then:
        result.config.process == 'proc-name'
        result.config.executor == 'exec-name'
        result.config.getCpus() == 10
        result.config.getMemory() == MemoryUnit.of('100 GB')
    }

    @Unroll
    def 'should validate inputs arity' () {
        given:
        def executor = Mock(Executor)
        def session = Mock(Session) {getFilePorter()>>Mock(FilePorter) }
        def processor = Spy(new TaskProcessor(session:session, executor:executor))
        and:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig())

        when:
        def param = new FileInParam(new Binding(), [])
                .setPathQualifier(true)
                .bind(FILE_NAME) as FileInParam
        if( ARITY )
            param.setArity(ARITY)

        processor.makeTaskContextStage2(task, [(param):FILE_VALUE], 0 )
        then:
        context.get(FILE_NAME) == EXPECTED

        where:
        FILE_NAME       | FILE_VALUE                                | ARITY     | EXPECTED
        'file.txt'      | '/some/file.txt'                          | null      | Path.of('/some/file.txt')
        'file.*'        | '/some/file.txt'                          | null      | Path.of('/some/file.txt')
        'file.*'        | ['/some/file1.txt','/some/file2.txt']     | null      | [Path.of('/some/file1.txt'), Path.of('/some/file2.txt')]
        '*'             | ['/some/file1.txt','/some/file2.txt']     | null      | [Path.of('/some/file1.txt'), Path.of('/some/file2.txt')]
        '*'             | []                                        | null      | []

        and:
        'file.txt'      | '/some/file.txt'                          | '1'      | Path.of('/some/file.txt')
        'f*'            | '/some/file.txt'                          | '1'      | Path.of('/some/file.txt')
        'f*'            | '/some/file.txt'                          | '1..2'   | [Path.of('/some/file.txt')]
        'f*'            | '/some/file.txt'                          | '1..*'   | [Path.of('/some/file.txt')]
        'f*'            | '/some/file.txt'                          | '1..*'   | [Path.of('/some/file.txt')]
        'f*'            | ['/some/file.txt']                        | '1..*'   | [Path.of('/some/file.txt')]
        'f*'            | ['/some/file1.txt', '/some/file2.txt']    | '1..*'   | [Path.of('/some/file1.txt'), Path.of('/some/file2.txt')]
    }

    def 'should throw an arity error' () {
        given:
        def executor = Mock(Executor)
        def session = Mock(Session) {getFilePorter()>>Mock(FilePorter) }
        def processor = Spy(new TaskProcessor(session:session, executor:executor))
        and:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig())

        when:
        def param = new FileInParam(new Binding(), [])
                .setPathQualifier(true)
                .bind(FILE_NAME) as FileInParam
        if( ARITY )
            param.setArity(ARITY)

        processor.makeTaskContextStage2(task, [(param):FILE_VALUE], 0 )
        then:
        def e = thrown(IllegalArityException)
        e.message == ERROR

        where:
        FILE_NAME       | FILE_VALUE                                | ARITY     | ERROR
        'file.txt'      | []                                        | '0'       | 'Path arity max value must be greater or equals to 1'
        'file.txt'      | []                                        | '1'       | 'Incorrect number of input files for process `foo` -- expected 1, found 0'
        'f*'            | []                                        | '1..*'    | 'Incorrect number of input files for process `foo` -- expected 1..*, found 0'
        'f*'            | '/some/file.txt'                          | '2..*'    | 'Incorrect number of input files for process `foo` -- expected 2..*, found 1'
        'f*'            | ['/some/file.txt']                        | '2..*'    | 'Incorrect number of input files for process `foo` -- expected 2..*, found 1'
        'f*'            | ['/a','/b']                               | '3'       | 'Incorrect number of input files for process `foo` -- expected 3, found 2'
    }

    def 'should validate collect output files' () {
        given:
        def executor = Mock(Executor)
        def session = Mock(Session) {getFilePorter()>>Mock(FilePorter) }
        def processor = Spy(new TaskProcessor(session:session, executor:executor))
        and:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig())
        and:
        def workDir = Path.of('/work')

        when:
        def param = new FileOutParam(new Binding(), [])
                .setPathQualifier(true)
                .optional(OPTIONAL)
                .bind(FILE_NAME) as FileOutParam
        if( ARITY )
            param.setArity(ARITY)
        and:
        processor.collectOutFiles(task, param, workDir, context)
        then:
        processor.fetchResultFiles(_,_,_) >> RESULTS
        processor.checkFileExists(_,_) >> EXISTS
        and:
        task.getOutputs().get(param) == EXPECTED

        where:
        FILE_NAME       | RESULTS                                   | EXISTS    | OPTIONAL  | ARITY         | EXPECTED
        'file.txt'      | null                                      | true      | false     | null          | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | true      | false     | null          | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | true      | false     | null          | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | true      | true      | null          | []
        and:
        'file.txt'      | null                                      | true      | false     | '1'           | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | true      | false     | '1'           | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | true      | false     | '1..*'        | [Path.of('/work/file.txt')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | true      | false     | '2'           | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | true      | false     | '1..*'        | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | true      | false     | '0..*'        | []
    }

    @Unroll
    def 'should report output error' () {
        given:
        def executor = Mock(Executor)
        def session = Mock(Session) {getFilePorter()>>Mock(FilePorter) }
        def processor = Spy(new TaskProcessor(session:session, executor:executor))
        and:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig())
        and:
        def workDir = Path.of('/work')

        when:
        def param = new FileOutParam(new Binding(), [])
                .setPathQualifier(true)
                .optional(OPTIONAL)
                .bind(FILE_NAME) as FileOutParam
        if( ARITY )
            param.setArity(ARITY)
        and:
        processor.collectOutFiles(task, param, workDir, context)
        then:
        processor.fetchResultFiles(_,_,_) >> RESULTS
        processor.checkFileExists(_,_) >> EXISTS
        and:
        def e = thrown(EXCEPTION)
        e.message == ERROR

        where:
        FILE_NAME       | RESULTS                                   | EXISTS    | OPTIONAL  | ARITY         | EXCEPTION             | ERROR
        'file.txt'      | null                                      | false     | false     | null          | MissingFileException  | "Missing output file(s) `file.txt` expected by process `foo`"
        '*'             | []                                        | true      | false     | null          | MissingFileException  | "Missing output file(s) `*` expected by process `foo`"
        and:
        'file.txt'      | null                                      | true      | false     | '2'           | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2, found 1"
        '*'             | [Path.of('/work/file.txt')]               | true      | false     | '2'           | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2, found 1"
        '*'             | [Path.of('/work/file.txt')]               | true      | false     | '2..*'        | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2..*, found 1"
        '*'             | []                                        | true      | true      | '1..*'        | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 1..*, found 0"

    }

}
