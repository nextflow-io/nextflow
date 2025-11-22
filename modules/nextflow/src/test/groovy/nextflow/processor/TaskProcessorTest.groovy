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

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ExecutorService

import com.google.common.hash.HashCode
import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.exception.IllegalArityException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.executor.NopeExecutor
import nextflow.file.FileHolder
import nextflow.file.FilePorter
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ProcessConfigV1
import nextflow.script.ScriptType
import nextflow.script.bundle.ResourcesBundle
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {

    def createProcessor(String name, Session session) {
        return new DummyProcessor(name, session, Mock(BaseScript), new ProcessConfig([:]))
    }

    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(String name, Session session, BaseScript script, ProcessConfig config) {
            super(name, new NopeExecutor(session: session), session, script, config, new BodyDef({}, '..'))
        }

        @Override
        protected void createOperator() { }
    }


    def 'should evaluate environment variables'() {

        setup:
        def home = Files.createTempDirectory('test')
        def binFolder = home.resolve('bin')
        binFolder.mkdirs()

        when:
        def session = new Session([env: [X:"1", Y:"2"]])
        session.setBaseDir(home)
        def processor = createProcessor('task1', session)
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
        processor = createProcessor('task1', session)
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
        def task3 = Mock(TaskRun) { getIndex()>>3 }
        and:
        def task2 = Mock(TaskRun) { getIndex()>>2 }
        and:
        def task1 = Mock(TaskRun) { getIndex()>>1 }
        and:
        def task5 = Mock(TaskRun) { getIndex()>>5 }
        and:
        def task4 = Mock(TaskRun) { getIndex()>>4 }

        when:
        processor.fairBindOutputs0(task3)
        then:
        processor.@fairBuffers[2] == task3
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(task2)
        then:
        processor.@fairBuffers[1] == task2
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(task5)
        then:
        processor.@fairBuffers[4] == task5
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(task1)
        then:
        1 * processor.bindOutputs0(task1)
        then:
        1 * processor.bindOutputs0(task2)
        then:
        1 * processor.bindOutputs0(task3)
        and:
        processor.@fairBuffers.size() == 2 
        processor.@fairBuffers[0] == null
        processor.@fairBuffers[1] == task5

        when:
        processor.fairBindOutputs0(task4)
        then:
        1 * processor.bindOutputs0(task4)
        then:
        1 * processor.bindOutputs0(task5)
        then:
        processor.@fairBuffers.size()==0
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
    def 'should apply input file arity' () {
        given:
        def executor = Mock(Executor)
        executor.isForeignFile(_) >> false
        def session = Mock(Session)
        def config = new ProcessConfigV1(Mock(BaseScript), null)
        def processor = Spy(new TaskProcessor(session:session, executor:executor, config:config))
        def foreignFiles = Mock(FilePorter.Batch)
        and:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig())

        when:
        def param = new FileInParam(config)
                .setPathQualifier(true)
                .bind(FILE_NAME)
        if( ARITY )
            param.setArity(ARITY)
        and:
        task.setInput(param)

        processor.resolveTaskInputs(task, [FILE_VALUE], foreignFiles )
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

    def 'should report input file arity error' () {
        given:
        def executor = Mock(Executor)
        executor.isForeignFile(_) >> false
        def session = Mock(Session)
        def config = new ProcessConfigV1(Mock(BaseScript), null)
        def processor = Spy(new TaskProcessor(session:session, executor:executor, config:config))
        def foreignFiles = Mock(FilePorter.Batch)
        and:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig())

        when:
        def param = new FileInParam(config)
                .setPathQualifier(true)
                .bind(FILE_NAME)
        if( ARITY )
            param.setArity(ARITY)
        and:
        task.setInput(param)

        processor.resolveTaskInputs(task, [FILE_VALUE], foreignFiles)
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

    def 'should collect output files' () {
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
        processor.collectOutFiles(task, param, workDir)
        then:
        processor.collectOutFiles0(_,_,_) >> RESULTS
        and:
        task.getOutputs().get(param) == EXPECTED

        where:
        FILE_NAME       | RESULTS                                   | OPTIONAL  | ARITY         | EXPECTED
        'file.txt'      | [Path.of('/work/file.txt')]               | false     | null          | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | false     | null          | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | false     | null          | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | true      | null          | []
        and:
        'file.txt'      | [Path.of('/work/file.txt')]               | false     | '1'           | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | false     | '1'           | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | false     | '1..*'        | [Path.of('/work/file.txt')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | false     | '2'           | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | false     | '1..*'        | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | false     | '0..*'        | []
    }

    @Unroll
    def 'should report output file arity error' () {
        given:
        def executor = Mock(Executor)
        def session = Mock(Session)
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
        processor.collectOutFiles(task, param, workDir)
        then:
        processor.collectOutFiles0(_,_,_) >> RESULTS
        and:
        def e = thrown(EXCEPTION)
        e.message == ERROR

        where:
        FILE_NAME       | RESULTS                                   | OPTIONAL  | ARITY         | EXCEPTION             | ERROR
        'file.txt'      | [Path.of('/work/file.txt')]               | false     | '2'           | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2, found 1"
        '*'             | [Path.of('/work/file.txt')]               | false     | '2'           | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2, found 1"
        '*'             | [Path.of('/work/file.txt')]               | false     | '2..*'        | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2..*, found 1"
        '*'             | []                                        | true      | '1..*'        | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 1..*, found 0"

    }

    def 'should submit a task' () {
        given:
        def exec = Mock(Executor)
        def proc = Spy(new TaskProcessor(executor: exec))
        and:
        def task = Mock(TaskRun)
        def hash = HashCode.fromString('0123456789abcdef')
        def workDir = Path.of('/work')

        when:
        proc.submitTask(task, hash, workDir)
        then:
        task.getConfig() >> new TaskConfig()
        and:
        1 * exec.submit(task)
    }

    def 'should collect a task' () {
        given:
        def exec = Mock(Executor)
        def collector = Mock(TaskArrayCollector)
        def proc = Spy(new TaskProcessor(executor: exec, arrayCollector: collector))
        and:
        def task = Mock(TaskRun)
        def hash = HashCode.fromString('0123456789abcdef')
        def workDir = Path.of('/work')

        when:
        proc.submitTask(task, hash, workDir)
        then:
        task.getConfig() >> new TaskConfig()
        and:
        1 * collector.collect(task)
        0 * exec.submit(task)

        when:
        proc.submitTask(task, hash, workDir)
        then:
        task.getConfig() >> new TaskConfig(attempt: 2)
        and:
        0 * collector.collect(task)
        1 * exec.submit(task)
    }

    def 'should compute eval outputs content deterministically'() {

        setup:
        def processor = createProcessor('test', Mock(Session))

        when:
        def result1 = processor.computeEvalOutputsContent([
            'nxf_out_eval_2': 'echo "value2"',
            'nxf_out_eval_1': 'echo "value1"',
            'nxf_out_eval_3': 'echo "value3"'
        ])
        
        def result2 = processor.computeEvalOutputsContent([
            'nxf_out_eval_3': 'echo "value3"',
            'nxf_out_eval_1': 'echo "value1"',
            'nxf_out_eval_2': 'echo "value2"'
        ])

        then:
        result1 == result2
        result1 == 'nxf_out_eval_1=echo "value1"\nnxf_out_eval_2=echo "value2"\nnxf_out_eval_3=echo "value3"'
    }
}
