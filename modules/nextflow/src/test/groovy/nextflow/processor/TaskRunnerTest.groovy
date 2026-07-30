/*
 * Copyright 2013-2026, Seqera Labs
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
import java.util.concurrent.ExecutorService
import java.util.function.Consumer

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.cache.CacheDB
import nextflow.exception.IllegalArityException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.file.FilePorter
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfigV1
import nextflow.script.ScriptBinding
import nextflow.script.ScriptType
import nextflow.script.params.FileInParam
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskRunnerTest extends Specification {

    private TaskRunner newRunner(TaskProcessor processor) {
        new TaskRunner(processor, { TaskRun task -> } as Consumer<TaskRun>)
    }

    def 'should execute a task without a dataflow network' () {
        given:
        def workDir = Files.createTempDirectory('test')
        def script = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        def body = new BodyDef({ 'echo hello' }, "echo hello")
        and:
        def executor = Mock(Executor) {
            getName() >> 'nope'
            getWorkDir() >> workDir
            getStageDir() >> workDir.resolve('stage')
        }
        def session = Mock(Session) {
            getUniqueId() >> UUID.randomUUID()
            getBinEntries() >> [:]
            getCache() >> Mock(CacheDB) { getTaskEntry(_,_) >> null }
            getFilePorter() >> Mock(FilePorter) { newBatch(_) >> Mock(FilePorter.Batch) }
        }
        and:
        // the process is *not* started, ie. no dataflow operator and no output channels
        def processor = new TaskProcessor(
                name: 'foo',
                session: session,
                executor: executor,
                ownerScript: script,
                taskBody: body,
                config: new ProcessConfigV1(script, null) )
        and:
        def completed = []
        def runner = new TaskRunner(processor, { TaskRun task -> completed << task } as Consumer<TaskRun>)

        when:
        runner.submit(new TaskStartParams(TaskId.of(1), 1), [])

        then:
        1 * executor.submit({ TaskRun task ->
            task.name == 'foo (1)' && task.script == 'echo hello' && task.workDir.startsWith(workDir)
        })
        and:
        processor.getOperator() == null
        completed.isEmpty()

        cleanup:
        workDir?.deleteDir()
    }

    def 'should return `ignore` strategy' () {

        given:
        def task
        def proc = newRunner([:] as TaskProcessor)
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
        def proc = newRunner([:] as TaskProcessor)
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

        def processor = [:] as TaskProcessor
        processor.session = session
        def proc = newRunner(processor)

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

    @Unroll
    def 'should apply input file arity' () {
        given:
        def executor = Mock(Executor)
        executor.isForeignFile(_) >> false
        def session = Mock(Session)
        def config = new ProcessConfigV1(Mock(BaseScript), null)
        def runner = newRunner(new TaskProcessor(session:session, executor:executor, config:config))
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

        runner.resolveTaskInputs(task, [FILE_VALUE], foreignFiles )
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
        def runner = newRunner(new TaskProcessor(session:session, executor:executor, config:config))
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

        runner.resolveTaskInputs(task, [FILE_VALUE], foreignFiles)
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

    def 'should submit a task' () {
        given:
        def exec = Mock(Executor)
        def proc = newRunner(new TaskProcessor(executor: exec))
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
        def proc = newRunner(new TaskProcessor(executor: exec, arrayCollector: collector))
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
}
