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
import java.nio.file.Paths
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
import nextflow.script.ProcessConfig
import nextflow.script.ProcessConfigV1
import nextflow.script.ScriptBinding
import nextflow.script.ScriptType
import nextflow.script.params.FileInParam
import nextflow.trace.TraceRecord
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

    /** A task hashable by the default {@link TaskHasher} without touching the file system. */
    private TaskRun hashableTask() {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            getBinEntries() >> [:]
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getOwnerScript() >> Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        }
        return Mock(TaskRun) {
            getSource() >> 'hello world'
            isContainerEnabled() >> false
            getConfig() >> Mock(TaskConfig)
            getProcessor() >> processor
            getGlobalVars(_) >> [:]
            getVariableNames() >> ([] as Set)
        }
    }

    def 'uses the hasher of the first TaskHasherFactory that answers, and TaskHasher when all abstain'() {
        given:
        def runner = newRunner(new TaskProcessor(session: Mock(Session), executor: Mock(Executor)))
        def task = hashableTask()
        def custom = Mock(TaskHasher)
        def abstaining = Mock(TaskHasherFactory) { create(_) >> null }
        def answering = Mock(TaskHasherFactory) { create(task) >> custom }

        when: 'factories are asked in order; the first non-null answer wins'
        runner.hasherFactories = [abstaining, answering]
        then:
        runner.createTaskHasher(task).is(custom)

        when: 'every factory abstains'
        runner.hasherFactories = [abstaining]
        def hasher = runner.createTaskHasher(task)
        then: 'the default hasher is used'
        hasher.getClass() == TaskHasher

        when: 'no factory is registered at all'
        runner.hasherFactories = []
        then:
        runner.createTaskHasher(task).getClass() == TaskHasher
    }

    def 'resolves the hasher factories lazily and defaults to TaskHasher when none is registered'() {
        given: 'a runner that has not hashed anything yet, in a JVM with no plugin system initialised'
        def runner = newRunner(new TaskProcessor(session: Mock(Session), executor: Mock(Executor)))
        def task = hashableTask()

        expect:
        runner.hasherFactories == null

        when:
        def hasher = runner.createTaskHasher(task)
        then: 'the registry was consulted once and yielded the default'
        runner.hasherFactories != null
        hasher.getClass() == TaskHasher
    }

    def 'should create a task run with the auto resource labels' () {
        given:
        def config = new ProcessConfig([resourceLabels: [team: 'genomics', 'nextflow.io/runName': 'custom']])
        def EXEC = Mock(Executor) { getName()>>'exec-name'}
        def BODY = Mock(BodyDef) { getType()>>ScriptType.SCRIPTLET }
        def SESS = Mock(Session) { getAutoResourceLabels() >> ['nextflow.io/runName': 'crazy_darwin', 'nextflow.io/sessionId': '1a2b3c'] }
        def runner = newRunner(new TaskProcessor(config: config, name: 'proc-name', executor: EXEC, taskBody: BODY, session: SESS))

        when:
        def task = runner.createTaskRun(new TaskStartParams(TaskId.of(1), 1))
        then:
        // the auto labels reach the task config, and the declared label still wins
        task.config.getResourceLabels() == [
                'nextflow.io/sessionId': '1a2b3c',
                'nextflow.io/runName': 'custom',
                team: 'genomics' ]
    }

    def 'should create a task run with no auto resource labels when the feature is off' () {
        given:
        def config = new ProcessConfig([resourceLabels: [team: 'genomics']])
        def EXEC = Mock(Executor) { getName()>>'exec-name'}
        def BODY = Mock(BodyDef) { getType()>>ScriptType.SCRIPTLET }
        def SESS = Mock(Session) { getAutoResourceLabels() >> [:] }
        def runner = newRunner(new TaskProcessor(config: config, name: 'proc-name', executor: EXEC, taskBody: BODY, session: SESS))

        when:
        def task = runner.createTaskRun(new TaskStartParams(TaskId.of(1), 1))
        then:
        task.config.getResourceLabels() == [team: 'genomics']
     }

    def 'dispatches the task resolution to the first enabled TaskCacheStrategy'() {
        given:
        def session = Mock(Session)
        def runner = newRunner(new TaskProcessor(session: session, executor: Mock(Executor)))
        def task = Mock(TaskRun)
        def hash = HashCode.fromInt(1)
        def disabled = Mock(TaskCacheStrategy)
        def enabled = Mock(TaskCacheStrategy)

        when: 'strategies are asked in priority order; the first one enabled for the session wins'
        runner.cacheStrategies = [disabled, enabled]
        runner.checkCachedOrLaunchTask(task, hash, true)
        then:
        1 * disabled.isEnabled(session) >> false
        1 * enabled.isEnabled(session) >> true
        1 * enabled.resolve(task, hash, true, { it instanceof TaskResolver })
        0 * disabled.resolve(*_)

        when: 'the choice is made once per runner'
        runner.checkCachedOrLaunchTask(task, hash, false)
        then:
        0 * _.isEnabled(_)
        1 * enabled.resolve(task, hash, false, runner.getTaskResolver())
    }

    def 'uses the default strategy when no registered strategy is enabled, and resolves them lazily'() {
        given: 'a runner that has not resolved anything yet, in a JVM with no plugin system initialised'
        def runner = newRunner(new TaskProcessor(session: Mock(Session), executor: Mock(Executor)))

        expect:
        runner.cacheStrategies == null

        when:
        def strategy = runner.getCacheStrategy()
        then: 'the registry was consulted once and yielded the default'
        runner.cacheStrategies != null
        strategy instanceof DefaultTaskCacheStrategy
        runner.getCacheStrategy().is(strategy)

        when: 'every registered strategy abstains'
        def other = newRunner(new TaskProcessor(session: Mock(Session), executor: Mock(Executor)))
        other.cacheStrategies = [ Mock(TaskCacheStrategy) { isEnabled(_) >> false } ]
        then:
        other.getCacheStrategy() instanceof DefaultTaskCacheStrategy
    }

    def 'the task resolver adapts the primitives of the runner'() {
        given:
        def cache = Mock(CacheDB)
        def session = Mock(Session) { getCache() >> cache }
        def exec = Mock(Executor) { getWorkDir() >> Paths.get('/work') }
        def processor = new TaskProcessor(session: session, executor: exec)
        def resolver = newRunner(processor).getTaskResolver()
        def hash = HashCode.fromString('0123456789abcdef')
        def entry = new TaskEntry(Mock(TraceRecord), null)
        def task = Mock(TaskRun) { getConfig() >> new TaskConfig() }

        when: 'an entry lookup'
        def found = resolver.entry(hash)
        then: 'goes to the session cache, on behalf of this processor'
        1 * cache.getTaskEntry(hash, processor) >> entry
        found.is(entry)

        expect: 'the work dir of a hash is the one the executor work root maps it to'
        resolver.workDirFor(hash) == Paths.get('/work/01/23456789abcdef')

        when: 'a launch'
        resolver.launch(task, hash, Paths.get('/work/01/23456789abcdef'))
        then: 'submits the task'
        1 * exec.submit(task)

        when: 'a resume of a work dir with no exit file'
        def folder = Files.createTempDirectory('wd')
        def copy = new TaskRun(type: ScriptType.SCRIPTLET, config: new TaskConfig(), name: 'foo')
        def cloned = 0
        def original = Mock(TaskRun) { clone() >> { cloned++; copy } }
        def resumed = resolver.resume(original, hash, folder, entry)
        then: 'is the cached-output check on a COPY of the task (the original stays launchable), which reports it as not resumable'
        cloned == 1
        !resumed

        cleanup:
        folder?.deleteDir()
    }
}
