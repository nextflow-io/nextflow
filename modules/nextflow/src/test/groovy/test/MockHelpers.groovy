package test

import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.ScriptRunner
import nextflow.script.ScriptType

class MockScriptRunner extends ScriptRunner {

    MockScriptRunner() {
        super(new MockSession())
    }

    MockScriptRunner setScript(String str) {
        def script = TestHelper.createInMemTempFile('main.nf', str)
        setScript(script)
        return this
    }

}

class MockSession extends Session {

    @Override
    Session start() {
        this.executorFactory = new MockExecutorFactory()
        return super.start()
    }

    MockSession() {
        super()
    }

}

class MockExecutorFactory extends ExecutorFactory {
    @Override
    protected Class<? extends Executor> getExecutorClass(String executorName) {
        return MockExecutor
    }

    @Override
    protected boolean isTypeSupported(ScriptType type, Object executor) {
        true
    }
}

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MockExecutor extends Executor {

    @Override
    void signal() { }

    protected TaskMonitor createTaskMonitor() {
        new MockMonitor()
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new  MockTaskHandler(task)
    }
}

class MockMonitor implements TaskMonitor {

    void schedule(TaskHandler handler) {
        handler.submit()
    }

    /**
     * Remove the {@code TaskHandler} instance from the queue of tasks to be processed
     *
     * @param handler A not null {@code TaskHandler} instance
     */
    boolean evict(TaskHandler handler) { }

    /**
     * Start the monitoring activity for the queued tasks
     * @return The instance itself, useful to chain methods invocation
     */
    TaskMonitor start() { }

    /**
     * Notify when a task terminates
     */
    void signal() { }
}

@Slf4j
class MockTaskHandler extends TaskHandler {

    protected MockTaskHandler(TaskRun task) {
        super(task)
    }

    @Override
    void submit() {
        log.info ">> launching mock task: ${task}"
        if( task.type == ScriptType.SCRIPTLET ) {
            task.workDir = Paths.get('.').complete()
            task.stdout = task.script
            task.exitStatus = 0
        }
        else {
            task.code.call()
        }
        status = TaskStatus.COMPLETED
        task.processor.finalizeTask(task)
    }

    @Override
    boolean checkIfRunning() {
        return false
    }

    @Override
    boolean checkIfCompleted() {
        true
    }

    @Override
    void kill() { }

}
