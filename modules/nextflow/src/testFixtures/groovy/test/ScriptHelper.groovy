/*
 * Copyright 2013-2025, Seqera Labs
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

package test

import java.nio.file.Path

import groovy.transform.InheritConstructors
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.Session
import nextflow.config.ConfigParserFactory
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ChannelOut
import nextflow.script.ProcessConfig
import nextflow.script.ProcessFactory
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import nextflow.script.ScriptLoaderFactory
import nextflow.script.ScriptType

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptHelper {

    /**
     * Load a config from source text.
     *
     * @param text
     */
    static Map loadConfig(String text) {
        return ConfigParserFactory.create().parse(text)
    }

    /**
     * Load a script from source text.
     *
     * This function compiles and executes the script without
     * running the pipeline. The entry workflow is executed
     * (i.e. to construct the workflow DAG) unless opts.module
     * is enabled.
     *
     * @param opts
     * @param text
     */
    static BaseScript loadScript(Map opts = [:], String text) {
        def session = opts.config ? new MockSession(opts.config as Map) : new MockSession()
        session.setBinding(new ScriptBinding())

        session.init(null)
        session.start()

        def loader = ScriptLoaderFactory.create(session)
        if( opts.module )
            loader.setModule(true)
        loader.parse(text)
        loader.runScript()

        return loader.getScript()
    }

    /**
     * Load a script from source text and return the TaskProcessor
     * for the last parsed process.
     *
     * The script should declare a single process.
     *
     * @param scriptText
     * @param binding
     */
    static TaskProcessor parseAndReturnProcess(String scriptText, Map binding = null) {
        def session = new TestSession()
        session.setBinding(binding ? new ScriptBinding(binding) : new ScriptBinding())

        session.init(null)

        ScriptLoaderFactory.create(session)
            .parse(scriptText)
            .runScript()

        session.fireDataflowNetwork()

        return TaskProcessor.currentProcessor()
    }

    @InheritConstructors
    private static class TestSession extends Session {

        TestProcessFactory processFactory

        @Override
        ProcessFactory newProcessFactory(BaseScript script) {
            return processFactory = new TestProcessFactory(script, this)
        }
    }

    private static class TestProcessFactory extends ProcessFactory {

        BaseScript script
        Session session

        TestProcessFactory(BaseScript script, Session session) {
            super(script, session)
            this.script = script
            this.session = session
        }

        @Override
        TaskProcessor newTaskProcessor(String name, Executor executor, ProcessConfig config, BodyDef taskBody) {
            new TestTaskProcessor(name, executor, session, script, config, taskBody)
        }
    }

    @InheritConstructors
    private static class TestTaskProcessor extends TaskProcessor {
        @Override
        def run () {
            // this is needed to mimic the out channels normalisation
            // made by the real 'run' method - check the superclass
            if( config.getOutputs().size() == 0 )
                config.fakeOutput()
        }
    }

    /**
     * Execute a script from source text, with an optional config map.
     *
     * This function compiles and executes the script, launches a
     * pipeline run, and waits for the run to finish. It returns the
     * last statement of the entry workflow, which can be used to pass
     * output channels to a test.
     *
     * @param text
     * @param config
     */
    static Object runScript(String text, Map config = null) {
        def session = config ? new MockSession(config) : new MockSession()
        session.setBinding(new ScriptBinding())

        session.init(null)
        session.start()

        def loader = ScriptLoaderFactory.create(session)
        loader.parse(text)
        loader.runScript()

        def result = normalizeResult(loader.getResult())

        session.fireDataflowNetwork()
        session.await()
        session.destroy()

        return result
    }

    /**
     * Execute a script from a source file, with an optional config map.
     *
     * This function compiles and executes the script, launches a
     * pipeline run, and waits for the run to finish. It returns the
     * last statement of the entry workflow, which can be used to pass
     * output channels to a test.
     *
     * @param opts
     * @param path
     */
    static Object runScript(Map opts = [:], Path path) {
        def session = opts.config ? new MockSession(opts.config) : new MockSession()
        session.setBinding(new ScriptBinding())

        session.init( new ScriptFile(path) )
        session.start()

        def loader = ScriptLoaderFactory.create(session)
        loader.parse(path)
        loader.runScript()

        def result = normalizeResult(loader.getResult())

        session.fireDataflowNetwork()
        session.await()
        session.destroy()

        return result
    }

    /**
     * Run a block of code that executes a dataflow network.
     *
     * This method ensures that the dataflow network is properly
     * ignited, and that any dataflow result is normalized into a
     * readable channel.
     *
     * @param action
     */
    static Object runDataflow(Map opts = [:], Closure action) {
        def session = opts.config ? new MockSession(opts.config) : new MockSession()

        def cl = (Closure)action.clone()
        cl.setDelegate(session.binding)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        
        def result = normalizeResult(cl.call())

        session.fireDataflowNetwork()

        return result
    }

    private static Object normalizeResult(value) {
        if( value instanceof ChannelOut ) {
            def result = new ArrayList<>(value.size())
            for( def el : value )
                result.add(normalizeResult0(el))
            return result.size() == 1 ? result[0] : result
        }

        if( value instanceof Object[] || value instanceof List ) {
            def result = new ArrayList<>(value.size())
            for( def el : value )
                result.add(normalizeResult0(el))
            return result
        }

        return normalizeResult0(value)
    }

    private static Object normalizeResult0(value) {
        if( value instanceof ChannelOut ) {
            def result = new ArrayList<>(value.size())
            for( def el : value )
                result.add(normalizeResult0(el))
            return result.size() == 1 ? result[0] : result
        }

        if( value instanceof DataflowBroadcast )
            return value.createReadChannel()

        return value
    }

}

class MockSession extends Session {

    MockSession(Map config) {
        super(config)
    }

    MockSession() {
        super()
    }

    @Override
    Session start() {
        this.executorFactory = new MockExecutorFactory()
        return super.start()
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

class MockExecutor extends Executor {

    @Override
    void signal() {}

    protected TaskMonitor createTaskMonitor() {
        new MockMonitor()
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new MockTaskHandler(task)
    }
}

class MockMonitor implements TaskMonitor {

    @Override void schedule(TaskHandler handler) { handler.submit() }
    @Override boolean evict(TaskHandler handler) {}
    @Override TaskMonitor start() {}
    @Override void signal() {}
}

class MockTaskHandler extends TaskHandler {

    protected MockTaskHandler(TaskRun task) {
        super(task)
    }

    @Override
    void submit() {
        if( task.type == ScriptType.SCRIPTLET ) {
            task.workDir = Path.of('.').complete()
            task.stdout = task.script
            task.exitStatus = 0
        }
        else {
            task.code.call()
        }
        status = TaskStatus.COMPLETED
        task.processor.finalizeTask(this)
    }

    @Override
    boolean checkIfRunning() { false }

    @Override
    boolean checkIfCompleted() { true }

    @Override
    protected void killTask() {}
}
