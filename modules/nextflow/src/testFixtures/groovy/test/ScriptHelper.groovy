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

    static Map loadConfig(String text) {
        return ConfigParserFactory.create().parse(text)
    }

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

    static TaskProcessor parseAndReturnProcess(String scriptText, Map binding = null) {
        def session = new TestSession()
        session.setBinding(binding ? new ScriptBinding(binding) : new ScriptBinding())

        session.init(null)
        // session.start()

        ScriptLoaderFactory.create(session)
            .parse(scriptText)
            .runScript()

        session.fireDataflowNetwork()
        // session.await()
        // session.destroy()

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
