/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.script


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import nextflow.processor.TaskProcessor
/**
 *  Factory class for {@TaskProcessor} instances
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessFactory {

    private Session session

    private Map config

    private BaseScript owner

    private ExecutorFactory executorFactory

    /* only for test -- do not use */
    protected ProcessFactory() { }

    ProcessFactory( BaseScript ownerScript, Session session ) {
        this.owner = ownerScript
        this.session = session
        this.config = session.config
        this.executorFactory = session.executorFactory
    }

    /**
     * Create a new task processor and initialise with the given parameters
     *
     * @param name The processor name
     * @param executor The executor object
     * @param session The session object
     * @param script The owner script
     * @param config The process configuration
     * @param taskBody The process task body
     * @return An instance of {@link nextflow.processor.TaskProcessor}
     */
    protected TaskProcessor newTaskProcessor(String name, Executor executor, ProcessConfig config, BodyDef taskBody ) {
        new TaskProcessor(name, executor, session, owner, config, taskBody)
    }

    /**
     * Create a task processor
     *
     * @param name
     *      The name of the process as defined in the script
     * @param body
     *      The process declarations provided by the user
     * @param options
     *      A map representing the named parameter specified after the process name eg:
     *          `process foo(bar: 'x') {  }`
     *      (not used)
     * @return
     *      The {@code Processor} instance
     */
    TaskProcessor createProcessor( String name, Closure<BodyDef> body ) {
        assert body
        assert config.process instanceof Map

        // -- the config object
        final processConfig = new ProcessConfig(owner, name)
        // Invoke the code block which will return the script closure to the executed.
        // As side effect will set all the property declarations in the 'taskConfig' object.
        processConfig.throwExceptionOnMissingProperty(true)
        final copy = (Closure<BodyDef>)body.clone()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.setDelegate(processConfig)
        final script = copy.call()
        processConfig.throwExceptionOnMissingProperty(false)
        if ( !script )
            throw new IllegalArgumentException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")

        // -- apply settings from config file to process config
        processConfig.applyConfigLegacy((Map)config.process, name)

        // -- get the executor for the given process config
        final execObj = executorFactory.getExecutor(name, processConfig, script, session)

        // -- create processor class
        newTaskProcessor( name, execObj, processConfig, script )
    }

}
