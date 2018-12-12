/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.aws.batch.AwsBatchExecutor
import nextflow.executor.CondorExecutor
import nextflow.executor.CrgExecutor
import nextflow.executor.Executor
import nextflow.executor.LocalExecutor
import nextflow.executor.LsfExecutor
import nextflow.executor.NopeExecutor
import nextflow.executor.NqsiiExecutor
import nextflow.executor.PbsExecutor
import nextflow.executor.PbsProExecutor
import nextflow.executor.SgeExecutor
import nextflow.executor.SlurmExecutor
import nextflow.executor.SupportedScriptTypes
import nextflow.k8s.K8sExecutor
import nextflow.script.BaseScript
import nextflow.script.ScriptType
import nextflow.script.TaskBody
import nextflow.util.ServiceDiscover
import nextflow.util.ServiceName
/**
 *  Factory class for {@TaskProcessor} instances
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessFactory {

    static public String DEFAULT_EXECUTOR = System.getenv('NXF_EXECUTOR') ?: 'local'

    /*
     * Map the executor class to its 'friendly' name
     */
    final protected Map executorsMap = [
            'nope': NopeExecutor,
            'local': LocalExecutor,
            'sge':  SgeExecutor,
            'oge':  SgeExecutor,
            'uge':  SgeExecutor,
            'lsf': LsfExecutor,
            'pbs': PbsExecutor,
            'pbspro': PbsProExecutor,
            'slurm': SlurmExecutor,
            'crg': CrgExecutor,
            'bsc': LsfExecutor,
            'condor': CondorExecutor,
            'k8s': K8sExecutor,
            'nqsii': NqsiiExecutor,
            'awsbatch': AwsBatchExecutor
    ]

    private final Session session

    private final Map config

    private BaseScript owner

    /* only for test -- do not use */
    protected ProcessFactory() {

    }

    ProcessFactory( BaseScript ownerScript, Session session, Map config = null ) {
        this.owner = ownerScript
        this.session = session
        this.config = config != null ? config : session.config

        // discover non-core executors
        for( Class<Executor> clazz : ServiceDiscover.load(Executor) ) {
            log.trace "Discovered executor class: ${clazz.name}"
            executorsMap.put(findNameByClass(clazz), clazz)
        }
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
     * @return An instance of {@link TaskProcessor}
     */
    protected TaskProcessor newTaskProcessor( String name, Executor executor, Session session, BaseScript script, ProcessConfig config, TaskBody taskBody ) {
        new TaskProcessor(name, executor, session, script, config, taskBody)
    }

    /**
     * Extract the executor name by using the annotation {@code ServiceName} or fallback to simple classname
     * if the annotation is not provided
     *
     * @param clazz
     * @return
     */
    static String findNameByClass( Class<Executor> clazz ) {
        def annotation = clazz.getAnnotation(ServiceName)
        if( annotation )
            return annotation.value()

        def name = clazz.getSimpleName().toLowerCase()
        if( name.endsWith('executor') ) {
            name = name.subSequence(0, name.size()-'executor'.length())
        }

        return name
    }

    protected Class<? extends Executor> loadExecutorClass(String executorName) {
        log.debug ">> processorType: '$executorName'"
        if( !executorName )
            return LocalExecutor

        def clazz =  executorsMap[executorName.toLowerCase()]
        if( !clazz )
            throw new IllegalArgumentException("Unknown executor name: $executorName")

        if( clazz instanceof Class )
            return clazz

        if( !(clazz instanceof String ) )
            throw new IllegalArgumentException("Not a valid executor class object: $clazz")

        // if the className is empty (because the 'processorType' does not map to any class, fallback to the 'processorType' itself)
        if( !clazz ) {
            clazz = executorName
        }

        log.debug "Loading executor class: ${clazz}"
        try {
            Thread.currentThread().getContextClassLoader().loadClass(clazz as String) as Class<Executor>
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Cannot find a valid class for specified executor: '${executorName}'")
        }

    }


    protected boolean isTypeSupported( ScriptType type, executor ) {

        if( executor instanceof Executor ) {
            executor = executor.class
        }

        if( executor instanceof Class ) {
            def annotation = executor.getAnnotation(SupportedScriptTypes)
            if( !annotation )
                throw new IllegalArgumentException("Specified argument is not a valid executor class: $executor -- Missing 'SupportedScriptTypes' annotation")

            return type in annotation.value()
        }

        throw new IllegalArgumentException("Specified argument is not a valid executor class: $executor")
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
    TaskProcessor createProcessor( String name, Closure body, Map options = null ) {
        assert body
        assert config.process instanceof Map

        /*
         * check if exists 'attributes' defined in the 'process' scope for this process, e.g.
         *
         * process.$name.attribute1 = xxx
         * process.$name.attribute2 = yyy
         *
         * NOTE: THIS HAS BEEN DEPRECATED AND WILL BE REMOVED
         */
        Map legacySettings = null
        if( config.process['$'+name] instanceof Map ) {
            legacySettings = (Map)config.process['$'+name]
            log.warn "Process configuration syntax \$processName has been deprecated -- Replace `process.\$$name = <value>` with a process selector"
        }

        // -- the config object
        final processConfig = new ProcessConfig(owner).setProcessName(name)

        // Invoke the code block which will return the script closure to the executed.
        // As side effect will set all the property declarations in the 'taskConfig' object.
        processConfig.throwExceptionOnMissingProperty(true)
        final copy = (Closure)body.clone()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.setDelegate(processConfig)
        final script = copy.call() as TaskBody
        processConfig.throwExceptionOnMissingProperty(false)
        if ( !script )
            throw new IllegalArgumentException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")

        // -- Apply the directives defined in the config object using the`withLabel:` syntax
        final processLabels = processConfig.getLabels() ?: ['']
        for( String lbl : processLabels ) {
            processConfig.applyConfigForLabel(config.process as Map, "withLabel:", lbl)
        }

        // -- apply setting for name
        processConfig.applyConfigForLabel(config.process as Map, "withName:", name)

        // -- Apply process specific setting defined using `process.$name` syntax
        //    NOTE: this is deprecated and will be removed
        if( legacySettings ) {
            processConfig.applyConfigSettings(legacySettings)
        }

        // -- Apply defaults
        processConfig.applyConfigDefaults( config.process as Map )

        // -- check for conflicting settings
        if( processConfig.scratch && processConfig.stageInMode == 'rellink' ) {
            log.warn("Directives `scratch` and `stageInMode=rellink` conflict each other -- Enforcing default stageInMode for process `$name`")
            processConfig.remove('stageInMode')
        }

        // -- load the executor to be used
        def execName = getExecutorName(processConfig) ?: DEFAULT_EXECUTOR
        def execClass = loadExecutorClass(execName)

        if( !isTypeSupported(script.type, execClass) ) {
            log.warn "Process '$name' cannot be executed by '$execName' executor -- Using 'local' executor instead"
            execName = 'local'
            execClass = LocalExecutor.class
        }

        def execObj = execClass.newInstance()
        // -- inject the task configuration into the executor instance
        execObj.session = session
        execObj.name = execName
        execObj.init()

        // -- create processor class
        newTaskProcessor( name, execObj, session, owner, processConfig, script )
    }


    /**
     * Find out the 'executor' to be used in the process definition or in teh session configuration object
     *
     * @param taskConfig
     */
    private getExecutorName(ProcessConfig taskConfig) {
        // create the processor object
        def result = taskConfig.executor?.toString()

        if( !result ) {
            if( session.config.executor instanceof String ) {
                result = session.config.executor
            }
            else if( session.config.executor?.name instanceof String ) {
                result = session.config.executor.name
            }
        }

        log.debug "<< taskConfig executor: $result"
        return result
    }
}
