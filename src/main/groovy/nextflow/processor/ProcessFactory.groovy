/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
        for( Class<Executor> exec : ServiceDiscover.load(Executor).iterator() ) {
            log.debug "Discovered executor class: ${exec.name}"
            executorsMap.put(findNameByClass(exec), exec.name)
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
    def static String findNameByClass( Class<Executor> clazz ) {
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

        def clazz =  executorsMap[executorName.toLowerCase().replace('-','')]
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
