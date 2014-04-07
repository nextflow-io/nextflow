/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.AbstractExecutor
import nextflow.executor.LocalExecutor
import nextflow.executor.LsfExecutor
import nextflow.executor.NopeExecutor
import nextflow.executor.SgeExecutor
import nextflow.executor.SlurmExecutor
import nextflow.executor.SupportedScriptTypes
import nextflow.script.BaseScript
import nextflow.script.ScriptType
import nextflow.script.TaskBody

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessFactory {


    /*
   * Map the executor class to its 'friendly' name
   */
    static executorsMap = [
            'nope': NopeExecutor.name,
            'local': LocalExecutor.name,
            'sge':  SgeExecutor.name,
            'oge':  SgeExecutor.name,
            'lsf': LsfExecutor.name,
            'slurm': SlurmExecutor.name,
            'dnanexus': 'nextflow.executor.DnaNexusExecutor',
            'hazelcast': 'nextflow.executor.HzExecutor',
            'gg': 'nextflow.executor.GgExecutor',
            'gridgain': 'nextflow.executor.GgExecutor'
    ]

    private final Session session

    private final Map config

    private BaseScript owner

    ProcessFactory( BaseScript ownerScript, Session session, Map config = null ) {
        this.owner = ownerScript
        this.session = session
        this.config = config != null ? config : session.config
    }

    @PackageScope
    static Class<? extends AbstractExecutor> loadExecutorClass(String executorName) {
        log.debug ">> processorType: $executorName"
        def className = executorName ? executorsMap[ executorName?.toLowerCase()  ] : LocalExecutor.name

        // if the className is empty (because the 'processorType' does not map to any class, fallback to the 'processorType' itself)
        if( !className ) {
            className = executorName
        }

        log.debug "Loading executor class: ${className}"
        try {
            Thread.currentThread().getContextClassLoader().loadClass(className) as Class<AbstractExecutor>
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Cannot find a valid class for specified executor: '${executorName}'")
        }

    }


    @PackageScope
    static boolean isTypeSupported( ScriptType type, executor ) {

        if( executor instanceof AbstractExecutor ) {
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
     * @param name The name of the process as defined in the script
     * @param body The process declarations provided by the user
     * @return The {@code Processor} instance
     */
    TaskProcessor createProcessor( String name, Closure body, Map options = null ) {
        assert body

        /*
         * check if exists 'attributes' defined in the 'process' scope for this process, e.g.
         *
         * process.$name.attribute1 = xxx
         * process.$name.attribute2 = yyy
         *
         */
        Map importantAttributes = null
        if( config.process instanceof Map && config.process['$'+name] instanceof Map ) {
            importantAttributes = (Map)config.process['$'+name]
        }

        // -- the config object
        def taskConfig = new TaskConfig(owner, importantAttributes)

        // -- set 'default' properties defined in the configuration file in the 'process' section
        if( config.process instanceof Map ) {
            config.process .each { String key, value ->
                if( key.startsWith('$')) return
                taskConfig.setProperty(key,value)
            }
        }

        /*
         * options that may be passed along the process declaration e.g.
         *
         * process name ( x: 1, ... ) {
         *
         * }
         */

        options?.each { String key, value -> taskConfig.setProperty(key,value)}

        // -- set the task name in the config object
        if( name ) {
            taskConfig.name = name
        }

        // Invoke the code block, which will return the script closure to the executed
        // As side effect will set all the properties declaration in the 'taskConfig' object
        // Note: the config object is wrapped by a TaskConfigWrapper because it is required
        // to raise a MissingPropertyException when some values is missing, so that the Closure
        // will try to fallback on the owner object
        def script = new TaskConfigWrapper(taskConfig).with ( body ) as TaskBody
        if ( !script )
            throw new IllegalArgumentException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")

        // load the executor to be used
        def execName = getExecutorName(taskConfig) ?: 'local'
        def execClass = loadExecutorClass(execName)

        if( !isTypeSupported(script.type, execClass) ) {
            log.warn "Process '$name' cannot be executed by '$execName' executor -- Using 'local' executor instead"
            execName = 'local'
            execClass = LocalExecutor.class
        }

        def execObj = execClass.newInstance()
        // inject the task configuration into the executor instance
        execObj.taskConfig = taskConfig
        execObj.session = session
        execObj.name = execName
        execObj.init()

        // create processor class
        def processorClass = taskConfig.merge ? MergeTaskProcessor : ParallelTaskProcessor
        def result = processorClass.newInstance( execObj, session, owner, taskConfig, script )
        return result

    }

    /**
     * Find out the 'executor' to be used in the process definition or in teh session configuration object
     *
     * @param taskConfig
     */
    private getExecutorName(Map taskConfig) {
        log.trace ">> taskConfig $taskConfig"

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
