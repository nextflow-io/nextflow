/*
 * Copyright (c) 2012, the authors.
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

package nextflow.script
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.AbstractExecutor
import nextflow.executor.LocalExecutor
import nextflow.executor.NopeExecutor
import nextflow.executor.SgeExecutor
import nextflow.executor.SlurmExecutor
import nextflow.processor.MergeTaskProcessor
import nextflow.processor.ParallelTaskProcessor
import nextflow.processor.TaskConfig
import nextflow.processor.TaskConfigWrapper
import nextflow.processor.TaskProcessor
import nextflow.util.CacheHelper
import nextflow.util.FileHelper
/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseScript extends Script {


    protected BaseScript(){ }

    protected BaseScript(Binding binding) {
        super(binding)
    }

    /*
     * The script execution session, declare it private to prevent the user script to be able to access it
     */
    @Lazy
    private Session session = { getBinding()?.getVariable('__$session') as Session } ()

    private final random = new Random()


    /**
     * Holds the configuration object which will used to execution the user tasks
     */
    @Lazy
    Map config = { session.config } ()

    @Lazy
    InputStream stdin = { System.in }()

    private TaskProcessor taskProcessor

    /** Access to the last *process* object -- only for testing purpose */
    @PackageScope
    TaskProcessor getTaskProcessor() { taskProcessor }

    private result

    /** Access to the last *process* result -- only for testing purpose */
    @PackageScope
    Object getResult() { result }

    /**
     * Enable disable task 'echo' configuration property
     * @param value
     */
    def void echo(boolean value = true) {
        config.process.echo = value
    }

    /**
     * Stop the current execution returning an error code and message
     *
     * @param exitCode The exit code to be returned
     * @param message The message that will be reported in the log file (optional)
     */
    def void exit(int exitCode, String message = null) {
        if ( exitCode && message ) {
            log.error message
        }
        else if ( message ) {
            log.info message
        }
        System.exit(exitCode)
    }

    /**
     * Stop the current execution returning a 0 error code and the specified message
     *
     * @param message The message that will be reported in the log file
     */
    def void exit( String message ) {
        exit(0, message)
    }

    /**
     * Create a folder for the given key. It guarantees to return the same folder name
     * the same provided object key.
     *
     * @param key An object to be used as cache-key creating the folder, it can be any object
     *          or an array or objects to use multi-objects key
     *
     * @return The {@code File} to the cached directory or a newly created folder for the specified key
     */
    File cacheableDir( Object key ) {
        assert key, "Please specify the 'key' argument on 'cacheableDir' method"

        def hash = CacheHelper.hasher([ session.uniqueId, key, session.cacheable ? 0 : random.nextInt() ]).hash()

        def file = FileHelper.createWorkFolder(session.workDir, hash)
        if( !file.exists() && !file.mkdirs() ) {
            throw new IOException("Unable to create folder: $file -- Check file system permission" )
        }

        return file
    }

    /**
     * Create a file for the given key. It guarantees to return the same file name
     * the same provided object key.
     *
     * @param key
     * @param name
     * @return
     */
    File cacheableFile( Object key, String name = null ) {

        // the cacheability is guaranteed by the folder
        def folder = cacheableDir(key)

        if( !name ) {
            name = key instanceof File ? key.name : key.toString()
        }

        return new File(folder, name)
    }

    /**
     * @return Create a temporary directory
     */
    File tempDir() {
        def file = FileHelper.createTempFolder(session.workDir)
        if( !file.exists() && !file.mkdirs() ) {
            throw new IOException("Unable to create folder: $file -- Check file system permission" )
        }
        return file
    }

    /**
     * @return Create a temporary file
     */
    File tempFile( String name = 'temp.file') {

        def folder = tempDir()
        return new File(folder, name)
    }

    /**
     * Create a task processor
     *
     * @param name The name used to label this processor
     * @param block The code block to be executed
     * @return The {@code Processor} instance
     */
    private createProcessor( Class<? extends TaskProcessor> processorClass, String name, Closure<String> block  ) {
        assert block

        def taskConfig = new TaskConfig(this)

        // set 'default' properties defined in the configuration file in the 'task' section
        if( config.process instanceof Map ) {
            config.process .each { String key, value -> taskConfig.setProperty(key,value) }
        }
        else if( config.task instanceof Map ) {
            log.warn "Note: configuration attribute 'task' has been deprecated -- replace it by using attribute 'process'"
            config.task.each { String key, value -> taskConfig.setProperty(key,value) }
        }

        // set the task name in the config object
        if( name ) {
            taskConfig.name = name
        }

        // Invoke the code block, which will return the script closure to the executed
        // As side effect will set all the properties declaration in the 'taskConfig' object
        // Note: the config object is wrapped by a TaskConfigWrapper because it is required
        // to raise a MissingPropertyException when some values is missing, so that the Closure
        // will try to fallback on the owner object
        def script = new TaskConfigWrapper(taskConfig).with ( block ) as Closure
        if ( !script ) throw new IllegalArgumentException("Missing script in the specified task block -- make sure it terminates with the script string to be executed")

        // load the executor to be used
        def executorClass = loadExecutorClass( getExecutorConfig(taskConfig) )

        def executor = executorClass.newInstance()
        def processor = processorClass.newInstance( executor, session, this, taskConfig, script )

        return taskProcessor = processor

    }

    /**
     * Find out the 'executor' to be used in the process definition or in teh session configuration object
     *
     * @param taskConfig
     */
    private getExecutorConfig(Map taskConfig) {
        log.debug ">> taskConfig $taskConfig"

        // create the processor object
        def result = taskConfig.executor?.toString()

        // fallback on deprecated attribute
        if( !result && taskConfig.processor ) {
            result = taskConfig.processor?.toString()
            log.warn "Note: configuration attribute 'processor' has been deprecated -- replace it by using the attribute 'executor'"
        }

        // fallback on config file definition
        if( !result ) {
            result = session.config.executor?.toString()
        }

        if( !result && session.config.processor ) {
            result = session.config.processor?.toString()
            log.warn "Note: configuration attribute 'processor' has been deprecated -- replace it by using the attribute 'executor' in the 'nextflow.conf' file"
        }

        log.debug "<< taskConfig result: $result"
        return result
    }


    def process( Map<String,?> args, String name, Closure code ) {
        log.debug "Create task: $name -- $args "
        def clazz = args.merge ? MergeTaskProcessor : ParallelTaskProcessor
        result = createProcessor(clazz, name, code).run()
    }

    def process( String name, Closure code ) {
        log.debug "Create task: $name -- noargs"
        result = createProcessor(ParallelTaskProcessor, name, code).run()
    }



//    def process(Map<String,?> args, String name ) {
//        log.debug "Create task: $name with: $args"
//        task(name, { return '' })
//    }

//    def process(String name ) {
//        log.debug "Create task: $name "
//        task(name, { return '' })
//    }


    /*
     * Map the executor class to its 'friendly' name
     */
    static executorsMap = [
            'local': LocalExecutor.name,
            'sge':  SgeExecutor.name,
            'oge':  SgeExecutor.name,
            'slurm': SlurmExecutor.name,
            'nope': NopeExecutor.name
    ]

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




}
