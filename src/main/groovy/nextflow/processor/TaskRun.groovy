/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessException
import nextflow.exception.ProcessMissingTemplateException
import nextflow.exception.ProcessNotRecoverableException
import nextflow.file.FileHolder
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.InParam
import nextflow.script.OutParam
import nextflow.script.ScriptType
import nextflow.script.StdInParam
import nextflow.script.TaskBody
import nextflow.script.ValueOutParam
import nextflow.container.ContainerScriptTokens
import nextflow.container.DockerBuilder
import nextflow.container.ShifterBuilder

/**
 * Models a task instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class TaskRun implements Cloneable {

    /**
     * Task unique id
     */
    def id

    /**
     * Task index within its execution group
     */
    def index

    /**
     * Task name
     */
    def String name

    /**
     * The unique hash code associated to this task
     */
    def HashCode hash

    /*
     * The processor that creates this 'task'
     */
    TaskProcessor processor

    /**
     * Holds the input value(s) for each task input parameter
     */
    Map<InParam,Object> inputs = [:]

    /**
     * Holds the output value(s) for each task output parameter
     */
    Map<OutParam,Object> outputs = [:]


    def void setInput( InParam param, Object value = null ) {
        assert param

        inputs[param] = value

        // copy the value to the task 'input' attribute
        // it will be used to pipe it to the process stdin
        if( param instanceof StdInParam) {
            stdin = value
        }
    }

    def void setOutput( OutParam param, Object value = null ) {
        assert param
        outputs[param] = value
    }


    /**
     * The value to be piped to the process stdin
     */
    def stdin

    /**
     * The exit code returned by executing the task script
     */
    Integer exitStatus = Integer.MAX_VALUE

    /**
     * Flag set when the bind stage is completed successfully
     */
    boolean canBind

    /**
     * Set {@code true} when the task executed is resumed from the cache
     */
    boolean cached

    /**
     * Task produced standard output
     */
    def stdout

    /**
     * Task produced standard error
     */
    def stderr

    /**
     * @return The task produced stdout result as string
     */
    def String getStdout() {

        if( stdout instanceof Path ) {
            try {
                return stdout.text
            }
            catch( NoSuchFileException e ) {
                return null
            }
        }

        if( stdout != null ) {
            return stdout.toString()
        }

        return null
    }

    def String getStderr() {

        if( stderr instanceof Path ) {
            try {
                return stderr.text
            }
            catch( NoSuchFileException e ) {
                return null
            }
        }

        if( stderr != null ) {
            return stderr.toString()
        }

        return null
    }

    /**
     * Print to the current system console the task produced stdout
     */
    void echoStdout() {

        // print the stdout
        if( stdout instanceof Path ) {
            try {
                synchronized (System.out) {
                    stdout.withReader {  reader ->
                        reader.eachLine { System.out.println(it) }
                    }
                }
            }
            catch( NoSuchFileException e ) {
                log.trace "Echo file does not exist: ${stdout}"
            }
            return
        }

        if( stdout != null ) {
            print stdout.toString()
        }

    }

    /**
     * Dump the first {@code max} line of the task {@code #stdout}
     *
     * @param message The message buffer that will hold the first lines
     * @param n The maximum number of lines to dump (default: 20)
     * @return The actual number of dumped lines into {@code message} buffer
     */
    List<String> dumpStdout(int n = 50) {

        try {
            return dumpObject(stdout,n)
        }
        catch( Exception e ) {
            log.debug "Unable to dump output of process '$name' -- Cause: ${e}"
            return []
        }
    }

    List<String> dumpStderr(int n = 50) {

        try {
            return dumpObject(stderr,n)
        }
        catch( Exception e ) {
            log.debug "Unable to dump error of process '$name' -- Cause: ${e}"
            return []
        }
    }


    protected List<String> dumpObject( obj, int max ) {
        List result = null

        if( obj instanceof Path )
            result = obj.tail(max).readLines()

        else if( obj != null ) {
            result = obj.toString().readLines()
            if( result.size()>max ) {
                result = result[-max..-1]
                result.add(0, '(more omitted..)')
            }
        }

        return result ?: []
    }

    /**
     * The directory used to run the task
     */
    Path workDir

    /**
     * The type of the task code: native or external system command
     */
    ScriptType type

    /**
     * The runtime exception when execution groovy native code
     */
    Throwable error

    /*
     * The closure implementing this task
     */
    Closure code

    /**
     * The script to be executed by the system
     */
    def script

    /**
     * The process source as entered by the user in the process definition
     */
    String source

    /**
     * The process template file
     */
    Path template

    /**
     * The number of times the execution of the task has failed
     */
    def volatile int failCount

    /**
     * Mark the task as failed
     */
    def volatile boolean failed

    def volatile boolean aborted

    def TaskConfig config

    def TaskContext context

    TaskProcessor.RunType runType = TaskProcessor.RunType.SUBMIT

    TaskRun clone() {
        final taskClone = (TaskRun)super.clone()
        taskClone.context = context.clone()
        taskClone.config = config.clone()
        taskClone.config.setContext(taskClone.context)
        return taskClone
    }

    String getName() {
        if( name )
            return name

        final baseName = processor.name
        if( config.containsKey('tag') )
            try {
                // -- look-up the 'sampleId' property, and if everything is fine
                //    cache this value in the 'name' attribute
                return name = "$baseName (${config.tag})"
            }
            catch( IllegalStateException e ) {
                log.debug "Cannot access `tag` property for task: $baseName ($index)"
            }

        // fallback on the current task index, however do not set the 'name' attribute
        // so it has a chance to recover the 'sampleId' at next invocation
        return "$baseName ($index)"
    }

    String getScript() {
        if( script instanceof Path ) {
            return script.text
        }
        else {
            return script?.toString()
        }
    }

    /**
     * Check whenever there are values to be cached
     */
    boolean hasCacheableValues() {

        if( config?.isDynamic() )
            return true

        for( OutParam it : outputs.keySet() ) {
            if( it.class == ValueOutParam ) return true
            if( it.class == FileOutParam && ((FileOutParam)it).isDynamic() ) return true
        }

        return false
    }

    Map<InParam,List<FileHolder>> getInputFiles() {
        (Map<InParam,List<FileHolder>>) getInputsByType( FileInParam )
    }

    /**
     * Return the list of all input files staged as inputs by this task execution
     */
    List<String> getStagedInputs()  {
        getInputFiles()
                .values()
                .flatten()
                .collect { it.stageName }
    }

    /**
     * @return A map object containing all the task input files as <stage name, store path> pairs
     */
    Map<String,Path> getInputFilesMap() {

        def result = [:]
        def allFiles = getInputFiles().values()
        for( List<FileHolder> entry : allFiles ) {
            if( entry ) for( FileHolder it : entry ) {
                result[ it.stageName ] = it.storePath
            }
        }

        return result
    }

    /**
     * It is made of two parts:
     *
     * 1) Look at the {@code nextflow.script.FileOutParam} which name is the expected
     *  output name
     *
     * 2) It looks shared file parameters, that being so are also output parameters
     *  The problem here is that we need the real file name as it has been staged in
     *  process execution folder.
     */
    @Memoized
    def List<String> getOutputFilesNames() {
        def result = []

        getOutputsByType(FileOutParam).keySet().each { FileOutParam param ->
            result.addAll( param.getFilePatterns(context, workDir) )
        }

        return result.unique()
    }


    /**
     * Get the map of *input* objects by the given {@code InParam} type
     *
     * @param types One ore more subclass of {@code InParam}
     * @return An associative array containing all the objects for the specified type
     */
    def <T extends InParam> Map<T,Object> getInputsByType( Class<T>... types ) {

        def result = [:]
        inputs.findAll() { types.contains(it.key.class) }.each { result << it }
        return result
    }

    /**
     * Get the map of *output* objects by the given {@code InParam} type
     *
     * @param types One ore more subclass of {@code InParam}
     * @return An associative array containing all the objects for the specified type
     */
    public <T extends OutParam> Map<T,Object> getOutputsByType( Class<T>... types ) {
        def result = [:]
        outputs.each {
            if( types.contains(it.key.class) ) {
                result << it
            }
        }
        return result
    }

    /**
     * @return A map containing the task environment defined as input declaration by this task
     */
    Map<String,String> getInputEnvironment() {
        final Map<String,String> environment = [:]
        getInputsByType( EnvInParam ).each { param, value ->
            environment.put( param.name, value?.toString() )
        }
        return environment
    }


    Path getTargetDir() {
        config.getStoreDir() ?: workDir
    }

    def getScratch() {
        config.scratch
    }

    static final public String CMD_LOG = '.command.log'
    static final public String CMD_ENV = '.command.env'
    static final public String CMD_SCRIPT = '.command.sh'
    static final public String CMD_INFILE = '.command.in'
    static final public String CMD_OUTFILE = '.command.out'
    static final public String CMD_ERRFILE = '.command.err'
    static final public String CMD_EXIT = '.exitcode'
    static final public String CMD_START = '.command.begin'
    static final public String CMD_RUN = '.command.run'
    static final public String CMD_STUB = '.command.run.1'
    static final public String CMD_CONTEXT = '.command.val'
    static final public String CMD_TRACE = '.command.trace'



    String toString( ) {

        "id: $id; name: $name; type: $type; status: $exitStatus; error: $error; workDirectory: $workDir"

    }


    /**
     * Given an {@code HashCode} instance renders only the first 8 characters using the format showed below:
     *      <pre>
     *          2f/0075a6
     *     </pre>
     *
     *
     * @param hash An {@code HashCode} object
     * @return The short representation of the specified hash code as string
     */
    String getHashLog() {
        if( !hash ) return null
        def str = hash.toString()
        def result = new StringBuilder()
        result << str[0]
        result << str[1]
        result << '/'
        for( int i=2; i<8 && i<str.size(); i++ ) {
            result << str[i]
        }
        return result.toString()
    }

    /**
     * The name of a docker container where the task is supposed to run when provided
     */
    String getContainer() {
        // set the docker container to be used
        String imageName
        if( isContainerExecutable() ) {
            imageName = ContainerScriptTokens.parse(script.toString()).image
        }
        else {
            imageName = config.container as String
        }

        if( getShifterConfig().enabled ) {
            ShifterBuilder.normalizeImageName(imageName, getShifterConfig())
        }
        else {
            DockerBuilder.normalizeImageName(imageName, getDockerConfig())
        }
    }

    /**
     * @return A {@link Map} object holding the configuration attributes for the Docker engine
     */
    Map getDockerConfig() {
        def result = processor.getSession().config?.docker as Map
        result != null ? result : Collections.emptyMap()
    }

    Map getShifterConfig() {
        def result = processor.getSession().config?.shifter as Map
        result != null ? result : Collections.emptyMap()
    }

    /**
     * @return {@true} when the process must run within a container and the docker engine is enabled
     */
    boolean isDockerEnabled() {
        if( isContainerExecutable() )
            return true

        return getDockerConfig()?.enabled?.toString() == 'true'
    }

    /**
     * @return {@true} when the process runs an *executable* container
     */
    boolean isContainerExecutable() {
        config.container == true
    }

    boolean isSuccess( status = exitStatus ) {
        if( status == null )
            return false

        if( status instanceof String )
            status = status.toInteger()

        if( status == Integer.MAX_VALUE )
            return false

        return status in config.getValidExitStatus()
    }

    /**
     * Given the task body token declared in the process definition:
     * 1) assign to the task run instance the code closure to execute
     * 2) extract the process code `source`
     * 3) assign the `script` code to execute
     *
     * @param body A {@code TaskBody} object instance
     */
    @PackageScope void resolve( TaskBody body ) {

        // -- initialize the task code to be executed
        this.code = body.closure.clone() as Closure
        this.code.delegate = this.context
        this.code.setResolveStrategy(Closure.DELEGATE_ONLY)

        // -- set the task source
        // note: this may be overwritten when a template file is used
        this.source = body.source

        if( body.type != ScriptType.SCRIPTLET )
            return

        // Important!
        // when the task is implemented by a script string
        // Invoke the closure which returns the script with all the variables replaced with the actual values
        try {
            def result = code.call()
            if ( result instanceof Path ) {
                script = renderTemplate(result, body.isShell)
            }
            else if( result != null && body.isShell ) {
                script = renderScript(result)
            }
            else {
                // note `toString` transform GString to plain string object resolving
                // interpolated variables
                script = result.toString()
            }
        }
        catch( ProcessException e ) {
            throw e
        }
        catch( Throwable e ) {
            throw new ProcessNotRecoverableException("Process `$name` script contains error(s)", e)
        }

    }


    /**
     * Given a template script file and a binding, returns the rendered script content
     *
     * @param script
     * @param binding
     * @return
     */
    final protected String renderTemplate( Path template, boolean shell=false ) {
        try {
            // read the template and override the task `source` with it
            this.source = template.text
            // keep track of template file
            this.template = template
            // parse the template
            final engine = new TaskTemplateEngine(processor.grengine)
            if( shell ) {
                engine.setPlaceholder(placeholderChar())
            }
            return engine.render(source, context)
        }
        catch( NoSuchFileException e ) {
            throw new ProcessMissingTemplateException("Process `${processor.name}` can't find template file: $template")
        }
    }

    final protected String renderScript( script ) {

        new TaskTemplateEngine(processor.grengine)
                .setPlaceholder(placeholderChar())
                .setEnableShortNotation(false)
                .render(script.toString(), context)
    }

    protected placeholderChar() {
        (config.placeholder ?: '!') as char
    }

}

