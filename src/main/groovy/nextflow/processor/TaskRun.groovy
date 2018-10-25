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

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.conda.CondaCache
import nextflow.conda.CondaConfig
import nextflow.container.ContainerConfig
import nextflow.container.ContainerHandler
import nextflow.container.ContainerScriptTokens
import nextflow.exception.ProcessException
import nextflow.exception.ProcessTemplateException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.file.FileHelper
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
    TaskId id

    /**
     * Task index within its execution group
     */
    def index

    /**
     * Task name
     */
    String name

    /**
     * The unique hash code associated to this task
     */
    HashCode hash

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


    void setInput( InParam param, Object value = null ) {
        assert param

        inputs[param] = value

        // copy the value to the task 'input' attribute
        // it will be used to pipe it to the process stdin
        if( param instanceof StdInParam) {
            stdin = value
        }
    }

    void setOutput( OutParam param, Object value = null ) {
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
    String getStdout() {

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

    String getStderr() {

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
                Files.copy(stdout, System.out)
            }
            catch( NoSuchFileException e ) {
                log.trace "Echo file does not exist: ${stdout}"
            }
            catch( Exception e ) {
                log.error "Unable to echo process output -- check the log file for details", e
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

    List<String> dumpLogFile(int n = 50) {
        try {
            return dumpObject(workDir.resolve(CMD_LOG),n)
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
    volatile int failCount

    /**
     * Mark the task as failed
     */
    volatile boolean failed

    /**
     * Mark the task as aborted
     */
    volatile boolean aborted

    /**
     * The action {@link ErrorStrategy} action applied if task has failed
     */
    volatile ErrorStrategy errorAction

    TaskConfig config

    TaskContext context

    TaskProcessor.RunType runType = TaskProcessor.RunType.SUBMIT

    TaskRun clone() {
        final taskClone = (TaskRun)super.clone()
        taskClone.context = context.clone()
        taskClone.config = config.clone()
        taskClone.config.setContext(taskClone.context)
        return taskClone
    }

    TaskRun makeCopy() {
        def copy = this.clone()
        // -- reset the error condition (if any)
        copy.id = TaskId.next()
        copy.name = null // <-- force to re-evaluate the name that can include a dynamic tag
        copy.error = null
        copy.exitStatus = Integer.MAX_VALUE
        return copy
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
        return processor.singleton ? baseName : "$baseName ($index)"
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
     * Look at the {@code nextflow.script.FileOutParam} which name is the expected
     *  output name
     *
     */
    @Memoized
    List<String> getOutputFilesNames() {
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
    def <T extends OutParam> Map<T,Object> getOutputsByType( Class<T>... types ) {
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
    protected Map<String,String> getInputEnvironment() {
        final Map<String,String> environment = [:]
        getInputsByType( EnvInParam ).each { param, value ->
            environment.put( param.name, value?.toString() )
        }
        return environment
    }

    /**
     * @return A map representing the task execution environment
     */
    Map<String,String> getEnvironment() {
        // note: create a copy of the process environment to avoid concurrent
        // process executions override each others
        // IMPORTANT: when copying the environment map a LinkedHashMap must be used to preserve
        // the insertion order of the env entries (ie. export FOO=1; export BAR=$FOO)
        final result = new LinkedHashMap( getProcessor().getProcessEnvironment() )
        result.putAll( getInputEnvironment() )
        return result
    }

    String getEnvironmentStr() {
        def env = getEnvironment()
        if( !env ) return null
        def result = new StringBuilder()
        env.each { k,v -> result.append(k).append('=').append(v).append('\n') }
        result.toString()
    }

    Path getTargetDir() {
        config.getStoreDir() ?: workDir
    }

    def getScratch() {
        config.scratch
    }

    String getWorkDirStr() {
        if( !workDir )
            return null
        return workDir.toUriString()
    }

    Path getWorkDirFor(HashCode hash) {
        FileHelper.getWorkFolder(processor.executor.getWorkDir(), hash)
    }

    static final public String CMD_LOG = '.command.log'
    static final public String CMD_SCRIPT = '.command.sh'
    static final public String CMD_INFILE = '.command.in'
    static final public String CMD_OUTFILE = '.command.out'
    static final public String CMD_ERRFILE = '.command.err'
    static final public String CMD_EXIT = '.exitcode'
    static final public String CMD_START = '.command.begin'
    static final public String CMD_RUN = '.command.run'
    static final public String CMD_STUB = '.command.stub'
    static final public String CMD_TRACE = '.command.trace'


    String toString( ) {
        "id: $id; name: $name; type: $type; exit: ${exitStatus==Integer.MAX_VALUE ? '-' : exitStatus}; error: $error; workDir: $workDir"
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

    @Memoized
    Path getCondaEnv() {
        if( !config.conda )
            return null

        final cfg = processor.session.config.conda as Map ?: Collections.emptyMap()
        final cache = new CondaCache(new CondaConfig(cfg))
        cache.getCachePathFor(config.conda as String)
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
        else if( !config.container ) {
            return null
        }
        else {
            imageName = config.container as String
        }

        final cfg = getContainerConfig()
        final handler = new ContainerHandler(cfg)
        handler.normalizeImageName(imageName)
    }

    /**
     * @return The {@link ContainerConfig} object associated to this task
     */
    ContainerConfig getContainerConfig() {
        processor.getSession().getContainerConfig()
    }


    /**
     * @return {@true} when the process must run within a container and the docker engine is enabled
     */
    boolean isDockerEnabled() {
        if( isContainerExecutable() )
            return true

        def config = getContainerConfig()
        return config && config.engine == 'docker' && config.enabled
    }

    /**
     * @return {@true} when the process runs an *executable* container
     */
    boolean isContainerExecutable() {
        getConfig().container == true
    }

    boolean isContainerNative() {
        processor.executor?.isContainerNative() ?: false
    }

    boolean isContainerEnabled() {
        getConfig().container && (isContainerExecutable() || getContainerConfig().enabled || isContainerNative())
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
            throw new ProcessUnrecoverableException("Process `$name` script contains error(s)", e)
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
            throw new ProcessTemplateException("Process `${processor.name}` can't find template file: $template")
        }
        catch( MissingPropertyException e ) {
            throw new ProcessTemplateException("No such property `$e.property` -- check template file: $template")
        }
        catch( Exception e ) {
            def msg = (e.message ?: "Unexpected error") + " -- check template file: $template"
            throw new ProcessTemplateException(msg, e)
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

