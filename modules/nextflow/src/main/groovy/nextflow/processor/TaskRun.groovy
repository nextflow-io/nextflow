/*
 * Copyright 2013-2024, Seqera Labs
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

import java.nio.file.FileSystems
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap
import java.util.function.Function

import com.google.common.hash.HashCode
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.conda.CondaCache
import nextflow.conda.CondaConfig
import nextflow.container.ContainerConfig
import nextflow.container.DockerConfig
import nextflow.container.resolver.ContainerInfo
import nextflow.container.resolver.ContainerMeta
import nextflow.container.resolver.ContainerResolver
import nextflow.container.resolver.ContainerResolverProvider
import nextflow.exception.ProcessException
import nextflow.exception.ProcessTemplateException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.script.BodyDef
import nextflow.script.ScriptType
import nextflow.script.TaskClosure
import nextflow.script.bundle.ResourcesBundle
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdInParam
import nextflow.script.params.ValueOutParam
import nextflow.packages.PackageManager
import nextflow.packages.PackageSpec
import nextflow.spack.SpackCache
/**
 * Models a task instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class TaskRun implements Cloneable {

    final private ConcurrentHashMap<String,?> cache0 = new ConcurrentHashMap()

    /**
     * Task unique id
     */
    TaskId id

    /**
     * Task index within its execution group
     */
    Integer index

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
    void echoStdout(Session session) {

        // print the stdout
        if( stdout instanceof Path ) {
            try {
                session.printConsole(stdout)
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
            session.printConsole(stdout.toString())
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
            return Collections.<String>emptyList()
        }
    }

    List<String> dumpStderr(int n = 50) {

        try {
            return dumpObject(stderr,n)
        }
        catch( Exception e ) {
            log.debug "Unable to dump error of process '$name' -- Cause: ${e}"
            return Collections.<String>emptyList()
        }
    }

    List<String> dumpLogFile(int n = 50) {
        if( !workDir || workDir.fileSystem!=FileSystems.default )
            return Collections.<String>emptyList()
        try {
            return dumpObject(workDir.resolve(CMD_LOG),n)
        }
        catch( Exception e ) {
            log.debug "Unable to dump error of process '$name' -- Cause: ${e}"
            return Collections.<String>emptyList()
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

        return result ?: Collections.<String>emptyList()
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
     * The number of times the submit of the task has been retried
     */
    volatile int submitRetries

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

    /**
     * Unique key for the container used by this task
     */
    volatile String containerKey

    TaskConfig config

    TaskContext context

    Set<String> templateVars

    TaskProcessor.RunType runType = TaskProcessor.RunType.SUBMIT

    BodyDef body

    TaskRun clone() {
        final taskClone = (TaskRun)super.clone()
        taskClone.context = context.clone()
        taskClone.config = config.clone()
        taskClone.config.setContext(taskClone.context)
        taskClone.cache0.clear()
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

    String lazyName() {
        if( name )
            return name
        // fallback on the current task index, however do not set the 'name' attribute
        // so it has a chance to recover the 'sampleId' at next invocation
        return processor.singleton ? processor.name : "$processor.name ($index)"
    }

    String getName() {
        if( name )
            return name

        final baseName = processor.name
        if( config.containsKey('tag') )
            try {
                // -- look-up the 'sampleId' property, and if everything is fine
                //    cache this value in the 'name' attribute
                return name = "$baseName (${String.valueOf(config.tag).trim()})"
            }
            catch( IllegalStateException e ) {
                log.debug "Cannot access `tag` property for task: $baseName ($index)"
            }
            catch( Exception e ) {
                log.debug "Unable to evaluate `tag` property for task: $baseName ($index)", e
            }

        return lazyName()
    }

    String getScript() {
        if( script instanceof Path ) {
            return script.text
        }
        else {
            return script?.toString()
        }
    }

    String getTraceScript() {
        return template!=null && body?.source
            ? body.source
            : getScript()
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

        final allFiles = getInputFiles().values()
        final result = new HashMap<String,Path>(allFiles.size())
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
     */
    List<String> getOutputFilesNames() {
        // note: use an explicit function instead of a closure or lambda syntax, otherwise
        // when calling this method from a subclass it will result into a MissingMethodExeception
        // see  https://issues.apache.org/jira/browse/GROOVY-2433
        cache0.computeIfAbsent('outputFileNames', new Function<String,List<String>>() {
            @Override
            List<String> apply(String s) {
                return getOutputFilesNames0()
            }})
    }

    private List<String> getOutputFilesNames0() {
        def result = []

        for( FileOutParam param : getOutputsByType(FileOutParam).keySet() ) {
            result.addAll( param.getFilePatterns(context, workDir) )
        }

        return result.unique()
    }

    /**
     * Get the map of *input* objects by the given {@code InParam} type
     *
     * @param types One or more subclass of {@code InParam}
     * @return An associative array containing all the objects for the specified type
     */
    def <T extends InParam> Map<T,Object> getInputsByType( Class<T>... types ) {

        def result = [:]
        for( def it : inputs ) {
            if( types.contains(it.key.class) )
                result << it
        }
        return result
    }

    /**
     * Get the map of *output* objects by the given {@code InParam} type
     *
     * @param types One or more subclass of {@code InParam}
     * @return An associative array containing all the objects for the specified type
     */
    def <T extends OutParam> Map<T,Object> getOutputsByType( Class<T>... types ) {
        def result = [:]
        for( def it : outputs ) {
            if( types.contains(it.key.class) )
                result << it
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
    static final public String CMD_STAGE = '.command.stage'
    static final public String CMD_TRACE = '.command.trace'
    static final public String CMD_ENV = '.command.env'


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

    List<String> getOutputEnvNames() {
        final items = getOutputsByType(EnvOutParam)
        if( !items )
            return List.<String>of()
        final result = new ArrayList<String>(items.size())
        for( EnvOutParam it : items.keySet() ) {
            if( !it.name ) throw new IllegalStateException("Missing output environment name - offending parameter: $it")
            result.add(it.name)
        }
        return result
    }

    /**
     * @return A {@link Map} instance holding a collection of key-pairs
     * where the key represents a environment variable name holding the command
     * output and the value the command the executed.
     */
    Map<String,String> getOutputEvals() {
        final items = getOutputsByType(CmdEvalParam)
        final result = new LinkedHashMap(items.size())
        for( CmdEvalParam it : items.keySet() ) {
            if( !it.name ) throw new IllegalStateException("Missing output eval name - offending parameter: $it")
            result.put(it.name, it.getTarget(context))
        }
        return result
    }

    Path getCondaEnv() {
        // note: use an explicit function instead of a closure or lambda syntax, otherwise
        // when calling this method from a subclass it will result into a MissingMethodExeception
        // see  https://issues.apache.org/jira/browse/GROOVY-2433
        cache0.computeIfAbsent('condaEnv', new Function<String,Path>() {
            @Override
            Path apply(String it) {
                return getCondaEnv0()
            }})
    }

    private Path getCondaEnv0() {
        if( !config.conda || !getCondaConfig().isEnabled() )
            return null

        // Show deprecation warning if new package system is enabled
        if (PackageManager.isEnabled(processor.session)) {
            log.warn "The 'conda' directive is deprecated when preview.package is enabled. Use 'package \"${config.conda}\", provider: \"conda\"' instead"
        }

        final cache = new CondaCache(getCondaConfig())
        cache.getCachePathFor(config.conda as String)
    }


    Path getSpackEnv() {
        // note: use an explicit function instead of a closure or lambda syntax, otherwise
        // when calling this method from a subclass it will result into a MissingMethodExeception
        // see  https://issues.apache.org/jira/browse/GROOVY-2433
        cache0.computeIfAbsent('spackEnv', new Function<String,Path>() {
            @Override
            Path apply(String it) {
                return getSpackEnv0()
            }})
    }

    private Path getSpackEnv0() {
        if( !config.spack || !processor.session.getSpackConfig().isEnabled() )
            return null

        final String arch = config.getArchitecture()?.spackArch

        final cache = new SpackCache(processor.session.getSpackConfig())
        cache.getCachePathFor(config.spack as String, arch)
    }

    PackageSpec getPackageSpec() {
        // note: use an explicit function instead of a closure or lambda syntax
        cache0.computeIfAbsent('packageSpec', new Function<String,PackageSpec>() {
            @Override
            PackageSpec apply(String it) {
                return getPackageSpec0()
            }})
    }

    private PackageSpec getPackageSpec0() {
        if (!PackageManager.isEnabled(processor.session))
            return null

        if (!config.package)
            return null

        def packageManager = new PackageManager(processor.session)
        
        // Parse the package configuration
        def packageDef = config.package
        def defaultProvider = processor.session.config.navigate('packages.provider', 'conda') as String
        
        try {
            return PackageManager.parseSpec(packageDef, defaultProvider)
        } catch (Exception e) {
            log.warn "Failed to parse package specification: ${e.message}"
            return null
        }
    }

    protected ContainerInfo containerInfo() {
        // note: use an explicit function instead of a closure or lambda syntax, otherwise
        // when calling this method from a subclass it will result into a MissingMethodException
        // see  https://issues.apache.org/jira/browse/GROOVY-2433
        cache0.computeIfAbsent('containerInfo', new Function<String,ContainerInfo>() {
            @Override
            ContainerInfo apply(String s) {
                return containerInfo0()
            }})
    }

    @Memoized
    protected ContainerResolver containerResolver() {
        ContainerResolverProvider.load()
    }

    private ContainerInfo containerInfo0() {
        // fetch the container image from the config
        def configImage = config.getContainer()
        // the boolean `false` literal can be provided
        // to signal the absence of the container
        if( configImage == false )
            return null
        if( !configImage )
            configImage = null

        final info = containerResolver().resolveImage(this, configImage as String)
        // track the key of the container used
        if( info!=null )
            this.containerKey = info.hashKey
        // return the info
        return info
    }

    /**
     * The name of a docker container where the task is supposed to run when provided
     */
    String getContainer() {
        final info = containerInfo()
        return info?.target
    }

    String getContainerFingerprint() {
        final info = containerInfo()
        return info?.hashKey
    }

    boolean isContainerReady() {
       return containerKey
            ? containerResolver().isContainerReady(containerKey)
            : true
    }

    ContainerMeta containerMeta() {
        return containerKey
            ? containerResolver().getContainerMeta(containerKey)
            : null
    }

    String getContainerPlatform() {
        final result = config.getArchitecture()
        return result ? result.dockerArch : containerResolver().defaultContainerPlatform()
    }

    ResourcesBundle getModuleBundle() {
        return this.getProcessor().getModuleBundle()
    }

    /**
     * @return The {@link ContainerConfig} object associated to this task
     */
    ContainerConfig getContainerConfig() {
        // get the container engine expected to be used by this executor
        final sess = this.getProcessor().getSession()
        final exe = this.getProcessor().getExecutor()
        final eng = exe.containerConfigEngine()
        // when 'eng' is null the setting for the current engine marked as 'enabled' will be used
        final result
                = sess.getContainerConfig(eng)
                ?: new DockerConfig([:])
        // if a configuration is found is expected to enabled by default
        if( exe.isContainerNative() ) {
            result.setEnabled(true)
        }
        return result
    }

    /**
     * @return {@true} when the process must run within a container and the docker engine is enabled
     */
    boolean isDockerEnabled() {
        final config = getContainerConfig()
        return config && config.engine == 'docker' && config.enabled
    }

    boolean isContainerNative() {
        return processor.executor?.isContainerNative() ?: false
    }

    boolean isArray() {
        return false
    }

    boolean isContainerEnabled() {
        return getContainerConfig().isEnabled() && getContainer()!=null
    }

    boolean isSecretNative() {
        return processor.executor?.isSecretNative() ?: false
    }

    boolean isSuccess( status = exitStatus ) {
        if( status == null )
            return false

        if( status instanceof String )
            status = status.toInteger()

        if( status == Integer.MAX_VALUE )
            return false

        return status == TaskConfig.EXIT_ZERO
    }

    /**
     * Given the task body token declared in the process definition:
     * 1) assign to the task run instance the code closure to execute
     * 2) extract the process code `source`
     * 3) assign the `script` code to execute
     *
     * @param body A {@code BodyDef} object instance
     */
    void resolve(BodyDef body)  {
        processor.session.stubRun && config.getStubBlock()
            ? resolveStub(config.getStubBlock())
            : resolveBody(body)
    }

    protected void resolveBody(BodyDef body) {

        // -- initialize the task code to be executed
        this.code = body.closure.clone() as Closure
        this.code.delegate = this.context
        this.code.setResolveStrategy(Closure.DELEGATE_ONLY)

        // -- set the task source
        this.body = body
        // note: this may be overwritten when a template file is used
        this.source = body.source

        if( body.type != ScriptType.SCRIPTLET )
            return

        // Important!
        // when the task is implemented by a script string
        // Invoke the closure which returns the script with all the variables replaced with the actual values
        try {
            final result = code.call()
            if ( result instanceof Path ) {
                script = renderTemplate(result, body.isShell)
            }
            else if( result != null && body.isShell ) {
                script = renderScript(result)
            }
            else {
                script = result.toString()
            }
        }
        catch( ProcessException e ) {
            throw e
        }
        catch( Throwable e ) {
            throw new ProcessUnrecoverableException("Process `${getName()}` script contains error(s)", e)
        }

        // make sure the task script is not empty
        if( !script )
            throw new ProcessUnrecoverableException("Process `${getName()}` script is empty")
    }

    protected void resolveStub(TaskClosure block) {
        this.code = block.clone() as Closure
        this.code.delegate = this.context
        this.code.setResolveStrategy(Closure.DELEGATE_ONLY)

        // -- set the task source
        // note: this may be overwritten when a template file is used
        this.source = block.getSource()

        try {
            final result = code.call()
            if ( result instanceof Path ) {
                script = renderTemplate(result)
            }
            else {
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
            final engine = new TaskTemplateEngine(processor.@grengine)
            if( shell ) {
                engine.setPlaceholder(placeholderChar())
            }
            // eval the template string
            engine.eval(source, context)
            templateVars = engine.getVariableNames()
            return engine.result
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

        final engine = new TaskTemplateEngine(processor.@grengine)
                .setPlaceholder(placeholderChar())
                .setEnableShortNotation(false)
                .eval(script.toString(), context)
        // fetch the template vars
        templateVars = engine.getVariableNames()
        // finally return the evaluated string
        return engine.result
    }

    protected char placeholderChar() {
        (config.placeholder ?: '!') as char
    }

    protected Set<String> getVariableNames() {
        if( templateVars!=null )
            return templateVars - processor.getDeclaredNames()
        else
            return context.getVariableNames()
    }

    /**
     * @param variableNames The collection of variables referenced in the task script
     * @param binding The script global binding
     * @param context The task variable context
     * @return The set of task variables accessed in global script context and not declared as input/output
     */
    Map<String,Object> getGlobalVars(Binding binding) {
        final variableNames = getVariableNames()
        final result = new HashMap(variableNames.size())
        final processName = name

        def itr = variableNames.iterator()
        while( itr.hasNext() ) {
            final varName = itr.next()

            final p = varName.indexOf('.')
            final baseName = p !=- 1 ? varName.substring(0,p) : varName

            if( context.isLocalVar(baseName) ) {
                // when the variable belong to the task local context just ignore it,
                // because it must has been provided as an input parameter
                continue
            }

            if( binding.hasVariable(baseName) ) {
                def value
                try {
                    if( p == -1 ) {
                        value = binding.getVariable(varName)
                    }
                    else {
                        value = processor.grengine.run(varName, binding)
                    }
                }
                catch( MissingPropertyException | NullPointerException e ) {
                    value = null
                    log.trace "Process `${processName}` cannot access global variable `$varName` -- Cause: ${e.message}"
                }

                // value for 'workDir' and 'baseDir' folders are added always as string
                // in order to avoid to invalid the cache key when resuming the execution
                if( varName=='workDir' || varName=='baseDir' || varName=='projectDir' )
                    value = value.toString()

                result.put( varName, value )
            }

        }
        return result
    }

    TaskBean toTaskBean() {
        return new TaskBean(this)
    }

    CondaConfig getCondaConfig() {
        return processor.session.getCondaConfig()
    }


    String getStubSource() {
        return config?.getStubBlock()?.source
    }
}

