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

package nextflow.executor

import static java.nio.file.StandardOpenOption.*

import java.nio.file.FileSystemException
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.container.ContainerBuilder
import nextflow.container.DockerBuilder
import nextflow.container.SingularityBuilder
import nextflow.exception.ProcessException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.secret.SecretsLoader
import nextflow.util.Escape
import nextflow.util.MemoryUnit
/**
 * Builder to create the Bash script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BashWrapperBuilder {

    private static MemoryUnit DEFAULT_STAGE_FILE_THRESHOLD = MemoryUnit.of('1 MB')
    private static int DEFAULT_WRITE_BACK_OFF_BASE = 3
    private static int DEFAULT_WRITE_BACK_OFF_DELAY = 250
    private static int DEFAULT_WRITE_MAX_ATTEMPTS = 5

    private MemoryUnit stageFileThreshold = SysEnv.get('NXF_WRAPPER_STAGE_FILE_THRESHOLD') as MemoryUnit ?: DEFAULT_STAGE_FILE_THRESHOLD
    private int writeBackOffBase = SysEnv.get('NXF_WRAPPER_BACK_OFF_BASE') as Integer ?: DEFAULT_WRITE_BACK_OFF_BASE
    private int writeBackOffDelay = SysEnv.get('NXF_WRAPPER_BACK_OFF_DELAY') as Integer ?: DEFAULT_WRITE_BACK_OFF_DELAY
    private int writeMaxAttempts = SysEnv.get('NXF_WRAPPER_MAX_ATTEMPTS') as Integer ?: DEFAULT_WRITE_MAX_ATTEMPTS

    static final public KILL_CMD = '[[ "$pid" ]] && nxf_kill $pid'

    static final private ENDL = '\n'

    static final public List<String> BASH

    static private int level

    @PackageScope
    static public String systemOsName = System.getProperty('os.name')

    static {
        /*
         * Env variable `NXF_DEBUG` is used to control debug options in executed BASH scripts
         * - 0: no debug
         * - 1: dump current environment in the `.command.log` file and trace the execution of user script
         * - 2: trace the execution of wrapper scripts
         */
        def str = System.getenv('NXF_DEBUG')
        try {
            level = str ? str as int : 0
        }
        catch( Exception e ) {
            log.warn "Invalid value for `NXF_DEBUG` variable: $str -- See http://www.nextflow.io/docs/latest/config.html#environment-variables"
        }
        BASH = Collections.unmodifiableList( level > 0 ? ['/bin/bash','-uex'] : ['/bin/bash','-ue'] )

    }


    @Delegate
    ScriptFileCopyStrategy copyStrategy

    @Delegate
    private TaskBean bean

    private boolean runWithContainer

    private ContainerBuilder containerBuilder

    private Path scriptFile

    private Path inputFile

    private Path startedFile

    private Path exitedFile

    private Path wrapperFile

    private Path stageFile

    private String stageScript

    private BashTemplateEngine engine = new BashTemplateEngine()

    BashWrapperBuilder( TaskRun task ) {
        this(new TaskBean(task))
    }

    BashWrapperBuilder( TaskBean bean, ScriptFileCopyStrategy strategy = null ) {
        this.bean = bean
        this.copyStrategy = strategy ?: new SimpleFileCopyStrategy(bean)
    }

    /** only for testing -- do not use */
    protected BashWrapperBuilder() { }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String getScratchDirectoryCommand() {

        // convert to string for safety
        final scratchStr = scratch?.toString()

        if( scratchStr == null || scratchStr == 'false' ) {
            return null
        }

        /*
         * when 'scratch' is defined as a bool value
         * try to use the 'TMP' variable, if does not exist fallback to a tmp folder
         */
        if( scratchStr == 'true' ) {
            return 'NXF_SCRATCH="$(set +u; nxf_mktemp $TMPDIR)"'
        }

        if( scratchStr.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NXF_SCRATCH="$(nxf_mktemp /dev/shm)"'
        }

        return "NXF_SCRATCH=\"\$(set +u; nxf_mktemp $scratchStr)\""
    }

    protected boolean shouldUnstageOutputs() {
        return targetDir && workDir!=targetDir
    }

    protected boolean fixOwnership() {
        systemOsName == 'Linux' && containerConfig?.fixOwnership && runWithContainer && containerConfig.engine == 'docker' // <-- note: only for docker (other container runtimes are not affected)
    }

    protected isMacOS() {
        systemOsName.startsWith('Mac')
    }

    @PackageScope String buildNew0() {
        final template = BashWrapperBuilder.class.getResourceAsStream('command-run.txt')
        try {
            return buildNew0(template.newReader())
        }
        finally {
            template.close()
        }
    }

    /**
     * Generate a Bash script to be appended to the task command script
     * that takes care of capturing the process output environment variables
     * and evaluation commands
     *
     * @param outEnvs
     *      The list of environment variables names whose value need to be captured
     * @param outEvals
     *      The set of commands to be evaluated to determine the output value to be captured
     * @return
     *      The Bash script to capture the output environment and eval commands
     */
    protected String getOutputEnvCaptureSnippet(List<String> outEnvs, Map<String,String> outEvals) {
        // load the env template
        final template = BashWrapperBuilder.class
            .getResourceAsStream('command-env.txt')
            .newReader()
        final binding = Map.of('env_file', TaskRun.CMD_ENV)
        final result = new StringBuilder()
        result.append( engine.render(template, binding) )
        appendOutEnv(result, outEnvs)
        appendOutEval(result, outEvals)
        return result.toString()
    }

    /**
     * Render a Bash script to capture the one or more env variables
     *
     * @param result A {@link StringBuilder} instance to which append the result Bash script
     * @param outEnvs The environment variables to be captured
     */
    protected void appendOutEnv(StringBuilder result, List<String> outEnvs) {
        if( outEnvs==null )
            outEnvs = List.<String>of()
        // out env
        for( String key : outEnvs ) {
            result << "#\n"
            result << "echo $key=\"\${$key[@]}\" >> ${TaskRun.CMD_ENV}\n"
            result << "echo /$key/ >> ${TaskRun.CMD_ENV}\n"
        }
    }

    /**
     * Render a Bash script to capture the result of one or more commands
     * evaluated in the task script context
     *
     * @param result
     *      A {@link StringBuilder} instance to which append the result Bash script
     * @param outEvals
     *      A {@link Map} of key-value pairs modeling the commands to be evaluated;
     *      where the key represents the environment variable (name) holding the
     *      resulting output, and the pair value represent the Bash command to be
     *      evaluated.
     */
    protected void appendOutEval(StringBuilder result, Map<String,String> outEvals) {
        if( outEvals==null )
            outEvals = Map.<String,String>of()
        // out eval
        for( Map.Entry<String,String> eval : outEvals ) {
            result << "#\n"
            result <<"nxf_eval_cmd STDOUT STDERR bash -c \"${eval.value.replace('"','\\\"')}\"\n"
            result << 'status=$?\n'
            result << 'if [ $status -eq 0 ]; then\n'
            result << "  echo $eval.key=\"\$STDOUT\" >> ${TaskRun.CMD_ENV}\n"
            result << "  echo /$eval.key/=exit:0 >> ${TaskRun.CMD_ENV}\n"
            result << 'else\n'
            result << "  echo $eval.key=\"\$STDERR\" >> ${TaskRun.CMD_ENV}\n"
            result << "  echo /$eval.key/=exit:\$status >> ${TaskRun.CMD_ENV}\n"
            result << 'fi\n'
        }
    }

    protected String stageCommand(String stagingScript) {
        if( !stagingScript )
            return null

        final header = "# stage input files\n"
        // enable only when the stage uses the default file system, i.e. it's not a remote object storage file
        // see https://github.com/nextflow-io/nextflow/issues/4279
        if( stageFile.fileSystem == FileSystems.default && stagingScript.size() >= stageFileThreshold.bytes ) {
            stageScript = stagingScript
            return header + "bash ${stageFile}"
        }
        else
            return header + stagingScript
    }

    protected Map<String,String> makeBinding() {
        /*
         * initialise command files
         */
        scriptFile = workDir.resolve(TaskRun.CMD_SCRIPT)
        inputFile = workDir.resolve(TaskRun.CMD_INFILE)
        startedFile = workDir.resolve(TaskRun.CMD_START)
        exitedFile = workDir.resolve(TaskRun.CMD_EXIT)
        wrapperFile = workDir.resolve(TaskRun.CMD_RUN)
        stageFile = workDir.resolve(TaskRun.CMD_STAGE)

        // set true when running with through a container engine
        runWithContainer = containerEnabled && !containerNative

        // whenever it has to change to the scratch directory
        final changeDir = getScratchDirectoryCommand()

        /*
         * create the container launcher command if needed
         */
        containerBuilder = runWithContainer ? createContainerBuilder(changeDir) : null

        // ugly, it should be done somewhere else
        if( script )
            script = TaskProcessor.normalizeScript(script, shell)

        /*
         * fetch the script interpreter i.e. BASH, Perl, Python, etc
         */
        final interpreter = TaskProcessor.fetchInterpreter(script)

        /*
         * append to the command script a prolog to capture the declared
         * output environment (variable) and evaluation commands
         */
        if( outputEnvNames || outputEvals ) {
            if( !isBash(interpreter) && outputEnvNames )
                throw new IllegalArgumentException("Process output of type 'env' is only allowed with Bash process scripts -- Current interpreter: $interpreter")
            if( !isBash(interpreter) && outputEvals )
                throw new IllegalArgumentException("Process output of type 'eval' is only allowed with Bash process scripts -- Current interpreter: $interpreter")
            script += getOutputEnvCaptureSnippet(outputEnvNames, outputEvals)
        }

        final binding = new HashMap<String,String>(20)
        binding.header_script = headerScript
        binding.task_metadata = getTaskMetadata()
        binding.task_name = name
        binding.helpers_script = getHelpersScript()

        if( runWithContainer ) {
            binding.container_boxid = 'export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d \'\\r\\n\')"'
            binding.container_helpers = containerBuilder.getScriptHelpers()
            binding.kill_cmd = containerBuilder.getKillCommand()
        }
        else {
            binding.container_boxid = null
            binding.container_helpers = null
            binding.kill_cmd = KILL_CMD
        }

        binding.cleanup_cmd = getCleanupCmd(changeDir)
        binding.sync_cmd = getSyncCmd()
        binding.scratch_cmd = ( changeDir ?: "NXF_SCRATCH=''" )

        binding.exit_file = exitFile(exitedFile)
        binding.touch_file = touchFile(startedFile)

        binding.module_load = getModuleLoadSnippet()
        binding.before_script = getBeforeScriptSnippet()
        binding.conda_activate = getCondaActivateSnippet()
        binding.spack_activate = getSpackActivateSnippet()

        /*
         * add the task environment
         */
        final env = copyStrategy.getEnvScript(environment, runWithContainer)
        if( runWithContainer ) {
            binding.task_env = null
            binding.container_env = env
        }
        else {
            binding.task_env = env
            binding.container_env = null
        }

        /*
         * add the task secrets
         */
        if( !isSecretNative() ) {
            binding.secrets_env = getSecretsEnv()
        }
        else {
            binding.secrets_env = null
        }

        /*
         * staging input files when required
         */
        final stagingScript = copyStrategy.getStageInputFilesScript(inputFiles)
        binding.stage_inputs = stageCommand(stagingScript)

        binding.stdout_file = TaskRun.CMD_OUTFILE
        binding.stderr_file = TaskRun.CMD_ERRFILE
        binding.trace_file = TaskRun.CMD_TRACE

        binding.trace_cmd = getTraceCommand(interpreter)
        binding.launch_cmd = getLaunchCommand(interpreter,env)
        binding.stage_cmd = getStageCommand()
        binding.unstage_cmd = getUnstageCommand()
        binding.unstage_controls = changeDir || shouldUnstageOutputs() ? getUnstageControls() : null

        if( changeDir || shouldUnstageOutputs() ) {
            binding.unstage_outputs = copyStrategy.getUnstageOutputFilesScript(outputFiles,targetDir)
        }
        else {
            binding.unstage_outputs = null
        }

        binding.after_script = afterScript ? "# 'afterScript' directive\n$afterScript" : null

        // patch root ownership problem on files created with docker
        binding.fix_ownership = fixOwnership() ? "[ \${NXF_OWNER:=''} ] && (shopt -s extglob; GLOBIGNORE='..'; chown -fR --from root \$NXF_OWNER ${workDir}/{*,.*}) || true" : null

        binding.trace_script = isTraceRequired() ? getTraceScript(binding) : null
        
        return binding
    }

    protected String getSecretsEnv() {
        return SecretsLoader.isEnabled()
                ? SecretsLoader.instance.load() .getSecretsEnv(secretNames)
                : null
    }

    protected boolean isBash(String interpreter) {
        interpreter.tokenize(' /').contains('bash')
    }

    protected String getTraceScript(Map binding) {
        def res = BashWrapperBuilder.class.getResource('command-trace.txt')
        engine.render(res.newReader(), binding)
    }

    @PackageScope String buildNew0(BufferedReader template) {
        final binding = makeBinding()
        engine.render(template, binding)
    }

    /**
     * Build up the BASH wrapper script file which will launch the user provided script
     * @return The {@code Path} of the created wrapper script
     */
    Path build() {
        assert workDir, "Missing 'workDir' property in BashWrapperBuilder object"
        assert script, "Missing 'script' property in BashWrapperBuilder object"

        if( statsEnabled && isMacOS() && !isContainerEnabled() )
            log.warn1("Task runtime metrics are not reported when using macOS without a container engine")

        final wrapper = buildNew0()
        final result = write0(targetWrapperFile(), wrapper)
        write0(targetScriptFile(), script)
        if( input != null )
            write0(targetInputFile(), input.toString())
        if( stageScript != null )
            write0(targetStageFile(), stageScript)
        return result
    }

    protected Path targetWrapperFile() { return wrapperFile }

    protected Path targetScriptFile() { return scriptFile }

    protected Path targetInputFile() { return inputFile }

    protected Path targetStageFile() { return stageFile }

    private Path write0(Path path, String data) {
        int attempt=0
        while( true ) {
            try {
                try (BufferedWriter writer=Files.newBufferedWriter(path, CREATE,WRITE,TRUNCATE_EXISTING)) {
                    writer.write(data)
                }
                return path
            }
            catch (Exception e) {
                if( !isRetryable0(e) )
                    throw e
                final isLocalFS = path.getFileSystem()==FileSystems.default
                // the retry logic is needed for non-local file system such as S3.
                // when the file is local fail without retrying
                if( isLocalFS || ++attempt>=writeMaxAttempts )
                    throw new ProcessException("Unable to create file ${path.toUriString()}", e)
                // use an exponential delay before making another attempt
                final delay = (Math.pow(writeBackOffBase, attempt) as long) * writeBackOffDelay
                log.debug "Unexpected error writing '${path.toUriString()}'; attempt: $attempt - cause: ${e.message}"
                Thread.sleep(delay)
            }
        }
    }

    static protected boolean isRetryable0(Exception e) {
        if( e instanceof FileSystemException )
            return true
        if( e instanceof SocketException )
            return true
        if( e instanceof RuntimeException )
            return true
        if( e.class.getSimpleName() == 'HttpResponseException' )
            return true
        return false
    }

    protected String getTaskMetadata() {
        final lines = new StringBuilder()
        lines << '### ---\n'
        lines << "### name: '${bean.name}'\n"
        if( bean.arrayIndexName ) {
            lines << '### array:\n'
            lines << "###   index-name: ${bean.arrayIndexName}\n"
            lines << "###   index-start: ${bean.arrayIndexStart}\n"
            lines << "###   work-dirs:\n"
            for( Path it : bean.arrayWorkDirs )
                lines << "###   - ${Escape.path(FilesEx.toUriString(it))}\n"
        }

        if( containerConfig?.isEnabled() )
            lines << "### container: '${bean.containerImage}'\n"

        if( outputFiles.size() > 0 ) {
            lines << '### outputs:\n'
            for( final output : bean.outputFiles )
                lines << "### - '${output}'\n"
        }

        lines << '### ...\n'
    }

    protected String getHelpersScript() {
        def result = new StringBuilder()

        def s1 = copyStrategy.getBeforeStartScript()
        if( s1 )
            result.append(s1).append('\n')

        def s2 = moduleNames ? engine.render(BashWrapperBuilder.class.getResource('modules-env.txt').newReader(), Collections.emptyMap()) : null
        if( s2 )
            result.append(s2).append('\n')

        result.size() ? result.toString() : null
    }

    private String getBeforeScriptSnippet() {
        beforeScript ? "# beforeScript directive\n$beforeScript\n" : null
    }

    private String getModuleLoadSnippet() {
        if( !moduleNames )
            return null
        String result=''
        for( String it : moduleNames) {
            result += moduleLoad(it) + ENDL
        }
        return result
    }

    private String getCondaActivateSnippet() {
        if( !condaEnv )
            return null
        final command = useMicromamba
            ? 'eval "$(micromamba shell hook --shell bash)" && micromamba activate'
            : 'source $(conda info --json | awk \'/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }\')/bin/activate'
        return """\
            # conda environment
            ${command} ${Escape.path(condaEnv)}
            """.stripIndent()
    }

    private String getSpackActivateSnippet() {
        if( !spackEnv )
            return null
        def result = "# spack environment\n"
        result += 'spack env activate -d '
        result += "${Escape.path(spackEnv)}\n"
        return result
    }

    protected String getTraceCommand(String interpreter) {
        String result = "${interpreter} ${fileStr(scriptFile)}"
        if( input != null )
            result += pipeInputFile(inputFile)

        return result
    }

    protected boolean isTraceRequired() {
        statsEnabled || fixOwnership()
    }

    protected String shellPath() {
        // keep the shell path as "/bin/bash" when a non-custom "shell" attribute is specified
        // to not introduce unexpected changes due to the fact BASH is defined as "/bin/bash -eu" by default
        return shell.is(BASH)
            ? "/bin/bash"
            : shell.join(' ')
    }

    protected String getLaunchCommand(String interpreter, String env) {
        /*
        * process stats
        */
        String launcher

        // NOTE: the isTraceRequired() check must match the logic in launchers (i.e. AwsBatchScriptLauncher)
        // that determines when to stage the file.
        final traceWrapper = isTraceRequired()
        if( traceWrapper ) {
            // executes the stub which in turn executes the target command
            launcher = "${shellPath()} ${fileStr(wrapperFile)} nxf_trace"
        }
        else {
            launcher = "${interpreter} ${fileStr(scriptFile)}"
        }

        /*
         * create the container engine command when needed
         */
        if( containerBuilder ) {
            String cmd = env ? 'eval $(nxf_container_env); ' + launcher : launcher
            // wrap the command with an extra bash invocation either :
            // - to propagate the container environment or
            // - to change in the task work directory as required by singularity
            final needChangeTaskWorkDir = containerBuilder instanceof SingularityBuilder
            if( (env || needChangeTaskWorkDir) && !containerConfig.entrypointOverride() ) {
                if( needChangeTaskWorkDir )
                    cmd = 'cd $NXF_TASK_WORKDIR; ' + cmd
                cmd = "/bin/bash -c \"$cmd\""
            }
            launcher = containerBuilder.getRunCommand(cmd)
        }

        /*
         * pipe the input file on the command standard input
         */
        if( !traceWrapper && input != null ) {
            launcher += pipeInputFile(inputFile)
        }

        return launcher
    }


    private String copyFileToWorkDir(String fileName) {
        copyFile(fileName, workDir.resolve(fileName))
    }
    

    String getCleanupCmd(String scratch) {
        String result = ''
        // -- cleanup the scratch dir
        if( scratch && cleanup != false ) {
            result += (containerBuilder !instanceof DockerBuilder ? 'rm -rf $NXF_SCRATCH || true' : '(sudo -n true && sudo rm -rf "$NXF_SCRATCH" || rm -rf "$NXF_SCRATCH")&>/dev/null || true')
            result += '\n'
        }
        // -- remove the container in this way because 'docker run --rm'  fail in some cases -- see https://groups.google.com/d/msg/docker-user/0Ayim0wv2Ls/-mZ-ymGwg8EJ
        final remove = containerBuilder?.getRemoveCommand()
        if( remove ) {
            result += "${remove} &>/dev/null || true"
            result += '\n'
        }
        return result
    }

    String getExitScriptLegacy(String scratch) {
        def result = getCleanupCmd(scratch)
        result += 'exit $exit_status'
        result.readLines().join('\n  ')
    }

    String getSyncCmd() {
        if ( SysEnv.get( 'NXF_ENABLE_FS_SYNC' ) == "true" ) {
            return 'sync || true'
        }
        return null
    }

    @PackageScope
    ContainerBuilder createContainerBuilder0(String engine) {
        ContainerBuilder.create(engine, containerImage)
    }

    protected boolean getAllowContainerMounts() {
        return true
    }

    /**
     * Build a {@link DockerBuilder} object to handle Docker commands
     *
     * @param envFile A file containing environment configuration
     * @param changeDir String command to change to the working directory
     * @return A {@link DockerBuilder} instance
     */
    @PackageScope
    ContainerBuilder createContainerBuilder(String changeDir) {

        final engine = containerConfig.getEngine()
        ContainerBuilder builder = createContainerBuilder0(engine)

        /*
         * initialise the builder
         */
        // do not mount inputs when they are copied in the task work dir -- see #1105
        if( stageInMode != 'copy' && allowContainerMounts )
            builder.addMountForInputs(inputFiles)

        if( allowContainerMounts )
            builder.addMounts(binDirs)

        if(this.containerMount)
            builder.addMount(containerMount)

        // task work dir
        if( allowContainerMounts )
            builder.setWorkDir(workDir)

        // set the name
        builder.setName('$NXF_BOXID')

        if( this.containerMemory )
            builder.setMemory(containerMemory)

        if( this.containerCpus )
            builder.setCpus(containerCpus)

        if( this.containerCpuset )
            builder.addRunOptions(containerCpuset)

        // export task work directory
        builder.addEnv('NXF_TASK_WORKDIR')
        // export the nextflow script debug variable
        if( isTraceRequired() )
            builder.addEnv( 'NXF_DEBUG=${NXF_DEBUG:=0}')

        // add the user owner variable in order to patch root owned files problem
        if( fixOwnership() )
            builder.addEnv( 'NXF_OWNER=$(id -u):$(id -g)' )

        for( String var : containerConfig.getEnvWhitelist() ) {
            builder.addEnv(var)
        }

        // when secret are not managed by the execution platform natively
        // the secret names are added to the container env var white list
        if( !isSecretNative() && secretNames )  {
            for( String var : secretNames )
                builder.addEnv(var)
        }

        // set up run docker params
        builder.params(containerConfig)

        // extra rule for the 'auto' temp dir temp dir
        def temp = containerConfig.temp?.toString()
        if( temp == 'auto' || temp == 'true' ) {
            builder.setTemp( changeDir ? '$NXF_SCRATCH' : '$(nxf_mktemp)' )
        }

        if( containerConfig.containsKey('kill') )
            builder.params(kill: containerConfig.kill)

        if( containerConfig.writableInputMounts==false )
            builder.params(readOnlyInputs: true)

        if( this.containerConfig.entrypointOverride() )
            builder.params(entry: '/bin/bash')

        // give a chance to override any option with process specific `containerOptions`
        if( containerOptions ) {
            builder.addRunOptions(containerOptions)
        }

        // The current work directory should be mounted only when
        // the task is executed in a temporary scratch directory (ie changeDir != null)
        // See https://github.com/nextflow-io/nextflow/issues/1710
        builder.addMountWorkDir( changeDir as boolean || FileHelper.getWorkDirIsSymlink() )

        builder.build()
        return builder
    }

    String moduleLoad(String name) {
        int p = name.lastIndexOf('/')
        p != -1 ? "nxf_module_load ${name.substring(0,p)} ${name.substring(p+1)}" : "nxf_module_load ${name}"
    }

    protected String getStageCommand() { 'nxf_stage' }

    protected String getUnstageCommand() { 'nxf_unstage' }

    protected String getUnstageControls() {
        def result = copyFileToWorkDir(TaskRun.CMD_OUTFILE) + ' || true' + ENDL
        result += copyFileToWorkDir(TaskRun.CMD_ERRFILE) + ' || true' + ENDL
        if( statsEnabled )
            result += copyFileToWorkDir(TaskRun.CMD_TRACE) + ' || true' + ENDL
        if( outputEnvNames || outputEvals )
            result += copyFileToWorkDir(TaskRun.CMD_ENV) + ' || true' + ENDL
        return result
    }

}
