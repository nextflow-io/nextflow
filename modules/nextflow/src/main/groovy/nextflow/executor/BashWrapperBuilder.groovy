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

package nextflow.executor

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.container.ContainerBuilder
import nextflow.container.DockerBuilder
import nextflow.container.PodmanBuilder
import nextflow.container.ShifterBuilder
import nextflow.container.SingularityBuilder
import nextflow.container.UdockerBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Escape

/**
 * Builder to create the BASH script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BashWrapperBuilder {

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
         * - 1: dump current environment in the `.command.log` file
         * - 2: trace the execution of user script adding the `set -x` flag
         * - 3: trace the execution of wrapper scripts
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
        return workDir != targetDir
    }

    protected boolean fixOwnership() {
        systemOsName == 'Linux' && containerConfig?.fixOwnership && runWithContainer && containerConfig.engine == 'docker' // <-- note: only for docker (shifter is not affected)
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

    protected String getOutputEnvCaptureSnippet(List<String> names) {
        def result = new StringBuilder()
        result.append('\n')
        result.append('# capture process environment\n')
        result.append('set +u\n')
        for( int i=0; i<names.size(); i++) {
            final key = names[i]
            result.append "echo $key=\$$key "
            result.append( i==0 ? '> ' : '>> ' )
            result.append(TaskRun.CMD_ENV)
            result.append('\n')
        }
        result.toString()
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

        if( outputEnvNames ) {
            if( !isBash(interpreter) ) throw new IllegalArgumentException("Process output of type env is only allowed with Bash process command -- Current interpreter: $interpreter")
            script += getOutputEnvCaptureSnippet(outputEnvNames)
        }

        final binding = new HashMap<String,String>(20)
        binding.header_script = headerScript
        binding.task_name = name
        binding.helpers_script = getHelpersScript()

        if( runWithContainer ) {
            binding.container_boxid = 'export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"'
            binding.container_helpers = containerBuilder.getScriptHelpers()
            binding.kill_cmd = containerBuilder.getKillCommand()
        }
        else {
            binding.container_boxid = null
            binding.container_helpers = null
            binding.kill_cmd = KILL_CMD
        }

        binding.cleanup_cmd = getCleanupCmd(changeDir)
        binding.scratch_cmd = ( changeDir ?: "NXF_SCRATCH=''" )

        binding.exit_file = exitFile(exitedFile)
        binding.touch_file = touchFile(startedFile)

        binding.module_load = getModuleLoadSnippet()
        binding.before_script = getBeforeScriptSnippet()
        binding.conda_activate = getCondaActivateSnippet()

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
         * staging input files when required
         */
        final stagingScript = copyStrategy.getStageInputFilesScript(inputFiles)
        binding.stage_inputs = stagingScript ? "# stage input files\n${stagingScript}" : null

        binding.stdout_file = TaskRun.CMD_OUTFILE
        binding.stderr_file = TaskRun.CMD_ERRFILE
        binding.trace_file = TaskRun.CMD_TRACE

        binding.trace_cmd = getTraceCommand(interpreter)
        binding.launch_cmd = getLaunchCommand(interpreter,env)
        binding.stage_cmd = getStageCommand()
        binding.unstage_cmd = getUnstageCommand()
        binding.unstage_controls = changeDir ? getUnstageControls() : null

        if( changeDir || shouldUnstageOutputs() ) {
            binding.unstage_outputs = copyStrategy.getUnstageOutputFilesScript(outputFiles,targetDir)
        }
        else {
            binding.unstage_outputs = null
        }

        binding.after_script = afterScript ? "# 'afterScript' directive\n$afterScript" : null

        // patch root ownership problem on files created with docker
        binding.fix_ownership = fixOwnership() ? "[ \${NXF_OWNER:=''} ] && chown -fR --from root \$NXF_OWNER ${workDir}/{*,.*} || true" : null

        binding.trace_script = isTraceRequired() ? getTraceScript(binding) : null
        
        return binding
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
        Files.write(wrapperFile, wrapper.getBytes())
        Files.write(scriptFile, script.getBytes())
        if( input != null )
            Files.write(inputFile, input.toString().getBytes())
        return wrapperFile
    }

    private String getHelpersScript() {
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
        def result = "# conda environment\n"
        result += 'source $(conda info --json | awk \'/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }\')'
        result += "/bin/activate ${Escape.path(condaEnv)}\n"
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

    protected String getLaunchCommand(String interpreter, String env) {
        /*
        * process stats
        */
        String launcher
        final traceWrapper = isTraceRequired()
        if( traceWrapper ) {
            // executes the stub which in turn executes the target command
            launcher = "/bin/bash ${fileStr(wrapperFile)} nxf_trace"
        }
        else {
            launcher = "${interpreter} ${fileStr(scriptFile)}"
        }

        /*
         * create the container engine command when needed
         */
        if( containerBuilder ) {
            def cmd = env ? 'eval $(nxf_container_env); ' + launcher : launcher
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
            result += (!containerBuilder ? 'rm -rf $NXF_SCRATCH || true' : '(sudo -n true && sudo rm -rf "$NXF_SCRATCH" || rm -rf "$NXF_SCRATCH")&>/dev/null || true')
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

    @PackageScope
    ContainerBuilder createContainerBuilder0(String engine) {
        /*
         * create a builder instance given the container engine
         */
        if( engine == 'docker' )
            return new DockerBuilder(containerImage)
        if( engine == 'podman' )
            return new PodmanBuilder(containerImage)
        if( engine == 'singularity' )
            return new SingularityBuilder(containerImage)
        if( engine == 'udocker' )
            return new UdockerBuilder(containerImage)
        if( engine == 'shifter' )
            return new ShifterBuilder(containerImage)
        //
        throw new IllegalArgumentException("Unknown container engine: $engine")
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
        if( stageInMode != 'copy' )
            builder.addMountForInputs(inputFiles)

        builder.addMount(binDir)

        if(this.containerMount)
            builder.addMount(containerMount)

        // task work dir
        builder.setWorkDir(workDir)

        // set the name
        builder.setName('$NXF_BOXID')

        if( this.containerMemory )
            builder.setMemory(containerMemory)

        if( this.containerCpuset )
            builder.addRunOptions(containerCpuset)

        // export the nextflow script debug variable
        if( isTraceRequired() )
            builder.addEnv( 'NXF_DEBUG=${NXF_DEBUG:=0}')

        // add the user owner variable in order to patch root owned files problem
        if( fixOwnership() )
            builder.addEnv( 'NXF_OWNER=$(id -u):$(id -g)' )

        if( engine=='docker' && System.getenv('NXF_DOCKER_OPTS') ) {
            builder.addRunOptions(System.getenv('NXF_DOCKER_OPTS'))
        }

        for( String var : containerConfig.getEnvWhitelist() ) {
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

        builder.params(entry: '/bin/bash')

        // give a chance to override any option with process specific `containerOptions`
        if( containerOptions ) {
            builder.addRunOptions(containerOptions)
        }

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
        if(  outputEnvNames )
            result += copyFileToWorkDir(TaskRun.CMD_ENV) + ' || true' + ENDL
        return result
    }

}
