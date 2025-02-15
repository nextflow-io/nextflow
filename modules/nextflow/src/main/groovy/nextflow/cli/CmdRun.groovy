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

package nextflow.cli

import static org.fusesource.jansi.Ansi.*

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.regex.Pattern

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.IStringConverter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import groovyx.gpars.GParsConfig
import nextflow.BuildInfo
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.SysEnv
import nextflow.config.ConfigBuilder
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import nextflow.util.CustomPoolFactory
import nextflow.util.Duration
import nextflow.util.HistoryFile
import org.apache.commons.lang.StringUtils
import org.fusesource.jansi.AnsiConsole
import org.yaml.snakeyaml.Yaml
/**
 * CLI sub-command RUN
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Execute a pipeline project")
class CmdRun extends CmdBase implements HubOptions {

    static final public Pattern RUN_NAME_PATTERN = Pattern.compile(/^[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}$/, Pattern.CASE_INSENSITIVE)

    static final public List<String> VALID_PARAMS_FILE = ['json', 'yml', 'yaml']

    static final public String DSL2 = '2'
    static final public String DSL1 = '1'

    static {
        // install the custom pool factory for GPars threads
        GParsConfig.poolFactory = new CustomPoolFactory()
    }

    static class DurationConverter implements IStringConverter<Long> {
        @Override
        Long convert(String value) {
            if( !value ) throw new IllegalArgumentException()
            if( value.isLong() ) {  return value.toLong() }
            return Duration.of(value).toMillis()
        }
    }

    static final public String NAME = 'run'

    private Map<String,String> sysEnv = System.getenv()

    @Parameter(names=['-name'], description = 'Assign a mnemonic name to the a pipeline run')
    String runName

    @Parameter(names=['-lib'], description = 'Library extension path')
    String libPath

    @Parameter(names=['-cache'], description = 'Enable/disable processes caching', arity = 1)
    Boolean cacheable

    @Parameter(names=['-resume'], description = 'Execute the script using the cached results, useful to continue executions that was stopped by an error')
    String resume

    @Parameter(names=['-ps','-pool-size'], description = 'Number of threads in the execution pool', hidden = true)
    Integer poolSize

    @Parameter(names=['-pi','-poll-interval'], description = 'Executor poll interval (duration string ending with ms|s|m)', converter = DurationConverter, hidden = true)
    long pollInterval

    @Parameter(names=['-qs','-queue-size'], description = 'Max number of processes that can be executed in parallel by each executor')
    Integer queueSize

    @Parameter(names=['-test'], description = 'Test a script function with the name specified')
    String test

    @Parameter(names=['-o', '-output-dir'], description = 'Directory where workflow outputs are stored')
    String outputDir

    @Parameter(names=['-w', '-work-dir'], description = 'Directory where intermediate result files are stored')
    String workDir

    @Parameter(names=['-bucket-dir'], description = 'Remote bucket where intermediate result files are stored')
    String bucketDir

    @Parameter(names=['-with-cloudcache'], description = 'Enable the use of object storage bucket as storage for cache meta-data')
    String cloudCachePath

    /**
     * Defines the parameters to be passed to the pipeline script
     */
    @DynamicParameter(names = '--', description = 'Set a parameter used by the pipeline', hidden = true)
    Map<String,String> params = new LinkedHashMap<>()

    @Parameter(names='-params-file', description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @DynamicParameter(names = ['-process.'], description = 'Set process options' )
    Map<String,String> process = [:]

    @DynamicParameter(names = ['-e.'], description = 'Add the specified variable to execution environment')
    Map<String,String> env = [:]

    @Parameter(names = ['-E'], description = 'Exports all current system environment')
    boolean exportSysEnv

    @DynamicParameter(names = ['-executor.'], description = 'Set executor options', hidden = true )
    Map<String,String> executorOptions = [:]

    @Parameter(description = 'Project name or repository url')
    List<String> args

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names=['-d','-deep'], description = 'Create a shallow clone of the specified depth')
    Integer deep

    @Parameter(names=['-latest'], description = 'Pull latest changes before run')
    boolean latest

    @Parameter(names='-stdin', hidden = true)
    boolean stdin

    @Parameter(names = ['-ansi'], hidden = true, arity = 0)
    void setAnsi(boolean value) {
        launcher.options.ansiLog = value
    }

    @Parameter(names = ['-ansi-log'], description = 'Enable/disable ANSI console logging', arity = 1)
    void setAnsiLog(boolean value) {
        launcher.options.ansiLog = value
    }

    @Parameter(names = ['-with-tower'], description = 'Monitor workflow execution with Seqera Platform (formerly Tower Cloud)')
    String withTower

    @Parameter(names = ['-with-wave'], hidden = true)
    String withWave

    @Parameter(names = ['-with-fusion'], hidden = true)
    String withFusion

    @Parameter(names = ['-with-weblog'], description = 'Send workflow status messages via HTTP to target URL')
    String withWebLog

    @Parameter(names = ['-with-trace'], description = 'Create processes execution tracing file')
    String withTrace

    @Parameter(names = ['-with-report'], description = 'Create processes execution html report')
    String withReport

    @Parameter(names = ['-with-timeline'], description = 'Create processes execution timeline file')
    String withTimeline

    @Parameter(names = '-with-charliecloud', description = 'Enable process execution in a Charliecloud container runtime')
    def withCharliecloud

    @Parameter(names = '-with-singularity', description = 'Enable process execution in a Singularity container')
    def withSingularity

    @Parameter(names = '-with-apptainer', description = 'Enable process execution in a Apptainer container')
    def withApptainer

    @Parameter(names = '-with-podman', description = 'Enable process execution in a Podman container')
    def withPodman

    @Parameter(names = '-without-podman', description = 'Disable process execution in a Podman container')
    def withoutPodman

    @Parameter(names = '-with-docker', description = 'Enable process execution in a Docker container')
    def withDocker

    @Parameter(names = '-without-docker', description = 'Disable process execution with Docker', arity = 0)
    boolean withoutDocker

    @Parameter(names = '-with-mpi', hidden = true)
    boolean withMpi

    @Parameter(names = '-with-dag', description = 'Create pipeline DAG file')
    String withDag

    @Parameter(names = ['-bg'], arity = 0, hidden = true)
    void setBackground(boolean value) {
        launcher.options.background = value
    }

    @Parameter(names=['-c','-config'], hidden = true )
    List<String> runConfig

    @DynamicParameter(names = ['-cluster.'], description = 'Set cluster options', hidden = true )
    Map<String,String> clusterOptions = [:]

    @Parameter(names=['-profile'], description = 'Choose a configuration profile')
    String profile

    @Parameter(names=['-dump-hashes'], description = 'Dump task hash keys for debugging purpose')
    String dumpHashes

    @Parameter(names=['-dump-channels'], description = 'Dump channels for debugging purpose')
    String dumpChannels

    @Parameter(names=['-N','-with-notification'], description = 'Send a notification email on workflow completion to the specified recipients')
    String withNotification

    @Parameter(names=['-with-conda'], description = 'Use the specified Conda environment package or file (must end with .yml|.yaml suffix)')
    String withConda

    @Parameter(names=['-without-conda'], description = 'Disable the use of Conda environments')
    Boolean withoutConda

    @Parameter(names=['-with-spack'], description = 'Use the specified Spack environment package or file (must end with .yaml suffix)')
    String withSpack

    @Parameter(names=['-without-spack'], description = 'Disable the use of Spack environments')
    Boolean withoutSpack

    @Parameter(names=['-offline'], description = 'Do not check for remote project updates')
    boolean offline = System.getenv('NXF_OFFLINE')=='true'

    @Parameter(names=['-entry'], description = 'Entry workflow name to be executed', arity = 1)
    String entryName

    @Parameter(names=['-main-script'], description = 'The script file to be executed when launching a project directory or repository' )
    String mainScript

    @Parameter(names=['-stub-run','-stub'], description = 'Execute the workflow replacing process scripts with command stubs')
    boolean stubRun

    @Parameter(names=['-preview'], description = "Run the workflow script skipping the execution of all processes")
    boolean preview

    @Parameter(names=['-plugins'], description = 'Specify the plugins to be applied for this run e.g. nf-amazon,nf-tower')
    String plugins

    @Parameter(names=['-disable-jobs-cancellation'], description = 'Prevent the cancellation of child jobs on execution termination')
    Boolean disableJobsCancellation

    Boolean skipHistoryFile

    Boolean getDisableJobsCancellation() {
        return disableJobsCancellation!=null
                ?  disableJobsCancellation
                : sysEnv.get('NXF_DISABLE_JOBS_CANCELLATION') as boolean
    }

    /**
     * Optional closure modelling an action to be invoked when the preview mode is enabled
     */
    Closure<Void> previewAction

    @Override
    String getName() { NAME }

    String getParamsFile() {
        return paramsFile ?: sysEnv.get('NXF_PARAMS_FILE')
    }

    boolean hasParams() {
        return params || getParamsFile()
    }

    @Override
    void run() {
        final scriptArgs = (args?.size()>1 ? args[1..-1] : []) as List<String>
        final pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new AbortOperationException("No project name was specified")

        if( withPodman && withoutPodman )
            throw new AbortOperationException("Command line options `-with-podman` and `-without-podman` cannot be specified at the same time")

        if( withDocker && withoutDocker )
            throw new AbortOperationException("Command line options `-with-docker` and `-without-docker` cannot be specified at the same time")

        if( withConda && withoutConda )
            throw new AbortOperationException("Command line options `-with-conda` and `-without-conda` cannot be specified at the same time")

        if( withSpack && withoutSpack )
            throw new AbortOperationException("Command line options `-with-spack` and `-without-spack` cannot be specified at the same time")

        if( offline && latest )
            throw new AbortOperationException("Command line options `-latest` and `-offline` cannot be specified at the same time")

        checkRunName()

        printBanner()
        Plugins.init()

        // -- specify the arguments
        final scriptFile = getScriptFile(pipeline)

        // create the config object
        final builder = new ConfigBuilder()
                .setOptions(launcher.options)
                .setCmdRun(this)
                .setBaseDir(scriptFile.parent)
        final config = builder .build()

        // check DSL syntax in the config
        launchInfo(config, scriptFile)

        // check if NXF_ variables are set in nextflow.config
        checkConfigEnv(config)

        // -- load plugins
        final cfg = plugins ? [plugins: plugins.tokenize(',')] : config
        Plugins.load(cfg)

        // -- create a new runner instance
        final runner = new ScriptRunner(config)
        runner.setScript(scriptFile)
        runner.setPreview(this.preview, previewAction)
        runner.session.profile = profile
        runner.session.commandLine = launcher.cliString
        runner.session.ansiLog = launcher.options.ansiLog
        runner.session.debug = launcher.options.remoteDebug
        runner.session.disableJobsCancellation = getDisableJobsCancellation()

        final isTowerEnabled = config.navigate('tower.enabled') as Boolean
        if( isTowerEnabled || log.isTraceEnabled() )
            runner.session.resolvedConfig = ConfigBuilder.resolveConfig(scriptFile.parent, this)
        // note config files are collected during the build process
        // this line should be after `ConfigBuilder#build`
        runner.session.configFiles = builder.parsedConfigFiles
        // set the commit id (if any)
        runner.session.commitId = scriptFile.commitId
        if( this.test ) {
            runner.test(this.test, scriptArgs)
            return
        }

        def info = CmdInfo.status( log.isTraceEnabled() )
        log.debug( '\n'+info )

        // -- add this run to the local history
        if( !skipHistoryFile )  {
            runner.verifyAndTrackHistory(launcher.cliString, runName)
        }

        // -- run it!
        runner.execute(scriptArgs, this.entryName)
    }

    protected void printBanner() {
        if( launcher.options.ansiLog ){
            // Plain header for verbose log
            log.debug "N E X T F L O W  ~  version ${BuildInfo.version}"

            // Fancy coloured header for the ANSI console output
            def fmt = ansi()
            fmt.a("\n")
            // Use exact colour codes so that they render the same on every terminal,
            //   irrespective of terminal colour scheme.
            // Nextflow green RGB (13, 192, 157) and exact black text (0,0,0),
            //   Apple Terminal only supports 256 colours, so use the closest match:
            //   light sea green | #20B2AA | 38;5;0
            //   Don't use black for text as terminals mess with this in their colour schemes.
            //   Use very dark grey, which is more reliable.
            // Jansi library bundled in Jline can't do exact RGBs,
            //   so just do the ANSI codes manually
            final BACKGROUND = "\033[1m\033[38;5;232m\033[48;5;43m"
            fmt.a("$BACKGROUND N E X T F L O W ").reset()

            // Show Nextflow version
            fmt.a(Attribute.INTENSITY_FAINT).a("  ~  ").reset().a("version " + BuildInfo.version).reset()
            fmt.a("\n")
            AnsiConsole.out.println(fmt.eraseLine())
        }
        else {
            // Plain header to the console if ANSI is disabled
            log.info "N E X T F L O W  ~  version ${BuildInfo.version}"
        }
    }

    protected checkConfigEnv(ConfigMap config) {
        // Warn about setting NXF_ environment variables within env config scope
        final env = config.env as Map<String, String>
        for( String name : env.keySet() ) {
            if( name.startsWith('NXF_') && name!='NXF_DEBUG' ) {
                final msg = "Nextflow variables must be defined in the launching environment - The following variable set in the config file is going to be ignored: '$name'"
                log.warn(msg)
            }
        }
    }

    protected void launchInfo(ConfigMap config, ScriptFile scriptFile) {
        // -- determine strict mode
        detectStrictFeature(config, sysEnv)
        // -- determine moduleBinary
        detectModuleBinaryFeature(config)
        // -- determine dsl mode
        final dsl = detectDslMode(config, scriptFile.main.text, sysEnv)
        NextflowMeta.instance.enableDsl(dsl)
        // -- show launch info
        final ver = NF.dsl2 ? DSL2 : DSL1
        final repo = scriptFile.repository ?: scriptFile.source.toString()
        final head = preview ? "* PREVIEW * $scriptFile.repository" : "Launching `$repo`"
        final revision = scriptFile.repository
            ? scriptFile.revisionInfo.toString()
            : scriptFile.getScriptId()?.substring(0,10)
        printLaunchInfo(ver, repo, head, revision)
    }

    static void detectModuleBinaryFeature(ConfigMap config) {
        final moduleBinaries = config.navigate('nextflow.enable.moduleBinaries', false)
        if( moduleBinaries ) {
            log.debug "Enabling module binaries"
            NextflowMeta.instance.moduleBinaries(true)
        }
    }

    static void detectStrictFeature(ConfigMap config, Map sysEnv) {
        final defStrict = sysEnv.get('NXF_ENABLE_STRICT') ?: false
        log
        final strictMode = config.navigate('nextflow.enable.strict', defStrict)
        if( strictMode ) {
            log.debug "Enabling nextflow strict mode"
            NextflowMeta.instance.strictMode(true)
        }
    }

    protected void printLaunchInfo(String ver, String repo, String head, String revision) {
        if( launcher.options.ansiLog ){
            log.debug "${head} [$runName] DSL${ver} - revision: ${revision}"

            def fmt = ansi()
            fmt.a("Launching").fg(Color.MAGENTA).a(" `$repo` ").reset()
            fmt.a(Attribute.INTENSITY_FAINT).a("[").reset()
            fmt.bold().fg(Color.CYAN).a(runName).reset()
            fmt.a(Attribute.INTENSITY_FAINT).a("]")
            fmt.a(" DSL${ver} - ")
            fmt.fg(Color.CYAN).a("revision: ").reset()
            fmt.fg(Color.CYAN).a(revision).reset()
            fmt.a("\n")
            AnsiConsole.out().println(fmt.eraseLine())
        }
        else {
            log.info "${head} [$runName] DSL${ver} - revision: ${revision}"
        }
    }

    static String detectDslMode(ConfigMap config, String scriptText, Map sysEnv) {
        // -- try determine DSL version from config file

        final dsl = config.navigate('nextflow.enable.dsl') as String

        // -- script can still override the DSL version
        final scriptDsl = NextflowMeta.checkDslMode(scriptText)
        if( scriptDsl ) {
            log.debug("Applied DSL=$scriptDsl from script declaration")
            return scriptDsl
        }
        else if( dsl ) {
            log.debug("Applied DSL=$dsl from config declaration")
            return dsl
        }
        // -- if still unknown try probing for DSL1
        if( NextflowMeta.probeDsl1(scriptText) ) {
            log.debug "Applied DSL=1 by probing script field"
            return DSL1
        }

        final envDsl = sysEnv.get('NXF_DEFAULT_DSL')
        if( envDsl ) {
            log.debug "Applied DSL=$envDsl from NXF_DEFAULT_DSL variable"
            return envDsl
        }
        else {
            log.debug "Applied DSL=2 by global default"
            return DSL2
        }
    }

    protected void checkRunName() {
        if( runName == 'last' )
            throw new AbortOperationException("Not a valid run name: `last`")
        if( runName && !matchRunName(runName) )
            throw new AbortOperationException("Not a valid run name: `$runName` -- It must match the pattern $RUN_NAME_PATTERN")

        if( !runName ) {
            if( HistoryFile.disabled() )
                throw new AbortOperationException("Missing workflow run name")
            // -- make sure the generated name does not exist already
            runName = HistoryFile.DEFAULT.generateNextName()
        }

        else if( !HistoryFile.disabled() && HistoryFile.DEFAULT.checkExistsByName(runName) )
            throw new AbortOperationException("Run name `$runName` has been already used -- Specify a different one")
    }

    static protected boolean matchRunName(String name) {
        RUN_NAME_PATTERN.matcher(name).matches()
    }

    protected ScriptFile getScriptFile(String pipelineName) {
        try {
            getScriptFile0(pipelineName)
        }
        catch (IllegalArgumentException | AbortOperationException e) {
            if( e.message.startsWith("Not a valid project name:") && !guessIsRepo(pipelineName)) {
                throw new AbortOperationException("Cannot find script file: $pipelineName")
            }
            else
                throw e
        }
    }

    static protected boolean guessIsRepo(String name) {
        if( FileHelper.getUrlProtocol(name) != null )
            return true
        if( name.startsWith('/') )
            return false
        if( name.startsWith('./') || name.startsWith('../') )
            return false
        if( name.endsWith('.nf') )
            return false
        if( name.count('/') != 1 )
            return false
        return true
    }

    protected ScriptFile getScriptFile0(String pipelineName) {
        assert pipelineName

        /*
         * read from the stdin
         */
        if( pipelineName == '-' ) {
            def file = tryReadFromStdin()
            if( !file )
                throw new AbortOperationException("Cannot access `stdin` stream")

            if( revision )
                throw new AbortOperationException("Revision option cannot be used when running a script from stdin")

            return new ScriptFile(file)
        }

        /*
         * look for a file with the specified pipeline name
         */
        def script = new File(pipelineName)
        if( script.isDirectory()  ) {
            script = mainScript ? new File(mainScript) : new AssetManager().setLocalPath(script).getMainScriptFile()
        }

        if( script.exists() ) {
            if( revision )
                throw new AbortOperationException("Revision option cannot be used when running a local script")
            return new ScriptFile(script)
        }

        /*
         * try to look for a pipeline in the repository
         */
        def manager = new AssetManager(pipelineName, this)
        def repo = manager.getProject()

        boolean checkForUpdate = true
        if( !manager.isRunnable() || latest ) {
            if( offline )
                throw new AbortOperationException("Unknown project `$repo` -- NOTE: automatic download from remote repositories is disabled")
            log.info "Pulling $repo ..."
            def result = manager.download(revision,deep)
            if( result )
                log.info " $result"
            checkForUpdate = false
        }
        // checkout requested revision
        try {
            manager.checkout(revision)
            manager.updateModules()
            final scriptFile = manager.getScriptFile(mainScript)
            if( checkForUpdate && !offline )
                manager.checkRemoteStatus(scriptFile.revisionInfo)
            // return the script file
            return scriptFile
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new AbortOperationException("Unknown error accessing project `$repo` -- Repository may be corrupted: ${manager.localPath}", e)
        }

    }

    static protected File tryReadFromStdin() {
        if( !System.in.available() )
            return null

        getScriptFromStream(System.in)
    }

    static protected File getScriptFromStream( InputStream input, String name = 'nextflow' ) {
        input != null
        File result = File.createTempFile(name, null)
        result.deleteOnExit()
        input.withReader { Reader reader -> result << reader }
        return result
    }

    @Memoized  // <-- avoid parse multiple times the same file and params
    Map parsedParams(Map configVars) {

        final result = [:]
        final file = getParamsFile()
        if( file ) {
            def path = validateParamsFile(file)
            def type = path.extension.toLowerCase() ?: null
            if( type == 'json' )
                readJsonFile(path, configVars, result)
            else if( type == 'yml' || type == 'yaml' )
                readYamlFile(path, configVars, result)
        }

        // set the CLI params
        if( !params )
            return result

        for( Map.Entry<String,String> entry : params ) {
            addParam( result, entry.key, entry.value )
        }
        return result
    }


    static final private Pattern DOT_ESCAPED = ~/\\\./
    static final private Pattern DOT_NOT_ESCAPED = ~/(?<!\\)\./

    static protected void addParam(Map params, String key, String value, List path=[], String fullKey=null) {
        if( !fullKey )
            fullKey = key
        final m = DOT_NOT_ESCAPED.matcher(key)
        if( m.find() ) {
            final p = m.start()
            final root = key.substring(0, p)
            if( !root ) throw new AbortOperationException("Invalid parameter name: $fullKey")
            path.add(root)
            def nested = params.get(root)
            if( nested == null ) {
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            else if( nested !instanceof Map ) {
                log.warn "Command line parameter --${path.join('.')} is overwritten by --${fullKey}"
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            addParam((Map)nested, key.substring(p+1), value, path, fullKey)
        }
        else {
            addParam0(params, key.replaceAll(DOT_ESCAPED,'.'), parseParamValue(value))
        }
    }

    static protected void addParam0(Map params, String key, Object value) {
        if( key.contains('-') )
            key = kebabToCamelCase(key)
        params.put(key, value)
    }

    static protected String kebabToCamelCase(String str) {
        final result = new StringBuilder()
        str.split('-').eachWithIndex { String entry, int i ->
            result << (i>0 ? StringUtils.capitalize(entry) : entry )
        }
        return result.toString()
    }

    static protected parseParamValue(String str) {
        if ( SysEnv.get('NXF_DISABLE_PARAMS_TYPE_DETECTION') )
            return str

        if ( str == null ) return null

        if ( str.toLowerCase() == 'true') return Boolean.TRUE
        if ( str.toLowerCase() == 'false' ) return Boolean.FALSE

        if ( str==~/-?\d+(\.\d+)?/ && str.isInteger() ) return str.toInteger()
        if ( str==~/-?\d+(\.\d+)?/ && str.isLong() ) return str.toLong()
        if ( str==~/-?\d+(\.\d+)?/ && str.isDouble() ) return str.toDouble()

        return str
    }

    private Path validateParamsFile(String file) {

        def result = FileHelper.asPath(file)
        def ext = result.getExtension()
        if( !VALID_PARAMS_FILE.contains(ext) )
            throw new AbortOperationException("Not a valid params file extension: $file -- It must be one of the following: ${VALID_PARAMS_FILE.join(',')}")

        return result
    }

    static private Pattern PARAMS_VAR = ~/(?m)\$\{(\p{javaJavaIdentifierStart}\p{javaJavaIdentifierPart}*)}/

    protected String replaceVars0(String content, Map binding) {
        content.replaceAll(PARAMS_VAR) { List<String> matcher ->
            // - the regex matcher is represented as list
            // - the first element is the matching string ie. `${something}`
            // - the second element is the group content ie. `something`
            // - make sure the regex contains at least a group otherwise the closure
            // parameter is a string instead of a list of the call fail
            final placeholder = matcher.get(0)
            final key = matcher.get(1)

            if( !binding.containsKey(key) ) {
                final msg = "Missing params file variable: $placeholder"
                if(NF.strictMode)
                    throw new AbortOperationException(msg)
                log.warn msg
                return placeholder
            }

            return binding.get(key)
        }
    }

    private void readJsonFile(Path file, Map configVars, Map result) {
        try {
            def text = configVars ? replaceVars0(file.text, configVars) : file.text
            def json = (Map)new JsonSlurper().parseText(text)
            result.putAll(json)
        }
        catch (NoSuchFileException | FileNotFoundException e) {
            throw new AbortOperationException("Specified params file does not exist: ${file.toUriString()}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()} - Cause: ${e.message}", e)
        }
    }

    private void readYamlFile(Path file, Map configVars, Map result) {
        try {
            def text = configVars ? replaceVars0(file.text, configVars) : file.text
            def yaml = (Map)new Yaml().load(text)
            result.putAll(yaml)
        }
        catch (NoSuchFileException | FileNotFoundException e) {
            throw new AbortOperationException("Specified params file does not exist: ${file.toUriString()}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()}", e)
        }
    }

}
