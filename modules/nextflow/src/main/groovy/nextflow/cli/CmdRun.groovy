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

package nextflow.cli
import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.IStringConverter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.GParsConfig
import nextflow.Const
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.scm.AssetManager
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import nextflow.util.CustomPoolFactory
import nextflow.util.Duration
import nextflow.util.HistoryFile
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


    static List<String> VALID_PARAMS_FILE = ['json', 'yml', 'yaml']

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

    static final public NAME = 'run'

    @Parameter(names=['-name'], description = 'Assign a mnemonic name to the a pipeline run')
    String runName

    @Parameter(names=['-lib'], description = 'Library extension path')
    String libPath

    @Parameter(names=['-cache'], description = 'Enable/disable processes caching', arity = 1)
    boolean cacheable = true

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

    @Parameter(names=['-w', '-work-dir'], description = 'Directory where intermediate result files are stored')
    String workDir

    @Parameter(names=['-bucket-dir'], description = 'Remote bucket where intermediate result files are stored')
    String bucketDir

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

    @Parameter(names=['-latest'], description = 'Pull latest changes before run')
    boolean latest

    @Parameter(names='-stdin', hidden = true)
    boolean stdin

    @Parameter(names = ['-ansi'], hidden = true, arity = 0)
    boolean setAnsi(boolean value) {
        launcher.options.ansiLog = value
    }

    @Parameter(names = ['-ansi-log'], description = 'Enable/disable ANSI console logging', arity = 1)
    boolean setAnsiLog(boolean value) {
        launcher.options.ansiLog = value
    }

    @Parameter(names = ['-with-weblog'], description = 'Send workflow status messages via HTTP to target URL')
    String withWebLog

    @Parameter(names = ['-with-trace'], description = 'Create processes execution tracing file')
    String withTrace

    @Parameter(names = ['-with-report'], description = 'Create processes execution html report')
    String withReport

    @Parameter(names = ['-with-timeline'], description = 'Create processes execution timeline file')
    String withTimeline

    @Parameter(names = '-with-singularity', description = 'Enable process execution in a Singularity container')
    def withSingularity

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
    boolean dumpHashes

    @Parameter(names=['-dump-channels'], description = 'Dump channels for debugging purpose')
    String dumpChannels

    @Parameter(names=['-N','-with-notification'], description = 'Send a notification email on workflow completion to the specified recipients')
    String withNotification

    @Parameter(names=['-with-conda'], description = 'Use the specified Conda environment package or file (must end with .yml|.yaml suffix)')
    String withConda

    @Parameter(names=['-offline'], description = 'Do not check for remote project updates')
    boolean offline = System.getenv('NXF_OFFLINE') as boolean

    @Override
    String getName() { NAME }

    @Override
    void run() {
        final scriptArgs = (args?.size()>1 ? args[1..-1] : []) as List<String>
        final pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new AbortOperationException("No project name was specified")

        if( withDocker && withoutDocker )
            throw new AbortOperationException("Command line options `-with-docker` and `-without-docker` cannot be specified at the same time")

        if( offline && latest )
            throw new AbortOperationException("Command line options `-latest` and `-offline` cannot be specified at the same time")

        checkRunName()

        log.info "N E X T F L O W  ~  version ${Const.APP_VER}"

        // -- specify the arguments
        final scriptFile = getScriptFile(pipeline)

        // create the config object
        final config = new ConfigBuilder()
                        .setOptions(launcher.options)
                        .setCmdRun(this)
                        .setBaseDir(scriptFile.parent)

        // -- create a new runner instance
        final runner = new ScriptRunner(config)
        runner.script = scriptFile
        runner.profile = profile
        runner.session.ansiLog = launcher.options.ansiLog

        if( this.test ) {
            runner.test(this.test, scriptArgs)
            return
        }

        def info = CmdInfo.status( log.isTraceEnabled() )
        log.debug( '\n'+info )

        // -- add this run to the local history
        runner.verifyAndTrackHistory(launcher.cliString, runName)

        // -- run it!
        runner.execute(scriptArgs)
    }

    protected void checkRunName() {
        if( runName == 'last' )
            throw new AbortOperationException("Not a valid run name: `last`")

        if( !runName ) {
            // -- make sure the generated name does not exist already
            runName = HistoryFile.DEFAULT.generateNextName()
        }

        else if( HistoryFile.DEFAULT.checkExistsByName(runName) )
            throw new AbortOperationException("Run name `$runName` has been already used -- Specify a different one")
    }

    protected ScriptFile getScriptFile(String pipelineName) {
        assert pipelineName

        /*
         * read from the stdin
         */
        if( pipelineName == '-' ) {
            def file = tryReadFromStdin()
            if( !file )
                throw new AbortOperationException("Cannot access `stdin` stream")

            if( revision )
                throw new AbortOperationException("Revision option cannot be used running a local script")

            return new ScriptFile(file)
        }

        /*
         * look for a file with the specified pipeline name
         */
        def script = new File(pipelineName)
        if( script.isDirectory()  ) {
            script = new AssetManager().setLocalPath(script).getMainScriptFile()
        }

        if( script.exists() ) {
            if( revision )
                throw new AbortOperationException("Revision option cannot be used running a script")
            def result = new ScriptFile(script)
            log.info "Launching `$script` [$runName] - revision: ${result.getScriptId()?.substring(0,10)}"
            return result
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
            def result = manager.download()
            if( result )
                log.info " $result"
            checkForUpdate = false
        }
        // checkout requested revision
        try {
            manager.checkout(revision)
            manager.updateModules()
            def scriptFile = manager.getScriptFile()
            log.info "Launching `$repo` [$runName] - revision: ${scriptFile.revisionInfo}"
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

    Map getParsedParams() {

        def result = [:]

        if( paramsFile ) {
            def path = validateParamsFile(paramsFile)
            def ext = path.extension.toLowerCase() ?: null
            if( ext == 'json' )
                readJsonFile(path, result)
            else if( ext == 'yml' || ext == 'yaml' )
                readYamlFile(path, result)
        }

        // read the params file if any

        // set the CLI params
        params?.each { key, value ->
            result.put( key, parseParam(value) )
        }
        return result
    }

    static private parseParam( String str ) {

        if ( str == null ) return null

        if ( str.toLowerCase() == 'true') return Boolean.TRUE
        if ( str.toLowerCase() == 'false' ) return Boolean.FALSE

        if ( str.isInteger() ) return str.toInteger()
        if ( str.isLong() ) return str.toLong()
        if ( str.isDouble() ) return str.toDouble()

        return str
    }

    private Path validateParamsFile(String file) {

        def result = FileHelper.asPath(file)
        if( !result.exists() )
            throw new AbortOperationException("Specified params file does not exists: $file")

        def ext = result.getExtension()
        if( !VALID_PARAMS_FILE.contains(ext) )
            throw new AbortOperationException("Not a valid params file extension: $file -- It must be one of the following: ${VALID_PARAMS_FILE.join(',')}")

        return result
    }


    private void readJsonFile(Path file, Map result) {
        try {
            def json = (Map)new JsonSlurper().parse(Files.newInputStream(file))
            result.putAll(json)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: $file", e)
        }
    }

    private void readYamlFile(Path file, Map result) {
        try {
            def yaml = (Map)new Yaml().load(Files.newInputStream(file))
            result.putAll(yaml)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: $file", e)
        }
    }

}
