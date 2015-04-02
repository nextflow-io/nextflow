/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow.cli
import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.IStringConverter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.scm.AssetManager
import nextflow.script.ConfigBuilder
import nextflow.script.ScriptRunner
import nextflow.util.Duration
import nextflow.util.HistoryFile
/**
 * CLI sub-command RUN
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Launch a pipeline execution")
class CmdRun extends CmdBase implements HubOptions {

    static class DurationConverter implements IStringConverter<Long> {
        @Override
        Long convert(String value) {
            if( !value ) throw new IllegalArgumentException()
            if( value.isLong() ) {  return value.toLong() }
            return Duration.of(value).toMillis()
        }
    }

    static final NAME = 'run'

    @Parameter(names=['-lib'], description = 'Library extension path')
    String libPath

    @Parameter(names=['-cache'], description = 'enable/disable processes caching', arity = 1)
    boolean cacheable = true

    @Parameter(names=['-resume'], description = 'Execute the script using the cached results, useful to continue executions that stopped by an error')
    String resume

    @Parameter(names=['-ps','-pool-size'], description = 'Number of threads in the execution pool', hidden = true)
    Integer poolSize

    @Parameter(names=['-pi','-poll-interval'], description = 'Executor poll interval (duration string ending with ms|s|m)', converter = DurationConverter, hidden = true)
    long pollInterval

    @Parameter(names=['-qs','-queue-size'], description = 'Max number of processes that can be executed in parallel by each executor')
    Integer queueSize

    @Parameter(names=['-test'], description = 'Test function with the specified name')
    String test

    @Parameter(names=['-w', '-work-dir'], description = 'Directory where intermediate results are stored')
    String workDir = System.getenv('NXF_WORK') ?: 'work'

    /**
     * Defines the parameters to be passed to the pipeline script
     */
    @DynamicParameter(names = '--', description = 'Set a parameter used by the pipeline', hidden = true)
    Map<String,String> params = new LinkedHashMap<>()

    @DynamicParameter(names = ['-process.'], description = 'Set process default options' )
    Map<String,String> process = [:]

    @DynamicParameter(names = ['-e.'], description = 'Add the specified variable to execution environment')
    Map<String,String> env = [:]

    @Parameter(names = ['-E'], description = 'Exports all the current system environment')
    boolean exportSysEnv

    @DynamicParameter(names = ['-executor.'], description = 'Executor(s) options', hidden = true )
    Map<String,String> executorOptions = [:]

    @Parameter(description = 'name of pipeline to run')
    List<String> args

    @Parameter(names=['-r','-revision'], description = 'Revision of pipeline to run (either a branch, tag or commit SHA number)')
    String revision

    @Parameter(names=['-latest'], description = 'Pull latest changes before run')
    boolean latest

    @Parameter(names='-stdin', hidden = true)
    boolean stdin

    @Parameter(names = ['-with-extrae'], description = 'Trace execution by using BSC Extrae', arity = 0, hidden = true)
    boolean withExtrae

    @Parameter(names = ['-with-drmaa'], description = 'Enable DRMAA binding')
    String withDrmaa

    @Parameter(names = ['-with-trace'], description = 'Trace execution to the specified file')
    String withTrace

    @Parameter(names = '-with-docker', description = 'Enable process execution in a Docker container')
    def withDocker

    @Parameter(names = '-without-docker', description = 'Disable process execution with Docker', arity = 0)
    boolean withoutDocker

    @Parameter(names = ['-bg'], arity = 0, hidden = true)
    void setBackground(boolean value) {
        launcher.options.background = value
    }

    @Parameter(names=['-c','-config'], hidden = true )
    List<String> runConfig

    @DynamicParameter(names = ['-cluster.'], description = 'Define cluster options', hidden = true )
    Map<String,String> clusterOptions = [:]

    @Override
    final String getName() { NAME }

    @Override
    void run() {
        String pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new AbortOperationException("No pipeline to run was specified")

        if( withDocker && withoutDocker )
            throw new AbortOperationException("Command line options `-with-docker` and `-without-docker` cannot be specified at the same time")


        log.info "N E X T F L O W  ~  version ${Const.APP_VER}"

        // -- specify the arguments
        def scriptArgs = (args?.size()>1 ? args[1..-1] : []) as List<String>
        def scriptFile = getScriptFile(pipeline).toPath().complete()

        // create the config object
        def config = new ConfigBuilder()
                        .setOptions(launcher.options)
                        .setCmdRun(this)
                        .setBaseDir(scriptFile.parent)
                        .build()
                        .toMap()

        // -- create a new runner instance
        def runner = new ScriptRunner(config)
                            .setScript(scriptFile)
                            .setRunOptions(this)

        if( this.test ) {
            runner.test(this.test, scriptArgs )
        }
        else {
            def info = CmdInfo.status( log.isTraceEnabled() )
            log.debug( '\n'+info )

            // -- add this run to the local history
            trackHistory(runner.session.uniqueId)

            // -- run it!
            runner.execute(scriptArgs)
        }
    }

    protected void trackHistory(UUID uuid) {
        def cli = launcher.cliString
        if( !cli )
            return
        def p = cli.indexOf('nextflow ')
        if( p != -1 ) cli = 'nextflow ' + cli.substring(p+9)
        HistoryFile.history.write( uuid, cli )
    }

    protected File getScriptFile(String pipelineName) {
        assert pipelineName

        /*
         * read from the stdin
         */
        if( pipelineName == '-' ) {
            def file = tryReadFromStdin()
            if( !file )
                throw new AbortOperationException("Cannot pipeline script from stdin")

            if( revision )
                throw new AbortOperationException("Revision option cannot be used running a local pipeline script")

            return file
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
                throw new AbortOperationException("Revision option cannot be used running a local pipeline script")
            return script
        }

        /*
         * try to look for a pipeline in the repository
         */
        def manager = new AssetManager(hub: getHubProvider(), user: getHubUser(), pwd: getHubPassword())
        def repo = manager.resolveName(pipelineName)

        /*
         * set the required name
         */
        manager.setPipeline(repo as String)

        if( !manager.isLocal() || latest ) {
            log.info "Pulling $repo ..."
            def result = manager.download()
            if( result )
                log.info " $result"
        }
        // checkout requested revision
        manager.checkout(revision)

        // return the script file
        return manager.getMainScriptFile()
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


}
