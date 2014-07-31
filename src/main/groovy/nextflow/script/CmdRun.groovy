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

package nextflow.script
import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.IStringConverter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.util.logging.Slf4j
import nextflow.exception.MissingScriptException
import nextflow.share.PipelineManager
import nextflow.util.Duration
import nextflow.util.FileHelper
import nextflow.util.HistoryFile
/**
 * CLI sub-command RUN
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Parameters(commandDescription = "Run a pipeline")
class CmdRun implements CmdX {

    static class DurationConverter implements IStringConverter<Long> {
        @Override
        Long convert(String value) {
            if( !value ) throw new IllegalArgumentException()
            if( value.isLong() ) {  return value.toLong() }
            return Duration.of(value).toMillis()
        }
    }


    @Parameter(names=['-cache'], description = 'Enable/disable processes caching', arity = 1)
    boolean cacheable = true

    @Parameter(names=['-resume'], description = 'Execute the script using the cached results, useful to continue executions that stopped by an error')
    String resume

    @Parameter(names=['-c','-config'], description = 'Use the specified configuration file(s)')
    List<String> config

    @Parameter(names=['-ps','-pool-size'], description = 'The number of threads in the execution pool', hidden = true)
    Integer poolSize

    @Parameter(names=['-pi','-poll-interval'], description = 'The executor poll interval (duration string ending with ms|s|m)', converter = DurationConverter)
    long pollInterval

    @Parameter(names=['-qs','-queue-size'], description = 'The max number of processes that can be executed in parallel by each executor')
    Integer queueSize

    @Parameter(names=['-test'], description = 'Test the function with the name specified')
    String test

    @Parameter(names=['-w', '-work-dir'], description = 'Directory where intermediate results are stored')
    String workDir = 'work'

    /**
     * Defines the parameters to be passed to the pipeline script
     */
    @DynamicParameter(names = '--', description = 'Set a parameter used by the workflow' )
    Map<String,String> params = new LinkedHashMap<>()

    @DynamicParameter(names = ['-process.'], description = 'Set default process options' )
    Map<String,String> process = [:]

    @DynamicParameter(names = ['-e.'], description = 'Add the specified variable to execution environment')
    Map<String,String> env = [:]

    @Parameter(names = ['-E'], description = 'Exports all the current system environment')
    boolean exportSysEnv

    @DynamicParameter(names = ['-executor.'], description = 'Executor(s) options', hidden = true )
    Map<String,String> executorOptions = [:]

    @Parameter(names = ['-bg'], description = 'Launch the pipeline as a background job', arity = 0)
    boolean background

    @Parameter(required=true, description = 'The name of the pipeline to run')
    List<String> args

    @Parameter(names=['-r','-revision'], description = 'Specify the pipeline revision to run')
    String revision

    String pipeline

    @Override
    void run() {
        pipeline = args[0]

        // -- specify the arguments
        def scriptArgs = args?.size()>1 ? args[1..-1] : []
        def scriptFile = getScriptFile()

        // -- check file system providers
        FileHelper.checkFileSystemProviders()

        // create the config object
        def config = new ConfigBuilder().setOptions(options).setCmdRun(this).build()

        // -- create a new runner instance
        def runner = new CliRunner(config)
        runner.session.cacheable = this.cacheable
        runner.session.resumeMode = this.resume != null
        // note -- make sure to use 'FileHelper.asPath' since it guarantee to handle correctly non-standard file system e.g. 'dxfs'
        runner.session.workDir = FileHelper.asPath(this.workDir).toAbsolutePath()
        runner.session.baseDir = scriptFile?.canonicalFile?.parentFile
        runner.libPath = options.libPath

        log.debug "Script bin dir: ${runner.session.binDir}"

        if( this.test ) {
            runner.test(scriptFile, this.test, scriptArgs )
        }
        else {
            log.debug( '\n'+CmdInfo.getInfo() )
            // -- add this run to the local history
            HistoryFile.history.append( runner.session.uniqueId, args )
            // -- run it!
            runner.execute(scriptFile,scriptArgs)
        }
    }


    protected File getScriptFile() {
        assert pipeline

        /*
         * read from the stdin
         */
        if( pipeline == '-' ) {
            def file = tryReadFromStdin()
            if( !file )
                throw new MissingScriptException(scriptFile: 'stdin')

            if( revision )
                throw new IllegalArgumentException("Revision option can be used running a local pipeline script")

            return file
        }

        /*
         * look for a file with the specified pipeline name
         */
        def script = new File(pipeline)
        if( script.isDirectory()  ) {
            script = new File(pipeline, PipelineManager.MAIN_FILE_NAME)
        }

        if( script.exists() ) {
            if( revision )
                throw new IllegalArgumentException("Revision option can be used running a local pipeline script")
            return script
        }

        /*
         * try to look for a pipeline in the repository
         */
        def manager = new PipelineManager(pipeline)
        script = manager.getMainFile()

        if( !script.exists() ) {
            log.info "Pull $pipeline ..."
            manager.pull()
        }

        if( script.exists() ) {
            manager.checkout(revision)
            return script
        }

        throw new MissingScriptException(scriptFile: script)
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
        input.withReader { reader -> result << reader }
        return result
    }

    static protected File getScriptFromUrl( String urlOrPath ) {
        def lower = urlOrPath.toLowerCase()
        def isUrl = ['http','https','ftp'].any { lower.startsWith(it+'://') }

        if( isUrl ) {
            def url = new URL(urlOrPath)
            def fileName = new File(url.getPath()).getBaseName()
            return getScriptFromStream( url.newInputStream(), fileName )
        }

        def file = new File(urlOrPath)
        if( !file.exists() ) {
            throw new FileNotFoundException("File do not exist: $file")
        }
        return file
    }


}
