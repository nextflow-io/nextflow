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
import groovy.transform.PackageScope
import nextflow.util.Duration
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CliOptions {

    /**
     * The packages to debug
     */
    @Parameter(hidden = true, names='-debug')
    List<String> debug

    @Parameter(names=['-history'], description = 'Show history of executed commands')
    boolean history

    @Parameter(names=['-log'], description = 'Define the application log file')
    String logFile

    @Parameter(names=['-lib'], description = 'Library extension path')
    String libPath

    /**
     * The packages to trace
     */
    @Parameter(hidden = true, names='-trace')
    List<String> trace

    @Parameter(names=['-c','-config'], description = 'Use the specified configuration file(s)')
    List<String> config

    @Parameter(names=['-cache'], description = 'Enable/disable processes caching', arity = 1)
    boolean cacheable = true

    @Parameter(names=['-resume'], description = 'Execute the script using the cached results, useful to continue executions that stopped by an error')
    String resume

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
     * Print out the version number and exit
     */
    @Parameter(names = ['-v','-version'], description = 'Show the program version')
    boolean version

    @Parameter(names = ['-info'], description = 'The runtime info', hidden = true)
    boolean info

    /**
     * Print out the 'help' and exit
     */
    @Parameter(names = ['-h','-help'], description = 'Show this help', help = true)
    boolean help

    @Parameter(names = ['-q','-quiet'], description = 'Do not print information messages' )
    boolean quiet

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

    @Parameter(names = ['-d'], description = 'Start in cluster daemon mode', arity = 0)
    boolean daemon

    @DynamicParameter(names = ['-daemon.'], description = 'Starts in daemon mode and provides extra options', hidden = true )
    Map<String,String> daemonOptions = [:]

    @DynamicParameter(names = ['-executor.'], description = 'Executor(s) options', hidden = true )
    Map<String,String> executorOptions = [:]

    @Parameter(names = ['-bg'], description = 'Launch the pipeline as a background job', arity = 0)
    boolean background

    @Parameter(names = ['-with-extrae'], description = 'Trace the pipeline by using BSC Extrae', arity = 0, hidden = true)
    boolean withExtrae

    @Parameter(names = ['-trace-file'], description = 'Save the process execution time to the specified file')
    String traceFile


    /**
     * Extra parameters for the script execution
     */
    @Parameter(description = 'Script level arguments')
    List<String> arguments


    boolean isDaemon() {
        return daemon || daemonOptions.size() > 0
    }


    static class DurationConverter implements IStringConverter<Long> {

        @Override
        Long convert(String value) {
            if( !value ) throw new IllegalArgumentException()
            if( value.isLong() ) {  return value.toLong() }
            return Duration.of(value).toMillis()
        }
    }


    @PackageScope
    static List<String> normalizeArgs( String ... args ) {

        def normalized = []
        int i=0
        while( true ) {
            if( i==args.size() ) { break }

            def current = args[i++]
            normalized << current

            if( current == '-resume' ) {
                if( i<args.size() && !args[i].startsWith('-') && (args[i]=='last' || args[i] =~~ /[0-9a-f]{8}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{8}/) ) {
                    normalized << args[i++]
                }
                else {
                    normalized << 'last'
                }
            }
            else if( current == '-test' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '%all'
            }

            else if( current ==~ /^\-\-[a-zA-Z\d].*/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-process\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-daemon\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-executor\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }
        }

        return normalized
    }

}
