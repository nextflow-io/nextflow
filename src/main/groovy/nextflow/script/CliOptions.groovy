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
    String workDir = './work'

    /**
     * Print out the version number and exit
     */
    @Parameter(names = ['-v','-version'], description = 'Show the program version')
    boolean version

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
        }

        return normalized
    }

}
