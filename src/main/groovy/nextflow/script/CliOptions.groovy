package nextflow.script

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import nextflow.Const

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
    String logFile = ".${Const.APP_NAME}.log"

    @Parameter(names=['-lib'], description = 'Library extension path')
    String libPath

    /**
     * The packages to trace
     */
    @Parameter(hidden = true, names='-trace')
    List<String> trace

    @Parameter(names=['-c','-config'], description = 'Use the specified configuration file(s)')
    List<String> config

    @Parameter(names=['-cache'], description = 'Enable/disable task(s) caching', arity = 1)
    boolean cacheable = true

    @Parameter(names=['-resume'], description = 'Execute the script using the cached results, useful to continue executions that stopped by an error')
    String resume

    @Parameter(names=['-ps','-pool-size'], description = 'The number of threads in the executor pool')
    Integer poolSize

    @Parameter(names=['-test'], description = 'Test the function with the name specified')
    String test

    @Parameter(names=['-w', '-work-dir'], description = 'Directory where tasks results are stored')
    File workDir = new File('./work')

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
    Map<String,String> process = new LinkedHashMap<>()

    /**
     * Extra parameters for the script execution
     */
    @Parameter(description = 'Script level arguments')
    List<String> arguments
}
