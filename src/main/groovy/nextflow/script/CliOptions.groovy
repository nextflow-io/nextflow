package nextflow.script

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter

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

    /**
     * The packages to trace
     */
    @Parameter(hidden = true, names='-trace')
    List<String> trace

    @Parameter(names=['-c','-config'], description = 'Use the specified configuration file(s)')
    List<String> config

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

    @Parameter(names= ['-wd','-work-dir'], description =  'Directory where output files are stored')
    String workDirectory

    /**
     * Defines the parameters to be passed to the pipeline script
     */
    @DynamicParameter(names = "--", description = "Define a variable to be used by the workflow" )
    Map<String,String> params = new LinkedHashMap<>()

    /**
     * Extra parameters for the script execution
     */
    @Parameter(description = 'Script level arguments')
    List<String> arguments
}
