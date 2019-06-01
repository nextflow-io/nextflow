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

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import org.fusesource.jansi.Ansi

/**
 * Main application command line options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CliOptions {

    /**
     * The packages to debug
     */
    @Parameter(hidden = true, names='-debug')
    List<String> debug

    @Parameter(names=['-log'], description = 'Set nextflow log file path')
    String logFile

    @Parameter(names=['-c','-config'], description = 'Add the specified file to configuration set')
    List<String> userConfig

    @Parameter(names=['-C'], description = 'Use the specified configuration file(s) overriding any defaults')
    List<String> config

    /**
     * the packages to trace
     */
    @Parameter(names='-trace', hidden = true)
    List<String> trace

    /**
     * Enable syslog appender
     */
    @Parameter(names = ['-syslog'], description = 'Send logs to syslog server (eg. localhost:514)' )
    String syslog

    /**
     * Print out the version number and exit
     */
    @Parameter(names = ['-v','-version'], description = 'Print the program version')
    boolean version

    /**
     * Print out the 'help' and exit
     */
    @Parameter(names = ['-h'], description = 'Print this help', help = true)
    boolean help

    @Parameter(names = ['-q','-quiet'], description = 'Do not print information messages' )
    boolean quiet

    @Parameter(names = ['-bg'], description = 'Execute nextflow in background', arity = 0)
    boolean background

    @DynamicParameter(names = ['-D'], description = 'Set JVM properties' )
    Map<String,String> jvmOpts = [:]

    @Parameter(names = ['-self-update'], description = 'Update nextflow to the latest version', arity = 0, hidden = true)
    boolean selfUpdate

    @Parameter(names = ['-d','-dockerize'], description = 'Launch nextflow via Docker (experimental)', arity = 0)
    boolean dockerize

    Boolean ansiLog

    boolean getAnsiLog() {
        if( ansiLog && quiet )
            throw new AbortOperationException("Command line options `quiet` and `ansi-log` cannot be used together")

        if( ansiLog != null )
            return ansiLog

        if( background )
            return ansiLog = false

        if( quiet )
            return ansiLog = false

        final env = System.getenv('NXF_ANSI_LOG')
        if( env ) try {
            return Boolean.parseBoolean(env)
        }
        catch (Exception e) {
            log.warn "Invalid boolean value for variable NXF_ANSI_LOG: $env -- it must be 'true' or 'false'"
        }
        return Ansi.isEnabled()
    }

    boolean hasAnsiLogFlag() {
        ansiLog==true || System.getenv('NXF_ANSI_LOG')=='true'
    }

}
