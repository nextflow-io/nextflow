/*
 * Copyright 2020-2022, Seqera Labs
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

import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import org.fusesource.jansi.Ansi
import picocli.CommandLine.Option

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
    @Option(hidden = true, names = ['-debug'])
    List<String> debug

    @Option(names = ['-log'], description = 'Set nextflow log file path')
    String logFile

    @Option(names = ['-c','-config'], description = 'Add the specified file to configuration set')
    List<String> userConfig

    @Option(names = ['-config-ignore-includes'], description = 'Disable the parsing of config includes')
    boolean ignoreConfigIncludes

    @Option(names = ['-C'], description = 'Use the specified configuration file(s) overriding any defaults')
    List<String> config

    /**
     * the packages to trace
     */
    @Option(names = ['-trace'], description = 'Enable trace level logging for the specified package name - multiple packages can be provided separating them with a comma e.g. \'-trace nextflow,io.seqera\'')
    List<String> trace

    /**
     * Enable syslog appender
     */
    @Option(names = ['-syslog'], arity = '0..1', fallbackValue = 'localhost', description = 'Send logs to syslog server (eg. localhost:514)' )
    String syslog

    /**
     * Print out the version number and exit
     */
    @Option(names = ['-v'], description = 'Print the program version')
    boolean version

    @Option(names = ['-version'], description = 'Print the program version (full)')
    boolean fullVersion

    @Option(names = ['-q','-quiet'], description = 'Do not print information messages' )
    boolean quiet

    @Option(names = ['-bg'], arity = '0', description = 'Execute nextflow in background')
    boolean background

    @Option(names = ['-D'], description = 'Set JVM properties' )
    Map<String,String> jvmOpts

    @Option(names = ['-self-update'], arity = '0', description = 'Update nextflow to the latest version', hidden = true)
    boolean selfUpdate

    @Option(names = ['-d','-dockerize'], arity = '0', description = 'Launch nextflow via Docker (experimental)')
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
