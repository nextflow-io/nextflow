/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
import com.beust.jcommander.Parameter
/**
 * Main application command line options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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

}
