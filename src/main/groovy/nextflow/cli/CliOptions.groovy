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

    @Parameter(names=['-log'], description = 'Define the application log file')
    String logFile

    @Parameter(names=['-lib'], description = 'Library extension path')
    String libPath

    /**
     * The packages to trace
     */
    @Parameter(names='-trace', hidden = true)
    List<String> trace


    /**
     * Print out the version number and exit
     */
    @Parameter(names = ['-v','-version'], description = 'Show the program version')
    boolean version

    @Deprecated
    @Parameter(names = ['-info'], description = 'The runtime info', hidden = true)
    boolean info

    /**
     * Print out the 'help' and exit
     */
    @Parameter(names = ['-h','-help'], description = 'Show this help', help = true)
    boolean help

    @Parameter(names = ['-q','-quiet'], description = 'Do not print information messages' )
    boolean quiet

    @Parameter(names = ['-d'], description = 'Start in cluster daemon mode', arity = 0)
    boolean daemon

    @DynamicParameter(names = ['-daemon.'], description = 'Starts in daemon mode and provides extra options', hidden = true )
    Map<String,String> daemonOptions = [:]


    boolean isDaemon() {
        return daemon || daemonOptions.size() > 0
    }



}
