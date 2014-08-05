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

import java.lang.management.ManagementFactory

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.script.PipelineManager

/**
 * CLI sub-command INFO
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Show the system runtime information")
class CmdInfo implements CmdX {

    @Parameter(description = 'pipeline name')
    List<String> args

    @Override
    final String getName() { 'info' }

    @Override
    void run() {
        if( !args ) {
            println getInfo() + '\n'
            return
        }

        def repo = new PipelineManager().resolveName(args[0])
        if( !repo ) {
            throw new AbortOperationException("Unknown pipeline '${args[0]}'")
        }

        def pipeline = new PipelineManager(repo as String)

        println " repo name  : ${pipeline.name}"
        println " home page  : ${pipeline.homePage}"
        println " local path : ${pipeline.localPath}"
        println " main script: ${pipeline.mainScriptName}"
        println " default rev: ${pipeline.masterBranch}"
        println " current rev: ${pipeline.currentRevision}"
    }

    /**
     * @return A string containing some system runtime information
     */
    static String getInfo() {

"""\
      Version: ${Const.APP_VER} build ${Const.APP_BUILDNUM}
      Last modified: ${Const.APP_TIMESTAMP_UTC} ${Const.deltaLocal()}
      Os: ${System.getProperty('os.name')} ${System.getProperty('os.version')}
      Groovy: ${GroovySystem.getVersion()}
      Jvm: ${System.getProperty('java.vendor')} ${System.getProperty('java.runtime.version')}
      Opts: ${ManagementFactory.getRuntimeMXBean().getInputArguments().join(' ')}
      Encoding: ${System.getProperty('file.encoding')} (${System.getProperty('sun.jnu.encoding')})
      Address: ${getLocalNameAndAddress()}"""

    }

    /**
     * @return A string holding the local host name and address used for logging
     */
    static private getLocalNameAndAddress() {
        def host = InetAddress.getLocalHost()
        "${host.getHostName()} [${host.getHostAddress()}]"
    }
}
