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
import java.nio.file.spi.FileSystemProvider

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.script.AssetManager
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

    @Parameter(names='-d',description = 'Show detailed information', arity = 0)
    boolean detailed

    @Override
    final String getName() { 'info' }

    @Override
    void run() {
        if( !args ) {
            println getInfo(detailed) + '\n'
            return
        }

        def repo = new AssetManager().resolveName(args[0])
        if( !repo ) {
            throw new AbortOperationException("Unknown pipeline '${args[0]}'")
        }

        def pipeline = new AssetManager(repo as String)

        println " repo name  : ${pipeline.name}"
        println " home page  : ${pipeline.homePage}"
        println " local path : ${pipeline.localPath}"
        println " main script: ${pipeline.mainScriptName}"

        def revs = pipeline.getRevisions()
        if( revs.size() == 1 )
            println " revision   : ${revs[0]}"
        else {
            println " revisions  : "
            revs.each { println " $it" }
        }
    }

    /**
     * @return A string containing some system runtime information
     */
    static String getInfo(boolean detailed = false) {

        final BLANK = '  '
        final NEWLINE = '\n'
        def result = new StringBuilder()
        result << BLANK << "Version: ${Const.APP_VER} build ${Const.APP_BUILDNUM}" << NEWLINE
        result << BLANK << "Last modified: ${Const.APP_TIMESTAMP_UTC} ${Const.deltaLocal()}" << NEWLINE
        result << BLANK << "Os: ${System.getProperty('os.name')} ${System.getProperty('os.version')}" << NEWLINE
        result << BLANK << "Groovy: ${GroovySystem.getVersion()}" << NEWLINE
        result << BLANK << "Jvm: ${System.getProperty('java.vendor')} ${System.getProperty('java.runtime.version')}" << NEWLINE
        result << BLANK << "Encoding: ${System.getProperty('file.encoding')} (${System.getProperty('sun.jnu.encoding')})" << NEWLINE
        result << BLANK << "Address: ${getLocalNameAndAddress()}" << NEWLINE

        if( !detailed )
            return result.toString()

        List<String> capsule = []
        List<String> args = []
        ManagementFactory
                .getRuntimeMXBean()
                .getInputArguments()
                .each { String it ->
                        if( it.startsWith('-Dcapsule.'))
                            capsule << it.substring(2)
                        else
                            args << it
                    }

        String[] classPath = System.getProperty('java.class.path').split(':')

        // file system
        result << BLANK << "File systems: "
        result << FileSystemProvider.installedProviders().collect { it.scheme }.join(', ')
        result << NEWLINE

        // JVM options
        result << BLANK << "Opts:" << NEWLINE
        for( String x : args ) {
            result << BLANK << BLANK << x << NEWLINE
        }

        // Capsule options
        result << BLANK << "Capsule:" << NEWLINE
        for( String x : capsule ) {
            result << BLANK << BLANK << x << NEWLINE
        }

        // Env
        result << BLANK << "Environment:" << NEWLINE
        for( Map.Entry<String,String> entry : System.getenv().entrySet() ) {
            if( entry.key.startsWith('NXF_') ) {
                result << BLANK << BLANK << entry.key << '=' << entry.value << NEWLINE
            }
        }

        // Class path
        result << BLANK << "Class-path:" << NEWLINE
        for( String x : classPath ) {
            result << BLANK << BLANK << x << NEWLINE
        }

        return result.toString()
    }

    /**
     * @return A string holding the local host name and address used for logging
     */
    static private getLocalNameAndAddress() {
        def host = InetAddress.getLocalHost()
        "${host.getHostName()} [${host.getHostAddress()}]"
    }
}
