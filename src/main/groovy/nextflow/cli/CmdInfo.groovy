/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import java.lang.management.ManagementFactory
import java.nio.file.spi.FileSystemProvider

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import com.sun.management.OperatingSystemMXBean
import groovy.transform.CompileStatic
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.scm.AssetManager
import nextflow.util.MemoryUnit

/**
 * CLI sub-command INFO
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Print project and system runtime information")
class CmdInfo extends CmdBase {

    static final public NAME = 'info'

    @Parameter(description = 'project name')
    List<String> args

    @Parameter(names='-d',description = 'Show detailed information', arity = 0)
    boolean detailed

    @Parameter(names='-dd', hidden = true, arity = 0)
    boolean moreDetailed

    @Override
    final String getName() { NAME }

    @Override
    void run() {
        int level = moreDetailed ? 2 : ( detailed ? 1 : 0 )
        if( !args ) {
            println getInfo(level)
            return
        }

        final manager = new AssetManager(args[0])
        if( !manager.isLocal() )
            throw new AbortOperationException("Unknown project `${args[0]}`")

        final manifest = manager.getManifest()
        println " project name: ${manager.project}"
        println " repository  : ${manager.repositoryUrl}"
        println " local path  : ${manager.localPath}"
        println " main script : ${manager.mainScriptName}"
        if( manager.homePage && manager.homePage != manager.repositoryUrl )
            println " home page   : ${manager.homePage}"
        if( manifest.description )
            println " description : ${manifest.description}"
        if( manifest.author )
            println " author      : ${manifest.author}"

        def revs = manager.getRevisions(level)
        if( revs.size() == 1 )
            println " revision    : ${revs[0]}"
        else {
            println " revisions   : "
            revs.each { println " $it" }
        }

        def updates = manager.getUpdates(level)
        if( updates ) {
            if( updates.size() == 1 && revs.size() == 1 )
                println " updates     : ${updates[0]}"
            else {
                println " updates     : "
                updates.each { println " $it" }
            }
        }

    }

    final static private BLANK = '  '
    final static private NEWLINE = '\n'

    static String status(boolean detailed = false) {
        getInfo( detailed ? 1 : 0, true )
    }

    /**
     * @return A string containing some system runtime information
     */
    static String getInfo(int level, boolean printProc=false) {

        def props = System.getProperties()
        def result = new StringBuilder()
        result << BLANK << "Version: ${Const.APP_VER} build ${Const.APP_BUILDNUM}" << NEWLINE
        result << BLANK << "Modified: ${Const.APP_TIMESTAMP_UTC} ${Const.deltaLocal()}" << NEWLINE
        result << BLANK << "System: ${props['os.name']} ${props['os.version']}" << NEWLINE
        result << BLANK << "Runtime: Groovy ${GroovySystem.getVersion()} on ${System.getProperty('java.vm.name')} ${props['java.runtime.version']}" << NEWLINE
        result << BLANK << "Encoding: ${System.getProperty('file.encoding')} (${System.getProperty('sun.jnu.encoding')})" << NEWLINE

        if( printProc ) {
            def OS = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
            result << BLANK << "Process: ${ManagementFactory.getRuntimeMXBean().getName()} " << getLocalAddress() << NEWLINE
            result << BLANK << "CPUs: ${OS.availableProcessors} - Mem: ${new MemoryUnit(OS.totalPhysicalMemorySize)} (${new MemoryUnit(OS.freePhysicalMemorySize)}) - Swap: ${new MemoryUnit(OS.totalSwapSpaceSize)} (${new MemoryUnit(OS.freeSwapSpaceSize)})"
        }

        if( level == 0  )
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

        // file system
        result << BLANK << "File systems: "
        result << FileSystemProvider.installedProviders().collect { it.scheme }.join(', ')
        result << NEWLINE

        // JVM options
        result << BLANK << "JVM opts:" << NEWLINE
        for( String entry : args ) {
            int p = entry.indexOf('=')
            if( p!=-1 ) {
                def key = entry.substring(0,p)
                def value = entry.substring(p+1)
                dump(key,value,2, result)
            }
            else {
                dump(entry,null,2, result)
            }
        }

        // Capsule options
        dump("Capsule", capsule, 1, result)

        // Env
        result << BLANK << "Environment:" << NEWLINE
        def entries = System.getenv().keySet().sort()
        for( def key : entries ) {
            if( key.startsWith('NXF_') || level>1 ) {
                dump(  key, System.getenv().get(key), 2, result )
            }
        }

        // java properties
        if( level>1 ) {
            result << BLANK << "Properties:" << NEWLINE
            entries = System.getProperties().keySet().sort()
            for( String key : entries ) {
                if( key == 'java.class.path' ) continue
                dump( key, System.getProperties().get(key), 2, result )
            }
        }

        // Class path
        dump("Class-path" , System.getProperty('java.class.path'), 1, result)

        // final string
        return result.toString()
    }

    static private void dump(String key, def value, int indent, StringBuilder result) {

        if( value instanceof String && value.count(':') && (key=~/.*PATH$/ || key=='PERL5LIB' || key.contains('.path') || key.contains('-path') || key.contains('.dir')) ) {
            value = value.split(":") as Collection
        }
        else if( value?.class?.isArray() ) {
            value = value as Collection
        }

        if ( value instanceof Collection ) {
            blanks(indent, result)
            result << key << (indent==1 ? ':' : '=')
            for( String item : value ) {
                result << NEWLINE
                blanks(indent+1, result)
                result << item
            }
        }
        else if( value ) {
            blanks(indent, result)
            result << key << (indent==1 ? ':' : '=')
            if( value == '\n' ) {
                result << '\\n'
            }
            else if( value == '\r' ) {
                result << '\\r'
            }
            else if( value == '\n\r' ) {
                result << '\\n\\r'
            }
            else if( value == '\r\n' ) {
                result << '\\r\\n'
            }
            else {
                result << value
            }
        }
        else if( key ) {
            blanks(indent, result)
            result << key
        }
        result << NEWLINE
    }

    static private blanks( int n, StringBuilder result ) {
        for( int i=0; i<n; i++ ) {
            result << BLANK
        }
    }

    /**
     * @return A string holding the local host name and address used for logging
     */
    static private getLocalAddress() {
        try {
            return "[${InetAddress.getLocalHost().getHostAddress()}]"
        }
        catch(Exception e) {
            return null
        }
    }
}
