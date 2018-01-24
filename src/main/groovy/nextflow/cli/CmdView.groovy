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

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.scm.AssetManager

/**
 * CLI sub-command VIEW -- Print a pipeline script to console
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "View project script file(s)")
class CmdView extends CmdBase {

    static final public NAME = 'view'

    @Override
    String getName() { NAME }

    @Parameter(description = 'project name', required = true)
    List<String> args = []

    @Parameter(names = '-q', description = 'Hide header line', arity = 0)
    boolean quiet

    @Parameter(names = '-l', description = 'List repository content', arity = 0)
    boolean all

    @Override
    void run() {

        def manager = new AssetManager(args[0])
        if( !manager.isLocal() )
            throw new AbortOperationException("Unknown project name `${args[0]}`")

        if( all ) {
            if( !quiet )
                println "== content of path: ${manager.localPath}"

            manager.localPath.eachFile { File it ->
                println it.name
            }
        }

        else {
            /*
             * prints the script main file
             */
            final script = manager.getMainScriptFile()
            if( !script.exists() )
                throw new AbortOperationException("Missing script file: '${script}'")

            if( !quiet )
                println "== content of file: $script"

            script.readLines().each { println it }
        }

    }
}
