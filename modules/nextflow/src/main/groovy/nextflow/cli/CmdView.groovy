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
