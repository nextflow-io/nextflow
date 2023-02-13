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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager

/**
 * CLI `view` sub-command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ViewImpl {

    interface Options {
        String getPipeline()
        boolean getQuiet()
        boolean getAll()
    }

    @Delegate
    private Options options

    ViewImpl(Options options) {
        this.options = options
    }

    void run() {
        Plugins.init()
        def manager = new AssetManager(pipeline)
        if( !manager.isLocal() )
            throw new AbortOperationException("Unknown project name `${pipeline}`")

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
