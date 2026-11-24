/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli.pipeline

import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.pipeline.PipelineManager
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.util.TestOnly

/**
 * Pipeline remove subcommand
 */
@CompileStatic
@Parameters(commandDescription = "Remove an installed pipeline")
class CmdPipelineRemove extends CmdBase {

    @Parameter(description = "<pipeline>", required = true)
    List<String> args

    @Parameter(names = ["-force"], description = "Remove even with local modifications", arity = 0)
    boolean force = false

    @Parameter(names = ['-g', '-global'], description = 'Remove the pipeline from the user-level store used by `nextflow run`', arity = 0)
    boolean global

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'remove'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Unexpected number of parameters")
        PipelineCmdSupport.checkPipelineName(args[0], null)
        if( global ) {
            removeUserLevel(PipelineCmdSupport.userLevelName(this, args[0]))
            return
        }
        final manager = new PipelineManager(root ?: Path.of('.').toAbsolutePath().normalize(), null)
        if( manager.remove(args[0], force) ) {
            println "Pipeline ${args[0]} removed"
            return
        }
        println "Pipeline ${args[0]} is not installed"
        if( isInUserLevelStore(args[0]) )
            println "Pipeline ${args[0]} is installed in the user-level store -- use `nextflow pipeline remove -g ${args[0]}` to remove it"
    }

    private static boolean isInUserLevelStore(String name) {
        return AssetManager.list().contains(name)
    }

    private void removeUserLevel(String name) {
        Plugins.init()
        try (final manager = new AssetManager(name)) {
            if( manager.isNotInitialized() ) {
                println "Pipeline ${name} is not installed in the user-level store"
                return
            }
            manager.drop(null, force)
            println "Pipeline ${manager.project} removed"
        }
    }
}
