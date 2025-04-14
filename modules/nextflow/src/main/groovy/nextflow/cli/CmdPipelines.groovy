/*
 * Copyright 2013-2024, Seqera Labs
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
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.util.TestOnly

/**
 * CLI sub-command PIPELINES. Manage locally installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Manage pipeline projects")
class CmdPipelines extends CmdBase {

    private static final List<SubCmd> commands = [
        new CmdDrop(),
        new CmdList(),
        new CmdView()
    ]

    static interface SubCmd {
        String getName()
        String getDescription()
        void apply(String pipeline)
        String usage()
    }

    static class CmdDrop implements SubCmd {

        @Override
        String getName() { 'drop' }

        @Override
        String getDescription() { 'Delete the local copy of a project' }

        @Override
        void apply(String pipeline) {
            if( !pipeline )
                throw new AbortOperationException("No project name was specified")

            Plugins.init()

            final manager = new AssetManager(pipeline)
            if( !manager.localPath.exists() )
                throw new AbortOperationException("No match found for: ${pipeline}")

            if( !manager.isClean() )
                throw new AbortOperationException("Local project repository contains uncommitted changes -- won't drop it")

            manager.close()
            if( !manager.localPath.deleteDir() )
                throw new AbortOperationException("Unable to delete project `${manager.project}` -- Check access permissions for path: ${manager.localPath}")
        }

        @Override
        String usage() {
            "Usage: nextflow pipelines drop <pipeline>"
        }
    }

    static class CmdList implements SubCmd {

        @Override
        String getName() { 'list' }

        @Override
        String getDescription() { 'List all downloaded projects' }

        @Override
        void apply(String nope) {
            final all = AssetManager.list()
            if( !all ) {
                log.info '(none)'
                return
            }

            all.each { println it }
        }

        @Override
        String usage() {
            "Usage: nextflow pipelines list"
        }
    }

    static class CmdView implements SubCmd {

        @Override
        String getName() { 'view' }

        @Override
        String getDescription() { 'View project script file(s)' }

        @Override
        void apply(String pipeline) {
            if( !pipeline )
                throw new AbortOperationException("No project name was specified")

            Plugins.init()
            final manager = new AssetManager(pipeline)
            if( !manager.isLocal() )
                throw new AbortOperationException("Unknown project name `${pipeline}`")

            // list repository content
            println "== content of pipeline directory: ${manager.localPath}"
            manager.localPath.eachFile { File it ->
                println it.name
            }

            // print the main script
            final script = manager.getMainScriptFile()
            if( !script.exists() )
                throw new AbortOperationException("Missing script file: '${script}'")

            println()
            println "== content of main script: $script"
            script.readLines().each { println it }
        }

        @Override
        String usage() {
            "Usage: nextflow pipelines view <pipeline>"
        }
    }

    @Parameter
    List<String> args

    @Override
    final String getName() { 'pipelines' }

    @Override
    void run() {
        if( !args ) {
            usage()
            return
        }

        final cmd = findCmd(args.pop())
        if( !cmd ) {
            throw new AbortOperationException("Unknown pipelines sub-command: `$cmd`")
        }

        final pipeline = args ? args.first() : null
        cmd.apply(pipeline)
    }

    void usage() {
        final result = []
        result << 'Usage: nextflow pipelines <command>'
        result << ''
        result << 'Commands:'
        commands.each {
            result << "  ${it.name}\t${it.description}"
        }
        result << ''
        println result.join('\n').toString()
    }

    private SubCmd findCmd(String name) {
        commands.find { it.name == name }
    }

}
