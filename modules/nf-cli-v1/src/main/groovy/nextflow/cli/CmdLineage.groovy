/*
 * Copyright 2013-2025, Seqera Labs
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
 *
 */

package nextflow.cli

import java.nio.file.Paths

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.config.ConfigCmdAdapter
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * CID command line interface
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Explore workflows lineage metadata", commandNames = ['li'])
class CmdLineage extends CmdBase implements UsageAware {

    private static final String NAME = 'lineage'

    interface LinCommand extends ExtensionPoint {
        void list(ConfigMap config)
        void view(ConfigMap config, List<String> args)
        void render(ConfigMap config, List<String> args)
        void diff(ConfigMap config, List<String> args)
        void find(ConfigMap config, List<String> args)
    }

    interface SubCmd {
        String getName()
        String getDescription()
        void apply(List<String> args)
        void usage()
    }

    private List<SubCmd> commands = new ArrayList<>()

    private LinCommand operation

    private ConfigMap config

    CmdLineage() {
        commands << new CmdList()
        commands << new CmdView()
        commands << new CmdRender()
        commands << new CmdDiff()
        commands << new CmdFind()
    }

    @Parameter(hidden = true)
    List<String> args

    @Override
    String getName() {
        return NAME
    }

    @Override
    void run() {
        if( !args ) {
            usage(List.of())
            return
        }
        // setup the plugins system and load the secrets provider
        Plugins.init()
        // load the config
        this.config = new ConfigCmdAdapter()
            .setOptions(launcher.options)
            .setBaseDir(Paths.get('.'))
            .build()
        // init plugins
        Plugins.load(config)
        // load the command operations
        this.operation = Plugins.getExtension(LinCommand)
        if( !operation )
            throw new IllegalStateException("Unable to load lineage extensions.")
        // consume the first argument
        getCmd(args).apply(args.drop(1))
    }

    /**
     * Print the command usage help
     */
    void usage() {
        usage(args)
    }

    /**
     * Print the command usage help
     *
     * @param args The arguments as entered by the user
     */
    void usage(List<String> args) {
        if( !args ) {
            List<String> result = []
            result << this.getClass().getAnnotation(Parameters).commandDescription()
            result << "Usage: nextflow $NAME <sub-command> [options]".toString()
            result << ''
            result << 'Commands:'
            int len = 0
            commands.forEach {len = it.name.size() > len ? it.name.size() : len }
            commands.sort(){it.name}.each { result << "  ${it.name.padRight(len)}\t${it.description}".toString()  }
            result << ''
            println result.join('\n').toString()
        }
        else {
            def sub = commands.find { it.name == args[0] }
            if( sub )
                sub.usage()
            else {
                throw new AbortOperationException("Unknown $NAME sub-command: ${args[0]}")
            }
        }
    }

    protected SubCmd getCmd(List<String> args) {

        def cmd = commands.find { it.name == args[0] }
        if( cmd ) {
            return cmd
        }

        def matches = commands.collect{ it.name }.closest(args[0])
        def msg = "Unknown cloud sub-command: ${args[0]}"
        if( matches )
            msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
        throw new AbortOperationException(msg)
    }

    class CmdList implements SubCmd {

        @Override
        String getName() {
            return 'list'
        }

        @Override
        String getDescription() {
            return 'List the executions with lineage enabled'
        }

        @Override
        void apply(List<String> args) {
            if (args.size() != 0) {
                println("ERROR: Incorrect number of parameters")
                return
            }
            operation.list(config)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name"
        }
    }

    class CmdView implements SubCmd{

        @Override
        String getName() {
            return 'view'
        }

        @Override
        String getDescription() {
            return 'Print the description of a Lineage ID (lid)'
        }

        void apply(List<String> args) {
            if (args.size() != 1) {
                println("ERROR: Incorrect number of parameters")
                return
            }

            operation.view(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <lid> "
        }
    }

    class CmdRender implements SubCmd {

        @Override
        String getName() { 'render' }

        @Override
        String getDescription() {
            return 'Render the lineage graph for a workflow output'
        }

        void apply(List<String> args) {
            if (args.size() < 1 || args.size() > 2) {
                println("ERROR: Incorrect number of parameters")
                return
            }

            operation.render(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <workflow output lid> [<html output file>]"
        }

    }

    class CmdDiff implements SubCmd {

        @Override
        String getName() { 'diff' }

        @Override
        String getDescription() {
            return 'Show differences between two lineage descriptions'
        }

        void apply(List<String> args) {
            if (args.size() != 2) {
                println("ERROR: Incorrect number of parameters")
                return
            }
            operation.diff(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <lid-1> <lid-2>"
        }

    }

    class CmdFind implements SubCmd {

        @Override
        String getName() { 'find' }

        @Override
        String getDescription() {
            return 'Find lineage metadata descriptions matching with a query'
        }

        void apply(List<String> args) {
            if (args.size() < 1) {
                println("ERROR: Incorrect number of parameters")
                return
            }
            operation.find(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <query>"
        }

    }

}
