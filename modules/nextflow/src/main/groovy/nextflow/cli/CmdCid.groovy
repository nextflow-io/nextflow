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
import nextflow.config.ConfigBuilder
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Explore workflows CID metadata")
class CmdCid extends CmdBase implements UsageAware {

    private static final String NAME = 'cid'

    interface CidOperation {
        void log(ConfigMap config)
        void show(ConfigMap config, List<String> args)
        void lineage(ConfigMap config, List<String> args)
        void diff(ConfigMap config, List<String> args)
    }

    interface SubCmd {
        String getName()
        String getDescription()
        void apply(List<String> args)
        void usage()
    }

    private List<SubCmd> commands = new ArrayList<>()

    private CidOperation operations

    private ConfigMap config

    CmdCid() {
        commands << new CmdLog()
        commands << new CmdShow()
        commands << new CmdLineage()
        commands << new CmdDiff()
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
            return
        }
        // setup the plugins system and load the secrets provider
        Plugins.init()
        // load the config
        this.config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(Paths.get('.'))
            .build()
        // load the command operations
        this.operations = ServiceLoader.load(CidOperation.class).findFirst().orElse(null)
        if( !operations )
            throw new IllegalStateException("Unable to load CID plugin")
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

    class CmdLog implements SubCmd {

        @Override
        String getName() {
            return 'log'
        }

        @Override
        String getDescription() {
            return 'Print the CID execution log'
        }

        @Override
        void apply(List<String> args) {
            if (args.size() != 0) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }
            operations.log(config)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name"
        }
    }

    class CmdShow implements SubCmd{

        @Override
        String getName() {
            return 'show'
        }

        @Override
        String getDescription() {
            return 'Print the description of a CID reference'
        }

        void apply(List<String> args) {
            if (args.size() != 1) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }

            operations.show(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <CID reference> "
        }
    }

    class CmdLineage implements SubCmd {

        @Override
        String getName() { 'lineage' }

        @Override
        String getDescription() {
            return 'Render a lineage graph for a workflow output'
        }

        void apply(List<String> args) {
            if (args.size() != 2) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }

            operations.lineage(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <workflow output CID> <html output file>"
        }

    }

    class CmdDiff implements SubCmd {

        @Override
        String getName() { 'diff' }

        @Override
        String getDescription() {
            return 'Show differences between two CID descriptions'
        }

        void apply(List<String> args) {
            if (args.size() != 2) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }
            operations.diff(config, args)
        }

        @Override
        void usage() {
            println description
            println "Usage: nextflow $NAME $name <CID 1> <CID 2>"
        }

    }
}
