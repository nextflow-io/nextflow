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

package nextflow.cli

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.module.ModuleGenerateMeta
import nextflow.cli.module.ModuleInfo
import nextflow.cli.module.ModuleInstall
import nextflow.cli.module.ModuleList
import nextflow.cli.module.ModulePublish
import nextflow.cli.module.ModuleRemove
import nextflow.cli.module.ModuleRun
import nextflow.cli.module.ModuleSearch
import nextflow.exception.AbortOperationException

/**
 * Implements `module` command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Slf4j
@Parameters(commandDescription = "Manage Nextflow modules")
class CmdModule extends CmdBase implements UsageAware {

    static final public String NAME = 'module'

    private JCommander jCommander

    static final List<CmdBase> commands = new ArrayList<>()

    static {
        commands << new ModuleGenerateMeta()
        commands << new ModuleInstall()
        commands << new ModuleRun()
        commands << new ModuleList()
        commands << new ModuleRemove()
        commands << new ModuleSearch()
        commands << new ModuleInfo()
        commands << new ModulePublish()
    }

    protected JCommander commander() {
        if( !this.jCommander ) {
            this.jCommander = new JCommander(this)
            this.jCommander.setProgramName('nextflow module')
            // Register all subcommands
            commands.each { cmd ->
                cmd.launcher = this.launcher
                this.jCommander.addCommand(cmd.getName(), cmd, new String[0])
            }
        }
        return jCommander
    }

    @Parameter
    List<String> args

    @Override
    String getName() {
        return NAME
    }

    @Override
    void run() {


        try {
            if( !args ) {
                usage()
                return
            }
            final jc = commander()
            final moduleArgs = args + unknownOptions
            jc.parse(moduleArgs as String[])

            final parsedCommand = jc.getParsedCommand()
            if( !parsedCommand ) {
                jc.usage()
                return
            }

            // Get the parsed subcommand instance
            final subcommand = jc.getCommands()
                .get(parsedCommand)
                .getObjects()[0] as CmdBase

            // Execute with fields already populated by JCommander
            subcommand.run()

        } catch( ParameterException e ) {
            throw new AbortOperationException("${e.getMessage()} -- Check the available commands and options and syntax with 'nextflow module -h'")
        }
    }

    private CmdBase findCmd(String name) {
        commands.find { it.name == name }
    }

    /**
     * Print the command usage help
     */
    @Override
    void usage() {
        usage(args)
    }

    /**
     * Print the command usage help
     *
     * @param args The arguments as entered by the user
     */
    @Override
    void usage(List<String> args) {
        def result = []
        if( !args ) {
            result << 'Usage: nextflow module <command> [options]'
            result << ''
            result << 'Commands:'
            commands.each {
                def description = it.getClass().getAnnotation(Parameters)?.commandDescription()
                result << "  ${it.name.padRight(15)}${description}"
            }
            result << ''
            println result.join('\n').toString()
        } else {
            final sub = findCmd(args[0])
            if( sub ) {
                commander().usage(args[0])
            } else {
                throw new AbortOperationException("Unknown module sub-command: ${args[0]}")
            }
        }
    }
}
