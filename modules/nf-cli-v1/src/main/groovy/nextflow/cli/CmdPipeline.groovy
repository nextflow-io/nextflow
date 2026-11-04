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
import nextflow.cli.pipeline.CmdPipelineInstall
import nextflow.cli.pipeline.CmdPipelineList
import nextflow.cli.pipeline.CmdPipelineRemove
import nextflow.cli.pipeline.CmdPipelineSearch
import nextflow.cli.pipeline.CmdPipelineView
import nextflow.exception.AbortOperationException

/**
 * Implements `pipeline` command
 */
@CompileStatic
@Parameters(commandDescription = "Manage registry pipelines")
class CmdPipeline extends CmdBase implements UsageAware {

    static final public String NAME = 'pipeline'

    private JCommander jCommander

    static final List<CmdBase> commands = new ArrayList<>()

    static {
        commands << new CmdPipelineInstall()
        commands << new CmdPipelineList()
        commands << new CmdPipelineRemove()
        commands << new CmdPipelineSearch()
        commands << new CmdPipelineView()
    }

    protected JCommander commander() {
        if( !this.jCommander ) {
            this.jCommander = new JCommander(this)
            this.jCommander.setProgramName('nextflow pipeline')
            commands.each { cmd ->
                cmd.launcher = this.launcher
                this.jCommander.addCommand(cmd.getName(), cmd, cmd.getAliases() as String[])
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
            jc.parse((args + unknownOptions) as String[])

            final parsedCommand = jc.getParsedCommand()
            if( !parsedCommand ) {
                jc.usage()
                return
            }

            final subcommand = jc.getCommands()
                .get(parsedCommand)
                .getObjects()[0] as CmdBase
            subcommand.run()
        }
        catch( ParameterException e ) {
            throw new AbortOperationException("${e.getMessage()} -- Check the available commands and options and syntax with 'nextflow pipeline -h'")
        }
    }

    private CmdBase findCmd(String name) {
        commands.find { it.name == name || name in it.aliases }
    }

    @Override
    void usage() {
        usage(args)
    }

    @Override
    void usage(List<String> args) {
        def result = []
        if( !args ) {
            result << 'Usage: nextflow pipeline <command> [options]'
            result << ''
            result << 'Commands:'
            commands.each {
                def description = it.getClass().getAnnotation(Parameters)?.commandDescription()
                result << "  ${it.name.padRight(12)}${description}"
            }
            result << ''
            println result.join('\n').toString()
        }
        else {
            final sub = findCmd(args[0])
            if( sub ) {
                commander().usage(args[0])
            }
            else {
                throw new AbortOperationException("Unknown pipeline sub-command: ${args[0]}")
            }
        }
    }
}
