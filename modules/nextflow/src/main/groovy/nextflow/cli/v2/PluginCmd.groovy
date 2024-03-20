/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cli.v2

import groovy.transform.CompileStatic
import nextflow.cli.CliOptions
import nextflow.cli.CmdPlugin
import nextflow.exception.AbortOperationException
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import picocli.CommandLine.ParentCommand
import static nextflow.cli.PluginExecAware.CMD_SEP

/**
 * CLI `plugin` sub-command (v2)
 * 
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'plugin',
    description = 'Manage plugins and execute custom plugin commands'
)
class PluginCmd extends AbstractCmd {

    @ParentCommand
    private Launcher launcher

    @Parameters(index = '0')
    String command

    @Parameters(index = '1..*')
    List<String> args = []

    @Override
    void run() {
        // plugin install command
        if( command == 'install' ) {
            if( args.size()!=1 )
                throw new AbortOperationException("Missing plugin install target - usage: nextflow plugin install <pluginId,..>")
            CmdPlugin.install(args[0].tokenize(','))
        }

        // plugin run command
        else if( command.contains(CMD_SEP) ) {
            CmdPlugin.exec(command, args, launcher.options)
        }

        else {
            throw new AbortOperationException("Invalid plugin command: ${command}")
        }
    }

}
