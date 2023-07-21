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
 *
 */

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import static nextflow.cli.PluginExecAware.CMD_SEP

/**
 * CLI `plugin` sub-command
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdPlugin {

    @Parameters(commandDescription = 'Manage plugins and execute custom plugin commands')
    static class V1 extends CmdBase {

        @Parameter(hidden = true)
        List<String> args = []

        @Override
        String getName() { 'plugin' }

        @Override
        void run() {
            if( !args )
                throw new AbortOperationException("Missing plugin command - usage: nextflow plugin install <pluginId,..>")

            // plugin install command
            if( args[0] == 'install' ) {
                if( args.size()!=2 )
                    throw new AbortOperationException("Missing plugin install target - usage: nextflow plugin install <pluginId,..>")
                CmdPlugin.install(args[1].tokenize(','))
            }

            // plugin run command
            else if( args[0].contains(CMD_SEP) ) {
                CmdPlugin.exec(args.pop(), args, launcher.options)
            }

            else {
                throw new AbortOperationException("Invalid plugin command: ${args[0]}")
            }
        }

    }

    static void install(List<String> ids) {
        Plugins.setup()
        Plugins.pull(ids)
    }

    static void exec(String command, List<String> args, ILauncherOptions launcherOptions) {
        Plugins.setup()

        final items = command.tokenize(CMD_SEP)
        final target = items[0]
        final targetCmd = items[1] ? items[1..-1].join(CMD_SEP) : null

        // push back the command as the first item
        Plugins.start(target)
        final wrapper = Plugins.manager.getPlugin(target)
        if( !wrapper )
            throw new AbortOperationException("Cannot find target plugin: $target")
        final plugin = wrapper.getPlugin()
        if( plugin instanceof PluginExecAware ) {
            final ret = plugin.exec(launcherOptions, target, targetCmd, args)
            // use explicit exit to invoke the system shutdown hooks
            System.exit(ret)
        }
        else
            throw new AbortOperationException("Invalid target plugin: $target")
    }

}
