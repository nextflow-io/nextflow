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

package nextflow.cli.v1

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.ILauncherOptions
import nextflow.cli.PluginImpl
import nextflow.exception.AbortOperationException
import static nextflow.cli.PluginExecAware.CMD_SEP

/**
 * CLI `plugin` sub-command (v1)
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Manage plugins and execute custom plugin commands')
class PluginCmd extends AbstractCmd {

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
            PluginImpl.install(args[1].tokenize(','))
        }

        // plugin run command
        else if( args[0].contains(CMD_SEP) ) {
            PluginImpl.exec(args.pop(), args, launcher.options)
        }

        else {
            throw new AbortOperationException("Invalid plugin command: ${args[0]}")
        }
    }

}
