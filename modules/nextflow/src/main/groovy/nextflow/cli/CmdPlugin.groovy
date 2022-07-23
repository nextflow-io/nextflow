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
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins

/**
 * Plugin manager command
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdPlugin extends CmdBase {

    @Override
    String getName() {
        return 'plugin'
    }

    @Parameter(hidden = true)
    List<String> args

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("Missing plugin command - usage: nextflow plugin install <pluginId,..>")
        // setup plugins system
        Plugins.setup()
        // check for the plugins install
        if( args[0] == 'install' ) {
            if( args.size()!=2 )
                throw new AbortOperationException("Missing plugin install target - usage: nextflow plugin install <pluginId,..>")
            Plugins.pull(args[1].tokenize(','))
        }
        // plugin run command
        else if( args[0] == 'exec' ) {
            if( args.size()<2 )
                throw new AbortOperationException("Missing plugin run target - usage: nextflow plugin exec <plugin name>")
            final target = args[1]
            Plugins.start(target)
            final wrapper = Plugins.manager.getPlugin(target)
            if( !wrapper )
                throw new AbortOperationException("Cannot find target plugin: $target")
            final plugin = wrapper.getPlugin()
            if( plugin instanceof PluginExecAware ) {
                final ret = plugin.exec(getLauncher(), args.subList(2, args.size()))
                // use explicit exit to invoke the system shutdown hooks
                System.exit(ret)
            }
            else
                throw new AbortOperationException("Invalid target plugin: $target")
        }
        else {
            throw new AbortOperationException("Invalid plugin command")
        }
    }

}
