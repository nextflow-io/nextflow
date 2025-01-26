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
 *
 */

package nextflow.cli

import com.beust.jcommander.Parameter
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdCid extends CmdBase {

    private static final String NAME = 'cid'

    interface SubCmd {
        String getName()
        void apply(List<String> result)
        void usage(List<String> result)
    }

    private List<SubCmd> commands = new ArrayList<>()

    CmdCid() {

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

        getCmd(args).apply(args.drop(1))
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
}
