/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
/**
 * CLI sub-command HELP
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Print the usage help for a command")
class CmdHelp extends CmdBase {

    static final public NAME = 'help'

    @Override
    final String getName() { NAME }

    @Parameter(description = 'command name', arity = 1)
    List<String> args

    private UsageAware getUsage( List<String> args ) {
        def result = args ? launcher.findCommand(args[0]) : null
        result instanceof UsageAware ? result as UsageAware: null
    }

    @Override
    void run() {
        def cmd = getUsage(args)
        if( cmd ) {
            cmd.usage(args.size()>1 ? args[1..-1] : Collections.<String>emptyList())
        }
        else {
            launcher.usage(args ? args[0] : null)
        }
    }
}
