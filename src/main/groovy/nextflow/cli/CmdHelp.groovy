/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.cli
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import picocli.CommandLine

/**
 * CLI sub-command HELP
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Print the usage help for a command")
@CommandLine.Command
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
