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
 */

package nextflow.cli.v1

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.FsImpl
import nextflow.exception.AbortOperationException

/**
 * CLI `fs` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Perform basic filesystem operations')
class FsCmd extends AbstractCmd implements UsageAware {

    trait SubCmd {
        abstract String getName()
        abstract int getArity()
        abstract String getDescription()
        abstract FsImpl.Command getCommand()

        String usage() {
            "Usage: nextflow fs ${name} " + (arity==1 ? "<path>" : "source_file target_file")
        }
    }

    private List<SubCmd> commands = (List<SubCmd>)[
        new CmdCopy(),
        new CmdMove(),
        new ListImpl(),
        new CmdCat(),
        new CmdRemove()
    ]

    @Parameter(hidden = true)
    List<String> args = []

    @Override
    String getName() { 'fs' }

    @Override
    void run() {
        if( !args ) {
            usage()
            return
        }

        final cmd = commands.find { it.name == args[0] }
        if( !cmd ) {
            def matches = commands.collect{ it.name }.closest(args[0])
            def msg = "Unknown fs sub-command: ${args[0]}"
            if( matches )
                msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
            throw new AbortOperationException(msg)
        }

        if( args.size() - 1 != cmd.getArity() )
            throw new AbortOperationException(cmd.usage())

        new FsImpl().run(cmd.getCommand(), args.drop(1))
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

        def result = []
        if( !args ) {
            result << 'Usage: nextflow fs <command> [args]'
            result << ''
            result << 'Commands:'
            commands.each {
                result << "  ${it.name}\t${it.description}"
            }
            result << ''
            println result.join('\n').toString()
        }
        else {
            def sub = commands.find { it.name == args[0] }
            if( sub )
                println sub.usage()
            else {
                throw new AbortOperationException("Unknown fs sub-command: ${args[0]}")
            }
        }

    }

    class CmdCopy implements SubCmd {
        @Override String getName() { 'cp' }
        @Override int getArity() { 2 }
        @Override String getDescription() { 'Copy a file' }
        @Override FsImpl.Command getCommand() { FsImpl.Command.COPY }
    }

    class CmdMove implements SubCmd {
        @Override String getName() { 'mv' }
        @Override int getArity() { 2 }
        @Override String getDescription() { 'Move a file' }
        @Override FsImpl.Command getCommand() { FsImpl.Command.MOVE }
    }

    class ListImpl implements SubCmd {
        @Override String getName() { 'ls' }
        @Override int getArity() { 1 }
        @Override String getDescription() { 'List the contents of a folder' }
        @Override FsImpl.Command getCommand() { FsImpl.Command.LIST }
    }

    class CmdCat implements SubCmd {
        @Override String getName() { 'cat' }
        @Override int getArity() { 1 }
        @Override String getDescription() { 'Print a file to stdout' }
        @Override FsImpl.Command getCommand() { FsImpl.Command.CAT }
    }

    class CmdRemove implements SubCmd {
        @Override String getName() { 'rm' }
        @Override int getArity() { 1 }
        @Override String getDescription() { 'Remove a file' }
        @Override FsImpl.Command getCommand() { FsImpl.Command.REMOVE }
    }

}
