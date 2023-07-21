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
 */

package nextflow.cli
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import nextflow.plugin.Plugins

/**
 * Implements `fs` command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdFs {

    enum Command {
        COPY,
        MOVE,
        LIST,
        CAT,
        REMOVE
    }

    @Parameters(commandDescription = 'Perform basic filesystem operations')
    static class V1 extends CmdBase implements UsageAware {

        trait SubCmd {
            abstract String getName()
            abstract int getArity()
            abstract String getDescription()
            abstract Command getCommand()

            String usage() {
                "Usage: nextflow fs ${name} " + (arity==1 ? "<path>" : "source_file target_file")
            }
        }

        private List<SubCmd> commands = (List<SubCmd>)[
            new CmdCopy(),
            new CmdMove(),
            new CmdList(),
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

            new CmdFs().run(cmd.getCommand(), args.drop(1))
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

        static class CmdCopy implements SubCmd {
            @Override String getName() { 'cp' }
            @Override int getArity() { 2 }
            @Override String getDescription() { 'Copy a file' }
            @Override Command getCommand() { Command.COPY }
        }

        static class CmdMove implements SubCmd {
            @Override String getName() { 'mv' }
            @Override int getArity() { 2 }
            @Override String getDescription() { 'Move a file' }
            @Override Command getCommand() { Command.MOVE }
        }

        static class CmdList implements SubCmd {
            @Override String getName() { 'ls' }
            @Override int getArity() { 1 }
            @Override String getDescription() { 'List the contents of a folder' }
            @Override Command getCommand() { Command.LIST }
        }

        static class CmdCat implements SubCmd {
            @Override String getName() { 'cat' }
            @Override int getArity() { 1 }
            @Override String getDescription() { 'Print a file to stdout' }
            @Override Command getCommand() { Command.CAT }
        }

        static class CmdRemove implements SubCmd {
            @Override String getName() { 'rm' }
            @Override int getArity() { 1 }
            @Override String getDescription() { 'Remove a file' }
            @Override Command getCommand() { Command.REMOVE }
        }

    }

    void run(Command command, List<String> args) {
        Plugins.setup()

        try {
            switch( command ) {
                case COPY:
                    traverse(args[0]) { Path path -> copy(path, args[1] as Path) }
                    break
                case MOVE:
                    traverse(args[0]) { Path path -> move(path, args[1] as Path) }
                    break
                case LIST:
                    traverse(args[0]) { Path path -> list(path) }
                    break
                case CAT:
                    traverse(args[0]) { Path path -> cat(path) }
                    break
                case REMOVE:
                    traverse(args[0]) { Path path -> remove(path) }
                    break
            }
        }
        finally {
            Plugins.stop()
        }
    }

    private void traverse( String source, Closure op ) {
        // if it isn't a glob pattern simply return it a normalized absolute Path object
        def splitter = FilePatternSplitter.glob().parse(source)
        if( splitter.isPattern() ) {
            final scheme = splitter.scheme
            final folder = splitter.parent
            final pattern = splitter.fileName
            final fs = FileHelper.fileSystemForScheme(scheme)

            def opts = [:]
            opts.type = 'file'

            FileHelper.visitFiles(opts, fs.getPath(folder), pattern, op)
        }
        else {
            def normalised = splitter.strip(source)
            op.call(FileHelper.asPath(normalised))
        }
    }

    void copy(Path source, Path target) {
        FilesEx.copyTo(source, target)
    }

    void move(Path source, Path target) {
        FilesEx.moveTo(source, target)
    }

    void list(Path source) {
        println source.name
    }

    void cat(Path source) {
        String line
        def reader = Files.newBufferedReader(source, Charset.defaultCharset())
        while( line = reader.readLine() )
            println line
    }

    void remove(Path source) {
        Files.isDirectory(source) ? FilesEx.deleteDir(source) : FilesEx.delete(source)
    }

}
