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

import static nextflow.file.FileHelper.toCanonicalPath

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.attribute.BasicFileAttributes

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.config.ConfigBuilder
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
@Slf4j
class CmdFs {

    static final public NAME = 'fs'

    enum Command {
        COPY,
        MOVE,
        LIST,
        CAT,
        REMOVE,
        STAT
    }

    interface Options {
        CliOptions getLauncherOptions()
    }

    @Parameters(commandDescription = 'Perform basic filesystem operations')
    static class V1 extends CmdBase implements UsageAware, Options {

        trait SubCmd {
            abstract int getArity()

            abstract String getName()

            abstract String getDescription()

            abstract Command getCommand()

            String usage() {
                "Usage: nextflow fs ${name} " + (arity==1 ? "<path>" : "source_file target_file")
            }
        }

        static class CmdCopy implements SubCmd {
            @Override
            int getArity() { 2 }

            @Override
            String getName() { 'cp' }

            @Override
            String getDescription() { 'Copy a file' }

            @Override
            Command getCommand() { Command.COPY }
        }

        static class CmdMove implements SubCmd {
            @Override
            int getArity() { 2 }

            @Override
            String getName() { 'mv' }

            @Override
            String getDescription() { 'Move a file' }

            @Override
            Command getCommand() { Command.MOVE }
        }

        static class CmdList implements SubCmd {
            @Override
            int getArity() { 1 }

            @Override
            String getDescription() { 'List the content of a folder' }

            @Override
            String getName() { 'ls' }

            @Override
            Command getCommand() { Command.LIST }
        }

        static class CmdCat implements SubCmd {
            @Override
            int getArity() { 1 }

            @Override
            String getName() { 'cat' }

            @Override
            String getDescription() { 'Print a file to the stdout' }

            @Override
            Command getCommand() { Command.CAT }
        }

        static class CmdRemove implements SubCmd {
            @Override
            int getArity() { 1 }

            @Override
            String getName() { 'rm' }

            @Override
            String getDescription() { 'Remove a file' }

            @Override
            Command getCommand() { Command.REMOVE }
        }

        static class CmdStat implements SubCmd {
            @Override
            int getArity() { 1 }

            @Override
            String getName() { 'stat' }

            @Override
            String getDescription() { 'Print file metadata' }

            @Override
            Command getCommand() { Command.STAT }
        }

        private List<SubCmd> commands = (List<SubCmd>)[
            new CmdCopy(),
            new CmdMove(),
            new CmdList(),
            new CmdCat(),
            new CmdRemove(),
            new CmdStat()
        ]

        @Parameter
        List<String> args = []

        @Override
        CliOptions getLauncherOptions() {
            launcher.options
        }

        @Override
        String getName() {
            return NAME
        }

        @Override
        void run() {
            if( !args ) {
                usage()
                return
            }

            final cmd = findCmd(args[0])
            if( !cmd ) {
                throw new AbortOperationException("Unknown file system command: `$cmd`")
            }

            if( args.size() - 1 != cmd.getArity() )
                throw new AbortOperationException(cmd.usage())

            new CmdFs(this).run(cmd.getCommand(), args.drop(1))
        }

        private SubCmd findCmd( String name ) {
            commands.find { it.name == name }
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
                def sub = findCmd(args[0])
                if( sub )
                    println sub.usage()
                else {
                    throw new AbortOperationException("Unknown fs sub-command: ${args[0]}")
                }
            }
        }

    }

    @Delegate
    private Options options

    CmdFs(Options options) {
        this.options = options
    }

    private Session createSession() {
        // create the config
        final config = new ConfigBuilder()
                .setOptions(getLauncherOptions())
                .setBaseDir(Paths.get('.'))
                .build()

        return new Session(config)
    }

    void run(Command command, List<String> args) {
        Plugins.setup()
        final session = createSession()
        try {
            run0(command, args)
        }
        finally {
            try {
                session.destroy()
                Global.cleanUp()
                Plugins.stop()
            } catch (Throwable t) {
                log.warn "Unexpected error while destroying the session object - cause: ${t.message}"
            }
        }
    }

    private void run0(Command command, List<String> args) {
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
            case STAT:
                traverse(args[0]) { Path path -> stat(path) }
                break
        }
    }

    private void traverse( String source, Closure op ) {

        // if it isn't a glob pattern simply return it a normalized absolute Path object
        def splitter = FilePatternSplitter.glob().parse(source)
        if( splitter.isPattern() ) {
            final scheme = splitter.scheme
            final target = scheme ? "$scheme://$splitter.parent" : splitter.parent
            final folder = toCanonicalPath(target)
            final pattern = splitter.fileName

            def opts = [:]
            opts.type = 'any'

            FileHelper.visitFiles(opts, folder, pattern, op)
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

    void stat(Path source) {
        try {
            final attr = Files.readAttributes(source, BasicFileAttributes)
            print """\
                name          : ${source.name}
                size          : ${attr.size()}
                is directory  : ${attr.isDirectory()}
                last modified : ${attr.lastModifiedTime() ?: '-'}
                creation time : ${attr.creationTime() ?: '-'}                    
                """.stripIndent()
        }
        catch (IOException e) {
            log.warn "Unable to read attributes for file: ${source.toUriString()} - cause: $e.message", e
        }
    }

}
