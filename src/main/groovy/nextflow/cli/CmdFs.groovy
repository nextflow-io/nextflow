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
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.Parameter
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import picocli.CommandLine

/**
 * Implements `fs` command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CommandLine.Command(name = "Fs", description ="") //TODO description?
class CmdFs extends CmdBase implements UsageAware {

    static final public NAME = 'fs'

    static final List<SubCmd> commands = new ArrayList<>()

    static {
        commands << new CmdCopy()
        commands << new CmdMove()
        commands << new CmdList()
        commands << new CmdCat()
        commands << new CmdRemove()
    }

    trait SubCmd {
        abstract int getArity()

        abstract String getName()

        abstract String getDescription()

        abstract void apply(Path source, Path target)

        String usage() {
            "Usage: nextflow fs ${name} " + (arity==1 ? "<path>" : "source_file target_file")
        }
    }

    static class CmdCopy implements SubCmd {

        @Override
        int getArity() { 2 }

        @Override
        String getName() { 'cp' }

        String getDescription() { 'Copy a file' }

        @Override
        void apply(Path source, Path target) {
            FilesEx.copyTo(source, target)
        }

    }

    static class CmdMove implements SubCmd {

        @Override
        int getArity() { 2 }

        @Override
        String getName() { 'mv' }

        String getDescription() { 'Move a file' }

        @Override
        void apply(Path source, Path target) {
            FilesEx.moveTo(source, target)
        }

    }

    static class CmdList implements SubCmd {

        @Override
        int getArity() { 1 }

        String getDescription() { 'List the content of a folder' }

        @Override
        String getName() { 'ls' }

        @Override
        void apply(Path source, Path target) {
            println source.name
        }

    }

    static class CmdCat implements SubCmd {

        @Override
        int getArity() { 1 }

        @Override
        String getName() { 'cat' }

        @Override
        String getDescription() { 'Print a file to the stdout' }

        @Override
        void apply(Path source, Path target) {
            String line
            def reader = Files.newBufferedReader(source, Charset.defaultCharset())
            while( line = reader.readLine() )
                println line
        }

    }

    static class CmdRemove implements SubCmd {

        @Override
        int getArity() { 1 }

        @Override
        String getName() { 'rm' }

        @Override
        String getDescription() { 'Remove a file' }

        @Override
        void apply(Path source, Path target) {
            Files.isDirectory(source) ? FilesEx.deleteDir(source) : FilesEx.delete(source)
        }

    }


    //@Parameter
    @CommandLine.Parameters(description = "")    //TODO mandatory? description?
    List<String> args

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
            throw new AbortOperationException("Unknow file system command: `$cmd`")
        }

        Path target
        String source
        if( cmd.arity==1 ) {
            if( args.size() < 2 )
                throw new AbortOperationException(cmd.usage())
            source = args[1]
            target = null
        }
        else {
            if( args.size() < 3 )
                throw new AbortOperationException(cmd.usage())
            source = args[1]
            target = args[2] as Path
        }

        traverse(source) { Path path -> cmd.apply(path, target) }
    }

    private SubCmd findCmd( String name ) {
        commands.find { it.name == name }
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
            result << 'Usage: nextflow fs <command> [arg]'
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
                throw new AbortOperationException("Unknown cloud sub-command: ${args[0]}")
            }
        }

    }

}
