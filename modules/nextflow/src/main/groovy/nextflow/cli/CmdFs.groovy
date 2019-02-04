/*
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

package nextflow.cli
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.Parameter
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
/**
 * Implements `fs` command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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


    @Parameter
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
