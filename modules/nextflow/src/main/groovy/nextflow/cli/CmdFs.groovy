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
 */

package nextflow.cli

import static nextflow.file.FileHelper.toCanonicalPath

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.attribute.BasicFileAttributes

import com.beust.jcommander.Parameter
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
class CmdFs extends CmdBase implements UsageAware {

    static final public NAME = 'fs'

    static final List<SubCmd> commands = new ArrayList<>()

    static {
        commands << new CmdCopy()
        commands << new CmdMove()
        commands << new CmdList()
        commands << new CmdCat()
        commands << new CmdRemove()
        commands << new CmdStat()
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

    static class CmdStat implements SubCmd {
        @Override
        int getArity() { 1 }

        @Override
        String getName() { 'stat' }

        @Override
        String getDescription() { 'Print file to meta info' }

        @Override
        void apply(Path source, Path target) {
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

    private Session createSession() {
        // create the config
        final config = new ConfigBuilder()
                .setOptions(getLauncher().getOptions())
                .setBaseDir(Paths.get('.'))
                .build()

        return new Session(config)
    }

    @Override
    void run() {
        if( !args ) {
            usage()
            return
        }

        Plugins.init()
        final session = createSession()
        try {
            run0()
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

    private void run0() {
        final cmd = findCmd(args[0])
        if( !cmd ) {
            throw new AbortOperationException("Unknown fs sub-command: `$cmd`")
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
                throw new AbortOperationException("Unknown fs sub-command: ${args[0]}")
            }
        }

    }

}
