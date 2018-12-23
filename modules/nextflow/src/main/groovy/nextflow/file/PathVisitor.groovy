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

package nextflow.file

import java.nio.file.FileSystem
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.concurrent.CompletableFuture
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Global
import nextflow.Session
import nextflow.util.CustomThreadFactory
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import static nextflow.Channel.STOP
import static nextflow.file.FileHelper.fileSystemForScheme
import static nextflow.file.FileHelper.isGlobAllowed
import static nextflow.file.FileHelper.visitFiles
/**
 * Implements the logic for {@code Channel.fromPath} and {@code Channel.fromFilePairs}
 * factory methods
 */
@CompileStatic
class PathVisitor {

    private static Logger log = LoggerFactory.getLogger(PathVisitor)

    DataflowQueue target

    boolean closeChannelOnComplete = true

    Map opts

    def bindPayload

    DataflowQueue apply(Object filePattern) {
        if( opts == null )
            opts = [:]

        if( !target )
            target = new DataflowQueue()

        if( filePattern instanceof Pattern )
            applyRegexPattern0(filePattern)

        else
            applyGlobPattern0(filePattern as Path)

        return target
    }

    CompletableFuture applyAsync(final CompletableFuture future, final Object filePattern ) {
        future.thenRunAsync({ apply(filePattern) } as Runnable, executor)
    }

    private void applyRegexPattern0( Pattern filePattern ) {
        assert filePattern
        // split the folder and the pattern
        final splitter = FilePatternSplitter.regex().parse(filePattern.toString())
        final fs = fileSystemForScheme(splitter.scheme)
        pathImpl( 'regex', splitter.parent, splitter.fileName, fs )
    }


    private void applyGlobPattern0(Path filePattern) {

        final glob = opts?.containsKey('glob') ? opts.glob as boolean : isGlobAllowed(filePattern)
        if( !glob ) {
            target << FileHelper.checkIfExists(filePattern, opts)
            if( closeChannelOnComplete ) target << STOP
            return
        }

        final fs = filePattern.getFileSystem()
        final path = filePattern.toString()
        final splitter = FilePatternSplitter.glob().parse(path)

        if( !splitter.isPattern() ) {
            final result = fs.getPath( splitter.strip(path) )
            target << FileHelper.checkIfExists(result, opts)
            if( closeChannelOnComplete ) target << STOP
            return
        }

        final folder = splitter.parent
        final pattern = splitter.fileName
        pathImpl('glob', folder, pattern, fs)
    }


    /**
     * Implement the logic for files matching
     *
     * @param syntax The "syntax" to match file names, either {@code regex} or {@code glob}
     * @param folder The parent folder
     * @param pattern The file name pattern
     * @param skipHidden Whenever skip the hidden files
     * @return A dataflow channel instance emitting the file matching the specified criteria
     */
    private void pathImpl(String syntax, String folder, String pattern, FileSystem fs )  {
        assert syntax in ['regex','glob']
        log.debug "files for syntax: $syntax; folder: $folder; pattern: $pattern; options: ${opts}"

        // now apply glob file search
        final path = fs.getPath(folder).complete()

        if( opts == null )
            opts = [:]

        // set the 'matcher' syntax: 'regex' or 'glob'
        opts.syntax = syntax

        // set the 'matcher' type: 'file', 'dir' or 'any' (default: file)
        if( !opts.type )
            opts.type = 'file'

        int count=0
        try {
            visitFiles(opts, path, pattern) { Path file ->
                count++
                if( bindPayload == null )
                    target.bind(file)

                else {
                    def pair = new ArrayList<>(2)
                    pair[0] = file
                    pair[1] = bindPayload
                    target.bind(pair)
                }
            }
        }
        catch (NoSuchFileException e) {
            log.debug "No such file: $folder -- Skipping visit"
        }
        finally {
            if( !count && opts.checkIfExists as boolean )
                throw new IllegalArgumentException("No files match pattern `$pattern` at path: $folder")
            if( closeChannelOnComplete )
                target.bind(STOP)
        }

    }


    @PackageScope
    static ExecutorService getExecutor() {
        createExecutor(Global.session as Session)
    }

    // note: the memoized annotation guarantee that for the same session
    // it return the same ExecutorService instance
    @Memoized
    @PackageScope
    static ExecutorService createExecutor(Session session) {
        final result = Executors.newCachedThreadPool(new CustomThreadFactory('PathVisitor'))
        session?.onShutdown { result.shutdown() }
        return result
    }

}