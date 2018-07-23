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

package nextflow.file

import static nextflow.Channel.STOP
import static nextflow.file.FileHelper.fileSystemForScheme
import static nextflow.file.FileHelper.isGlobAllowed
import static nextflow.file.FileHelper.visitFiles

import java.nio.file.FileSystem
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Global
import nextflow.Session
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Implements the logic for {@code Channel.fromPath} and {@code Channel.fromFilePairs}
 * factory methods
 */
@CompileStatic
class PathVisitor {

    private static Logger log = LoggerFactory.getLogger(PathVisitor)

    @Lazy
    static ExecutorService singleThread = createExecutorService()

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
            target << filePattern.complete()
            if( closeChannelOnComplete ) target << STOP
            return
        }

        final fs = filePattern.getFileSystem()
        final path = filePattern.toString()
        final splitter = FilePatternSplitter.glob().parse(path)

        if( !splitter.isPattern() ) {
            target << fs.getPath( splitter.strip(path) ).complete()
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

        singleThread.execute {
            try {
                visitFiles(opts, path, pattern) { Path file ->
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
                if( closeChannelOnComplete )
                    target.bind(STOP)
            }
        }
    }

    static private ExecutorService createExecutorService() {
        final result = Executors.newSingleThreadExecutor()
        def session = Global.session as Session
        if( !session ) return result

        session.onShutdown {
            result.shutdown()
            result.awaitTermination(5,TimeUnit.MINUTES)
        }
        
        return result
    }

}