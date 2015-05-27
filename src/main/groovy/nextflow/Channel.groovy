/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow
import static java.nio.file.StandardWatchEventKinds.ENTRY_CREATE
import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE
import static java.nio.file.StandardWatchEventKinds.ENTRY_MODIFY
import static java.nio.file.StandardWatchEventKinds.OVERFLOW
import static nextflow.util.CheckHelper.checkParams

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService
import java.util.regex.Pattern

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.ControlMessage
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.file.FileHelper
import nextflow.util.Duration
import org.codehaus.groovy.runtime.NullObject
/**
 * Channel factory object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Channel  {

    static ControlMessage STOP = PoisonPill.getInstance()

    static NullObject VOID = NullObject.getNullObject()

    /**
     * Create an new channel
     *
     * @return The channel instance
     */
    static <T> DataflowChannel<T> create() { new DataflowQueue() }

    /**
     * Create a empty channel i.e. only emits a STOP signal
     *
     * @return The channel instance
     */
    static <T> DataflowChannel<T> empty() {
        def result = new DataflowQueue()
        result.bind(STOP)
        return result
    }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */
    static DataflowChannel from( Collection items ) { Nextflow.channel(items) }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */

    static DataflowChannel from( Object... items ) { Nextflow.channel(items) }

    /**
     * Convert an object into a *channel* variable that emits that object
     *
     * @param obj
     * @return
     */
    static DataflowVariable just( obj = null ) {
        value(obj)
    }

    static DataflowVariable value( obj = null ) {

        def result = new DataflowVariable()
        if( obj != null ) result.bind(obj)

        return result
    }

    /**
     * Creates a channel emitting a sequence of integers spaced by a given time interval
     *
     * @param duration
     * @return
     */
    static DataflowChannel interval(String duration) {

        interval( duration, { index -> index })

    }

    /**
     * Creates a channel emitting a sequence of value given by the closure and spaced by a given time interval.
     *
     * To stop the interval return the special value {@code #STOP}
     *
     * @param duration
     * @return
     */

    static DataflowChannel interval(String duration, Closure closure ) {

        def millis = Duration.of(duration).toMillis()
        def timer = new Timer()
        def channel = create()
        long index = 0

        def task = {
            def value = closure.call(index++)
            channel << value
            if( value == STOP ) {
                timer.cancel()
            }
        }

        timer.schedule( task as TimerTask, millis )

        return channel
    }

    static DataflowChannel<Path> fromPath( Map options = null, filePattern ) {
        assert filePattern

        if( filePattern instanceof Pattern )
            fromPath(options, filePattern)

        else if( filePattern instanceof Path )
            fromPath(options, filePattern)

        else
            fromPath(options, filePattern.toString() as Path)
    }

    static DataflowChannel<Path> fromPath( Map opts = null, Pattern filePattern ) {
        assert filePattern

        // split the folder and the pattern
        def ( String folder, String pattern, String scheme ) = FileHelper.getFolderAndPattern(filePattern.toString())
        def fs = FileHelper.fileSystemForScheme(scheme)
        pathImpl( 'regex', folder, pattern, opts, fs )
    }


    static DataflowChannel<Path> fromPath( Map options = null, Path path ) {
        assert path

        final fs = path.getFileSystem()
        def filePattern = path.toString()

        if( !FileHelper.isGlobPattern(filePattern)) {
            return from( path.complete() )
        }

        // split the folder and the pattern
        def ( String folder, String pattern ) = FileHelper.getFolderAndPattern(filePattern)

        pathImpl('glob', folder, pattern, options, fs)
    }

    @Deprecated
    static DataflowChannel<Path> path(Map options = null, filePattern ) {
        log.warn "Operator 'path' has been deprecated -- Use operator 'fromPath' instead"
        fromPath(options,filePattern)
    }

    @Deprecated
    static DataflowChannel<Path> path(Map options = null, Pattern filePattern ) {
        log.warn "Operator 'path' has been deprecated -- Use operator 'fromPath' instead"
        fromPath(options,filePattern)
    }

    @Deprecated
    static DataflowChannel<Path> path(Map options = null, String filePattern ) {
        log.warn "Operator 'path' has been deprecated -- Use operator 'fromPath' instead"
        fromPath(options,filePattern)
    }

    /*
     * valid parameters for fromPath operator
     */
    static private Map VALID_PATH_PARAMS = [
            type:['file','dir','any'],
            followLinks: [false, true],
            hidden: [false, true],
            maxDepth: Integer
            ]

    /**
     * Implement the logic for files matching
     *
     * @param syntax The "syntax" to match file names, either {@code regex} or {@code glob}
     * @param folder The parent folder
     * @param pattern The file name pattern
     * @param skipHidden Whenever skip the hidden files
     * @return A dataflow channel instance emitting the file matching the specified criteria
     */
    static private DataflowChannel<Path> pathImpl( String syntax, String folder, String pattern, Map opts, FileSystem fs )  {
        assert syntax in ['regex','glob']
        log.debug "files for syntax: $syntax; folder: $folder; pattern: $pattern; options: ${opts}"

        // verify that the 'type' parameter has a valid value
        checkParams( 'path', opts, VALID_PATH_PARAMS )

        // now apply glob file search
        final path = fs.getPath(folder).complete()
        final channel = new DataflowQueue<Path>()

        if( opts == null )
            opts = [:]

        // set the 'matcher' syntax: 'regex' or 'glob'
        opts.syntax = syntax

        // set the 'matcher' type: 'file', 'dir' or 'any' (default: file)
        if( !opts.type )
            opts.type = 'file'

        Thread.start {
            try {
                FileHelper.visitFiles(opts, path, pattern) { Path file -> channel.bind(file) }
            }
            catch (NoSuchFileException e) {
                log.debug "No such file: $folder -- Skipping visit"
            }
            finally {
                channel << STOP
            }
        }

        return channel
    }


    static private DataflowChannel<Path> watchImpl( String syntax, String folder, String pattern, boolean skipHidden, String events, FileSystem fs ) {
        assert syntax in ['regex','glob']
        log.debug "Watch service for path: $folder; syntax: $syntax; pattern: $pattern; skipHidden: $skipHidden; events: $events"

        // now apply glob file search
        final result = create()
        final path = fs.getPath(folder)
        if( !path.isDirectory() ) {
            log.warn "Cannot watch a not existing path: $path -- Make sure that path exists and it is a directory"
            result.bind(STOP)
            return result
        }

        final rule = "$syntax:${folder}${pattern}"
        final matcher = FileHelper.getPathMatcherFor(rule, path.fileSystem)
        final eventsToWatch = stringToWatchEvents(events)

        Thread.start {
            WatchService watcher = path.getFileSystem().newWatchService()
            path.register(watcher, eventsToWatch)

            while( true ) {
                // wait for key to be signaled
                try {
                    WatchKey key = watcher.take();

                    for (WatchEvent<?> event: key.pollEvents()) {
                        WatchEvent.Kind<?> kind = event.kind();

                        if( kind == OVERFLOW ) {
                            log.debug "Watcher on path > $path -- get a OVERFLOW event"
                            continue
                        }

                        // The filename is the context of the event.
                        Path fileName = (event as WatchEvent<Path>).context();
                        log.trace "Watcher path > $path -- event: $kind; fileName: $fileName"
                        Path target = path.resolve(fileName)

                        if (matcher.matches(target) && ( !skipHidden || !Files.isHidden(target) )) {
                            log.trace "File watcher: $target matching: $rule -- event: $kind"
                            result.bind(target.toAbsolutePath())
                        }

                    }

                    // Reset the key -- this step is critical if you want to
                    // receive further watch events.  If the key is no longer valid,
                    // the directory is inaccessible so exit the loop.
                    boolean valid = key.reset();
                    if (!valid) {
                        break;
                    }
                }
                catch (Exception e) {
                    log.debug "Exception while watching path: $path", e
                    return;
                }

            }
        }

        return result
    }


    /**
     * Watch the a folder for the specified events emitting the files that matches
     * the specified regular expression.
     *
     *
     * @param filePattern
     *          The file pattern to match e.g. /*.fasta/
     *
     * @param events
     *          The events to watch, a comma separated string of the following values:
     *          {@code create}, {@code modify}, {@code delete}
     *
     * @return  A dataflow channel that will emit the matching files
     *
     */
    static DataflowChannel<Path> watchPath( Pattern filePattern, String events = 'create' ) {
        assert filePattern
        // split the folder and the pattern
        def ( String folder, String pattern, String scheme ) = FileHelper.getFolderAndPattern(filePattern.toString())
        def fs = FileHelper.fileSystemForScheme(scheme)
        watchImpl( 'regex', folder, pattern, false, events, fs )
    }

    /**
     * Watch the a folder for the specified events emitting the files that matches
     * the specified {@code glob} pattern.
     *
     * @link http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
     *
     * @param filePattern
     *          The file pattern to match e.g. /some/path/*.fasta
     *
     * @param events
     *          The events to watch, a comma separated string of the following values:
     *          {@code create}, {@code modify}, {@code delete}
     *
     * @return  A dataflow channel that will emit the matching files
     *
     */
    static DataflowChannel<Path> watchPath( String filePattern, String events = 'create' ) {
        def ( String folder, String pattern, String scheme ) = FileHelper.getFolderAndPattern(filePattern)
        def fs = FileHelper.fileSystemForScheme(scheme)
        watchImpl('glob', folder, pattern, pattern.startsWith('*'), events, fs)
    }

    static DataflowChannel<Path> watchPath( Path path, String events = 'create' ) {
        def fs = path.getFileSystem()
        def ( String folder, String pattern ) = FileHelper.getFolderAndPattern(path.toString())
        watchImpl('glob', folder, pattern, pattern.startsWith('*'), events, fs)
    }


    static private EVENT_MAP = [
            'create':ENTRY_CREATE,
            'delete':ENTRY_DELETE,
            'modify':ENTRY_MODIFY
    ]

    /**
     * Converts a comma separated events string to the corresponding {@code WatchEvent.Kind} instances
     *
     * @param events the list of events to watch
     * @return
     */
    @PackageScope
    static WatchEvent.Kind<Path>[]  stringToWatchEvents(String events = null){
        def result = []
        if( !events )
            result << ENTRY_CREATE

        else {
            events.split(',').each {
                def ev = it.trim().toLowerCase()
                def val = EVENT_MAP[ev]
                if( !val )
                    throw new IllegalArgumentException("Invalid watch event: $it -- Valid values are: ${EVENT_MAP.keySet().join(', ')}")
                result << val
            }
        }

        result as WatchEvent.Kind<Path>[]

    }

}