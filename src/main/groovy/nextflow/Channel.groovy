/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import nextflow.exception.AbortOperationException
import nextflow.extension.GroupTupleOp
import nextflow.extension.MapOp
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import nextflow.util.Duration
import org.codehaus.groovy.runtime.NullObject
/**
 * Channel factory object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Channel  {

    static public ControlMessage STOP = PoisonPill.getInstance()

    static public NullObject VOID = NullObject.getNullObject()

    @Lazy
    static private Session session = { Global.session as Session }()

    /**
     * Create an new channel
     *
     * @return The channel instance
     */
    static <T> DataflowChannel<T> create() {
        new DataflowQueue()
    }

    /**
     * Create a empty channel i.e. only emits a STOP signal
     *
     * @return The channel instance
     */
    static <T> DataflowChannel<T> empty() {
        def result = new DataflowQueue()
        result.bind(STOP)
        session.dag.addSourceNode('Channel.empty', result)
        return result
    }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */
    static DataflowChannel from( Collection items ) {
        final result = Nextflow.channel(items)
        session.dag.addSourceNode('Channel.from', result)
        return result
    }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */

    static DataflowChannel from( Object... items ) {
        final result = Nextflow.channel(items)
        session.dag.addSourceNode('Channel.from', result)
        return result
    }

    /**
     * Convert an object into a *channel* variable that emits that object
     *
     * @param obj
     * @return
     */
    @Deprecated
    static DataflowVariable just( obj = null ) {
        log.warn "The operator `just` is deprecated -- Use `value` instead."
        def result = new DataflowVariable()
        if( obj != null ) result.bind(obj)
        return result
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

        final result = interval( duration, { index -> index })

        session.dag.addSourceNode('Channel.interval', result)
        return result
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
        def result = create()
        long index = 0

        def task = {
            def value = closure.call(index++)
            result << value
            if( value == STOP ) {
                timer.cancel()
            }
        }

        timer.schedule( task as TimerTask, millis )

        session.dag.addSourceNode('Channel.interval', result)
        return result
    }

    static DataflowChannel<Path> fromPath( Map opts = null, filePattern ) {
        if( !filePattern ) throw new AbortOperationException("Missing `fromPath` parameter")

        // verify that the 'type' parameter has a valid value
        checkParams( 'path', opts, VALID_PATH_PARAMS )

        final result = fromPath0(opts, filePattern)
        session.dag.addSourceNode('Channel.fromPath', result)
        return result
    }

    static private DataflowChannel<Path> fromPath0( Map options, filePattern ) {

        if( filePattern instanceof Pattern )
            return fromPathWithPattern(options, filePattern)

        else
            return fromPathWithMap(options, filePattern as Path)

    }

    static private DataflowChannel<Path> fromPathWithPattern( Map opts = null, Pattern filePattern ) {
        assert filePattern

        // split the folder and the pattern
        final splitter = FilePatternSplitter.regex().parse(filePattern.toString())
        final fs = FileHelper.fileSystemForScheme(splitter.scheme)
        pathImpl( 'regex', splitter.parent, splitter.fileName, opts, fs )
    }


    static private DataflowChannel<Path> fromPathWithMap( Map opts = null, Path filePattern ) {

        def glob = opts?.containsKey('glob') ? opts.glob as boolean : true
        if( !glob ) {
            def result = ( filePattern instanceof Path
                    ? filePattern.complete()
                    : FileHelper.asPath(filePattern.toString()).complete())
            return Nextflow.channel(result)
        }

        final fs = filePattern.getFileSystem()
        final path = filePattern.toString()
        final splitter = FilePatternSplitter.glob().parse(path)

        if( !splitter.isPattern() ) {
            return Nextflow.channel( fs.getPath( splitter.strip(path) ).complete() )
        }

        final folder = splitter.parent
        final pattern = splitter.fileName
        final result = pathImpl('glob', folder, pattern, opts, fs)
        return result
    }

    /*
     * valid parameters for fromPath operator
     */
    static private Map VALID_PATH_PARAMS = [
            type:['file','dir','any'],
            followLinks: [false, true],
            hidden: [false, true],
            maxDepth: Integer,
            exists: [false, true],
            glob: [false,true]
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
        final splitter = FilePatternSplitter.regex().parse(filePattern.toString())
        def fs = FileHelper.fileSystemForScheme(splitter.scheme)
        def result = watchImpl( 'regex', splitter.parent, splitter.fileName, false, events, fs )

        session.dag.addSourceNode('Channel.watchPath', result)
        return result
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

        final splitter = FilePatternSplitter.regex().parse(filePattern.toString())
        final fs = FileHelper.fileSystemForScheme(splitter.scheme)
        final folder = splitter.parent
        final pattern = splitter.fileName
        def result = watchImpl('glob', folder, pattern, pattern.startsWith('*'), events, fs)

        session.dag.addSourceNode('Channel.watchPath', result)
        return result
    }

    static DataflowChannel<Path> watchPath( Path path, String events = 'create' ) {
        final fs = path.getFileSystem()
        final splitter = FilePatternSplitter.regex().parse(path.toString())
        final folder = splitter.parent
        final pattern = splitter.fileName
        final result = watchImpl('glob', folder, pattern, pattern.startsWith('*'), events, fs)

        session.dag.addSourceNode('Channel.watchPath', result)
        return result
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


    static DataflowChannel fromFilePairs(Map options = null, filePattern) {
        def closure = { Path path -> readPrefix(path,filePattern) }
        fromFilePairs(options, filePattern, closure)
    }

    static DataflowChannel fromFilePairs(Map options = null, filePattern, Closure grouping) {
        if( !filePattern ) throw new AbortOperationException("Missing `fromFilePairs` parameter")
        if( !grouping ) throw new AbortOperationException("Missing `fromFilePairs` grouping parameter")

        // -- a channel from the path
        def fromOpts = fetchParams(VALID_PATH_PARAMS, options)
        def files = fromPath0(fromOpts,filePattern)

        // -- map the files to a tuple like ( ID, filePath )
        def mapper = { path ->
            def prefix = grouping.call(path)
            return [ prefix, path ]
        }
        def mapChannel = new MapOp(files, mapper).apply()

        // -- result the files having the same ID
        def size = (options?.size ?: 2)
        def groupOpts = [sort: true, size: size]
        def groupChannel = new GroupTupleOp(groupOpts, mapChannel).apply()

        // -- flat the group resulting tuples
        DataflowChannel result
        if( options?.flat == true )  {
            def makeFlat = {  id, List items ->
                def tuple = new ArrayList(items.size()+1);
                tuple.add(id)
                tuple.addAll(items)
                return tuple
            }
            result = new MapOp(groupChannel,makeFlat).apply()
        }
        else {
            result = groupChannel
        }

        session.dag.addSourceNode('Channel.fromFilePairs', result)
        return result
    }

    static private Map fetchParams( Map valid, Map actual ) {
        if( actual==null ) return null
        def result = [:]
        def keys = valid.keySet()
        keys.each {
            if( actual.containsKey(it) ) result.put(it, actual.get(it))
        }

        return result
    }

    /*
     * Helper function, given a file Path
     * returns the file name region matching a specified glob pattern
     * starting from the beginning of the name up to last matching group.
     *
     * For example:
     *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
     *
     * Returns:
     *   'file_alpha'
     */
    @PackageScope
    static String readPrefix( Path actual, template ) {

        final fileName = actual.getFileName().toString()

        def filePattern = template.toString()
        int p = filePattern.lastIndexOf('/')
        if( p != -1 )
            filePattern = filePattern.substring(p+1)
        if( !filePattern.contains('*') && !filePattern.contains('?') && !filePattern.contains('{') )
            filePattern = '*' + filePattern

        int q = filePattern.indexOf('{')
        if( q == -1 ) q = filePattern.indexOf('[')

        def regex = filePattern
                .replace('.','\\.')
                .replace('*','(.*)')
                .replace('?','(.?)')
                .replace('{','(?:')
                .replace('}',')')
                .replace(',','|')

        def matcher = (fileName =~ /$regex/)
        if( matcher.matches() ) {
            def c=matcher.groupCount()
            def end = c ? matcher.end(c) : ( q != -1 ? q : fileName.size() )
            def prefix = fileName.substring(0,end)
            while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
                prefix=prefix[0..-2]

            return prefix
        }

        return null
    }

}