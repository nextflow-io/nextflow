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

package nextflow

import static nextflow.util.CheckHelper.checkParams

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.ControlMessage
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.dag.NodeMarker
import nextflow.exception.AbortOperationException
import nextflow.extension.GroupTupleOp
import nextflow.extension.MapOp
import nextflow.file.DirWatcher
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import nextflow.file.PathVisitor
import nextflow.util.Duration
import org.codehaus.groovy.runtime.NullObject
/**
 * Channel factory object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Channel  {

    static public ControlMessage STOP = PoisonPill.getInstance()

    static public NullObject VOID = NullObject.getNullObject()

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
        NodeMarker.addSourceNode('Channel.empty', result)
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
        NodeMarker.addSourceNode('Channel.from', result)
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
        NodeMarker.addSourceNode('Channel.from', result)
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

        NodeMarker.addSourceNode('Channel.interval', result)
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

        NodeMarker.addSourceNode('Channel.interval', result)
        return result
    }

    /*
     * valid parameters for fromPath operator
     */
    static private Map VALID_FROM_PATH_PARAMS = [
            type:['file','dir','any'],
            followLinks: Boolean,
            hidden: Boolean,
            maxDepth: Integer,
            checkIfExists: Boolean,
            glob: Boolean,
            relative: Boolean
    ]

    /**
     * Implements {@code Channel.fromPath} factory method
     *
     * @param opts
     *      A map object holding the optional method parameters
     * @param pattern
     *      A file path or a glob pattern matching the required files.
     *      Multiple paths or patterns can be using a list object
     * @return
     *      A channel emitting the matching files
     */
    static DataflowChannel<Path> fromPath( Map opts = null, pattern ) {
        if( !pattern ) throw new AbortOperationException("Missing `fromPath` parameter")

        // verify that the 'type' parameter has a valid value
        checkParams( 'path', opts, VALID_FROM_PATH_PARAMS )

        final result = fromPath0(opts, pattern instanceof List ? pattern : [pattern])
        NodeMarker.addSourceNode('Channel.fromPath', result)
        return result
    }

    private static DataflowChannel<Path> fromPath0( Map opts, List allPatterns ) {

        final result = new DataflowQueue()
        for( int i=0; i<allPatterns.size(); i++ ) {
            def factory = new PathVisitor(target: result, opts: opts)
            factory.closeChannelOnComplete = i==allPatterns.size()-1
            factory.apply(allPatterns[i])
        }

        return result
    }


    static private DataflowChannel<Path> watchImpl( String syntax, String folder, String pattern, boolean skipHidden, String events, FileSystem fs ) {
        final result = create()

        new DirWatcher(syntax,folder,pattern,skipHidden,events, fs)
                .setOnComplete { result.bind(STOP) }
                .apply { Path file -> result.bind(file.toAbsolutePath()) }

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

        NodeMarker.addSourceNode('Channel.watchPath', result)
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

        if( filePattern.endsWith('/') )
            filePattern += '*'
        else if(Files.isDirectory(Paths.get(filePattern)))
            filePattern += '/*'

        final splitter = FilePatternSplitter.glob().parse(filePattern)
        final fs = FileHelper.fileSystemForScheme(splitter.scheme)
        final folder = splitter.parent
        final pattern = splitter.fileName
        def result = watchImpl('glob', folder, pattern, pattern.startsWith('*'), events, fs)

        NodeMarker.addSourceNode('Channel.watchPath', result)
        return result
    }

    static DataflowChannel<Path> watchPath( Path path, String events = 'create' ) {
        final fs = path.getFileSystem()
        final splitter = FilePatternSplitter.glob().parse(path.toString())
        final folder = splitter.parent
        final pattern = splitter.fileName
        final result = watchImpl('glob', folder, pattern, pattern.startsWith('*'), events, fs)

        NodeMarker.addSourceNode('Channel.watchPath', result)
        return result
    }

    /**
     * Implements the `fromFilePairs` channel factory method
     * 
     * @param options
     *      A {@link Map} holding the optional parameters
     *      - type: either `file`, `dir` or `any`
     *      - followLinks: Boolean
     *      - hidden: Boolean
     *      - maxDepth: Integer
     *      - glob: Boolean
     *      - relative: Boolean
     * @param pattern
     *      One or more path patterns eg. `/some/path/*_{1,2}.fastq`
     * @return
     *      A channel emitting the file pairs matching the specified pattern(s)
     */
    static DataflowChannel fromFilePairs(Map options = null, pattern) {
        final allPatterns = pattern instanceof List ? pattern : [pattern]
        final allGrouping = new ArrayList(allPatterns.size())
        for( int i=0; i<allPatterns.size(); i++ ) {
            final template = allPatterns[i]
            allGrouping[i] = { Path path -> readPrefix(path,template) }
        }

        final result = fromFilePairs0(options, allPatterns, allGrouping)
        NodeMarker.addSourceNode('Channel.fromFilePairs', result)
        return result
    }

    /**
     * Implements the `fromFilePairs` channel factory method
     *
     * @param options
     *      A {@link Map} holding the optional parameters
     *      - type: either `file`, `dir` or `any`
     *      - followLinks: Boolean
     *      - hidden: Boolean
     *      - maxDepth: Integer
     *      - glob: Boolean
     *      - relative: Boolean
     * @param pattern
     *      One or more path patterns eg. `/some/path/*_{1,2}.fastq`
     * @param grouping
     *      A closure implementing a pair grouping rule for the specified
     *      file patterns
     * @return
     *      A channel emitting the file pairs matching the specified pattern(s)
     */
    static DataflowChannel fromFilePairs(Map options = null, pattern, Closure grouping) {
        final allPatterns = pattern instanceof List ? pattern : [pattern]
        final allGrouping = new ArrayList(allPatterns.size())
        for( int i=0; i<allPatterns.size(); i++ ) {
            allGrouping[i] = grouping
        }

        final result = fromFilePairs0(options, allPatterns, allGrouping)
        NodeMarker.addSourceNode('Channel.fromFilePairs', result)
        return result
    }

    private static DataflowChannel fromFilePairs0(Map options, List allPatterns, List<Closure> grouping) {
        assert allPatterns.size() == grouping.size()
        if( !allPatterns ) throw new AbortOperationException("Missing `fromFilePairs` parameter")
        if( !grouping ) throw new AbortOperationException("Missing `fromFilePairs` grouping parameter")

        // -- a channel from the path
        final fromOpts = fetchParams(VALID_FROM_PATH_PARAMS, options)
        final files = new DataflowQueue()
        for( int index=0; index<allPatterns.size(); index++ )  {
            def factory = new PathVisitor(opts: fromOpts, target: files)
            factory.bindPayload = index
            factory.closeChannelOnComplete = index == allPatterns.size()-1
            factory.apply( allPatterns.get(index) )
        }

        // -- map the files to a tuple like ( ID, filePath )
        def mapper = { path, int index ->
            def prefix = grouping[index].call(path)
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

        final indexOfWildcards = filePattern.findIndexOf { it=='*' || it=='?' }
        final indexOfBrackets = filePattern.findIndexOf { it=='{' || it=='[' }
        if( indexOfWildcards==-1 && indexOfBrackets==-1 )
            filePattern = '*' + filePattern

        // count the `*` and `?` wildcard before any {} and [] glob pattern
        int groupCount = 0
        for( int i=0; i<filePattern.size(); i++ ) {
            def ch = filePattern[i]
            if( ch=='?' || ch=='*' )
                groupCount++
            else if( ch=='{' || ch=='[' )
                break
        }

        def regex = filePattern
                .replace('.','\\.')
                .replace('*','(.*)')
                .replace('?','(.?)')
                .replace('{','(?:')
                .replace('}',')')
                .replace(',','|')

        def matcher = (fileName =~ /$regex/)
        if( matcher.matches() ) {
            def c=Math.min(groupCount, matcher.groupCount())
            def end = c ? matcher.end(c) : ( indexOfBrackets != -1 ? indexOfBrackets : fileName.size() )
            def prefix = fileName.substring(0,end)
            while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
                prefix=prefix[0..-2]

            return prefix
        }

        return null
    }

}
