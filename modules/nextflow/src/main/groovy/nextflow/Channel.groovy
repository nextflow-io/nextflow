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

package nextflow

import static nextflow.util.CheckHelper.*

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.CompletableFuture
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ControlMessage
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.dag.NodeMarker
import nextflow.datasource.SraExplorer
import nextflow.exception.AbortOperationException
import nextflow.extension.CH
import nextflow.extension.GroupTupleOp
import nextflow.extension.MapOp
import nextflow.file.DirListener
import nextflow.file.DirWatcher
import nextflow.file.DirWatcherV2
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import nextflow.file.PathVisitor
import nextflow.plugin.extension.PluginExtensionProvider
import nextflow.util.Duration
import nextflow.util.TestOnly
import org.codehaus.groovy.runtime.InvokerHelper
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

    @TestOnly
    private static CompletableFuture fromPath0Future

    static private Session getSession() { Global.session as Session }

    /**
     * Allow the dynamic loading of plugin provided channel extension methods
     *
     * @param name The name of the method
     * @param args The method arguments
     * @return The method return value
     */
    static def $static_methodMissing(String name, Object args) {
        PluginExtensionProvider.INSTANCE().invokeFactoryExtensionMethod(name, InvokerHelper.asArray(args))
    }

    static Object[] array0(Object value) {
        if( value instanceof Object[] )
            return (Object[])value
        if( value instanceof Collection )
            return value.toArray()
        else
            return new Object[] {value}
    }

    static private checkNoChannels(String name, Object value) {
        if( value==null )
            return
        Object[] items = array0(value)
        if( items.size()==1 && CH.isChannel(items[0]))
            throw new IllegalArgumentException("Argument of '$name' method cannot be a channel object — Likely you can replace the use of '$name' with the channel object itself")
        for( int i=0; i<items.size(); i++ ) {
            if( CH.isChannel(items[i]) )
                throw new IllegalArgumentException("Argument ${nth(i+1)} of method '$name' cannot be a channel object")
        }
    }

    static private String nth(int i) {
        if( i==1 ) return '1-st'
        if( i==2 ) return '2-nd'
        if( i==3 ) return '3-rd'
        return "$i-th"
    }

    /**
     * Create an new channel
     *
     * @return The channel instance
     */
    @Deprecated
    static DataflowChannel create() {
        if( NF.isDsl2() )
            throw new DeprecationException("Channel `create` method is not supported any more")
        return CH.queue()
    }

    static DataflowWriteChannel topic(String name) {
        if( !NF.topicChannelEnabled )
            throw new IllegalStateException("Channel.topic() requires the `nextflow.preview.topic` feature flag")
        return CH.topic(name)
    }

    /**
     * Create a empty channel i.e. only emits a STOP signal
     *
     * @return The channel instance
     */
    static DataflowWriteChannel empty() {
        final result = CH.emit(CH.queue(), STOP)
        NodeMarker.addSourceNode('Channel.empty', result)
        return result
    }

    static DataflowWriteChannel of(Object ... items) {
        checkNoChannels('channel.of', items)
        final result = CH.create()
        final values = new ArrayList()
        if( items == null ) {
            values.add(null)
        }
        else {
            for( int i=0; i<items.size(); i++ )
                addToList0(values, items[i])
        }
        values.add(STOP)
        CH.emitValues(result, values)
        NodeMarker.addSourceNode('Channel.of', result)
        return result
    }

    static private void addToList0(List list, obj) {
        if( obj instanceof Range ) {
            for( def x : obj ) {
                list.add(x)
            }
        }
        else {
            list.add(obj)
        }
    }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */
    @Deprecated
    static DataflowWriteChannel from( Collection items ) {
        final result = from0(items)
        NodeMarker.addSourceNode('Channel.from', result)
        return result
    }

    static private DataflowWriteChannel from0( Collection items ) {
        final result = CH.create()
        if( items != null )
            CH.emitAndClose(result, items)
        return result
    }

    static DataflowWriteChannel fromList( Collection items ) {
        final result = CH.create()
        CH.emitAndClose(result, items as List)
        NodeMarker.addSourceNode('Channel.fromList', result)
        return result
    }
    
    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */
    @Deprecated
    static DataflowWriteChannel from( Object... items ) {
        checkNoChannels('channel.from', items)
        for( Object it : items ) if(CH.isChannel(it))
            throw new IllegalArgumentException("channel.from argument is already a channel object")

        final result = from0(items as List)
        NodeMarker.addSourceNode('Channel.from', result)
        return result
    }

    static DataflowVariable value( obj = null ) {
        checkNoChannels('channel.value', obj)
        obj != null ? CH.value(obj) : CH.value()
    }

    /**
     * Creates a channel emitting a sequence of integers spaced by a given time interval
     *
     * @param duration
     * @return
     */
    static DataflowWriteChannel interval(String duration) {
        final result = interval0( duration, { index -> index })
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

    static DataflowWriteChannel interval(String duration, Closure closure ) {
        final result = interval0(duration, closure)
        NodeMarker.addSourceNode('Channel.interval', result)
        return result
    }

    static private DataflowWriteChannel interval0(String duration, Closure closure) {
        def millis = Duration.of(duration).toMillis()
        def timer = new Timer()
        def result = CH.create()
        long index = 0

        def task = {
            def value = closure.call(index++)
            result << value
            if( value == STOP ) {
                timer.cancel()
            }
        }

        if( NF.isDsl2() )
            session.addIgniter { timer.schedule( task as TimerTask, 0, millis ) }  
        else 
            timer.schedule( task as TimerTask, 0, millis )
        
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
    static DataflowWriteChannel<Path> fromPath( Map opts = null, pattern ) {
        if( !pattern ) throw new AbortOperationException("Missing `fromPath` parameter")
        checkNoChannels('channel.fromPath', pattern)

        // verify that the 'type' parameter has a valid value
        checkParams( 'path', opts, VALID_FROM_PATH_PARAMS )

        final result = fromPath0(opts, pattern instanceof List ? pattern : [pattern])
        NodeMarker.addSourceNode('Channel.fromPath', result)
        return result
    }

    private static DataflowWriteChannel<Path> fromPath0( Map opts, List allPatterns ) {

        final result = CH.create()
        if( NF.isDsl2() ) {
            session.addIgniter { pumpFiles0(result, opts, allPatterns) }
        }
        else {
            pumpFiles0(result, opts, allPatterns)
        }
        return result
    }
    
    private static void pumpFiles0(DataflowWriteChannel result, Map opts, List allPatterns) {
        
        def future = CompletableFuture.completedFuture(null)
        for( int index=0; index<allPatterns.size(); index++ ) {
            def factory = new PathVisitor(target: result, opts: opts)
            factory.closeChannelOnComplete = index==allPatterns.size()-1
            future = factory.applyAsync(future, allPatterns[index])
        }
        
        // abort the execution when an exception is raised
        fromPath0Future = future.exceptionally(Channel.&handlerException)
    }

    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }

    static private DataflowWriteChannel watchImpl( String syntax, String folder, String pattern, boolean skipHidden, String events, FileSystem fs ) {

        final result = CH.create()
        final legacy = System.getenv('NXF_DIRWATCHER_LEGACY')
        final DirListener watcher = legacy=='true'
                        ? new DirWatcher(syntax,folder,pattern,skipHidden,events,fs)
                        : new DirWatcherV2(syntax,folder,pattern,skipHidden,events,fs)
        watcher.onComplete { result.bind(STOP) }
        Global.onCleanup(it -> watcher.terminate())

        if( NF.isDsl2() )  {
            session.addIgniter {
                watcher.apply { Path file -> result.bind(file.toAbsolutePath()) }   
            }   
        }
        else {
            watcher.apply { Path file -> result.bind(file.toAbsolutePath()) }
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
    static DataflowWriteChannel watchPath( Pattern filePattern, String events = 'create' ) {
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
    static DataflowWriteChannel watchPath( String filePattern, String events = 'create' ) {

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

    static DataflowWriteChannel watchPath( Path path, String events = 'create' ) {
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
    static DataflowWriteChannel fromFilePairs(Map options = null, pattern) {
        checkNoChannels('channel.fromFilePairs', pattern)
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
    static DataflowWriteChannel fromFilePairs(Map options = null, pattern, Closure grouping) {
        checkNoChannels('channel.fromFilePairs', pattern)
        final allPatterns = pattern instanceof List ? pattern : [pattern]
        final allGrouping = new ArrayList(allPatterns.size())
        for( int i=0; i<allPatterns.size(); i++ ) {
            allGrouping[i] = grouping
        }

        final result = fromFilePairs0(options, allPatterns, allGrouping)
        NodeMarker.addSourceNode('Channel.fromFilePairs', result)
        return result
    }

    private static DataflowWriteChannel fromFilePairs0(Map options, List allPatterns, List<Closure> grouping) {
        assert allPatterns.size() == grouping.size()
        if( !allPatterns ) throw new AbortOperationException("Missing `fromFilePairs` parameter")
        if( !grouping ) throw new AbortOperationException("Missing `fromFilePairs` grouping parameter")

        // -- a channel from the path
        final fromOpts = fetchParams0(VALID_FROM_PATH_PARAMS, options)
        final files = new DataflowQueue()
        if( NF.isDsl2() )
            session.addIgniter { pumpFilePairs0(files,fromOpts,allPatterns) }
        else 
            pumpFilePairs0(files,fromOpts,allPatterns)

        // -- map the files to a tuple like ( ID, filePath )
        final mapper = { path, int index ->
            def prefix = grouping[index].call(path)
            return [ prefix, path ]
        }
        final mapChannel = (DataflowReadChannel)new MapOp(files, mapper)
                            .setTarget(new DataflowQueue())
                            .apply()

        boolean anyPattern=false
        for( int index=0; index<allPatterns.size(); index++ )  {
            anyPattern |= FilePatternSplitter.isMatchingPattern(allPatterns.get(index))
        }
        
        // -- result the files having the same ID        
        def DEF_SIZE = anyPattern ? 2 : 1
        def size = (options?.size ?: DEF_SIZE)
        def isFlat = options?.flat == true
        def groupOpts = [sort: true, size: size]
        def groupChannel = isFlat ? new DataflowQueue<>() : CH.create()

        new GroupTupleOp(groupOpts, mapChannel)
                .setTarget(groupChannel)
                .apply()

        // -- flat the group resulting tuples
        DataflowWriteChannel result
        if( isFlat )  {
            def makeFlat = {  id, List items ->
                def tuple = new ArrayList(items.size()+1);
                tuple.add(id)
                tuple.addAll(items)
                return tuple
            }
            result = new MapOp((DataflowReadChannel)groupChannel,makeFlat).apply()
        }
        else {
            result = groupChannel
        }

        return result
    }

    static private void pumpFilePairs0(DataflowWriteChannel files, Map fromOpts, List allPatterns) {
        def future = CompletableFuture.completedFuture(null)
        for( int index=0; index<allPatterns.size(); index++ )  {
            def factory = new PathVisitor(opts: fromOpts, target: files)
            factory.bindPayload = index
            factory.closeChannelOnComplete = index == allPatterns.size()-1
            future = factory.applyAsync( future, allPatterns.get(index) )
        }
        // abort the execution when an exception is raised
        fromPath0Future = future.exceptionally(Channel.&handlerException)
    }
    
    static private Map fetchParams0(Map valid, Map actual ) {
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
        if( indexOfWildcards==-1 && indexOfBrackets==-1 ) {
            // when the pattern does not contain any glob characters
            // then it can only have been used only to filter files with exactly
            // that name, therefore if the name is matching return it (without suffix)
            if( fileName == filePattern )
                return actual.getSimpleName()
            throw new IllegalArgumentException("Not a valid file pair globbing pattern: pattern=$filePattern file=$fileName")
        }

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

    static DataflowWriteChannel fromSRA(query) {
        checkNoChannels('channel.fromSRA', query)
        fromSRA( Collections.emptyMap(), query )
    }

    static DataflowWriteChannel fromSRA(Map opts, query) {
        checkParams('fromSRA', opts, SraExplorer.PARAMS)
        checkNoChannels('channel.fromSRA', query)

        def target = new DataflowQueue()
        def explorer = new SraExplorer(target, opts).setQuery(query)
        if( NF.isDsl2() ) {
            session.addIgniter { fetchSraFiles0(explorer) }
        }
        else {
            fetchSraFiles0(explorer)
        }

        NodeMarker.addSourceNode('Channel.fromSRA', target)
        return target
    }

    static private void fetchSraFiles0(SraExplorer explorer) {
        def future = CompletableFuture.runAsync ({ explorer.apply() } as Runnable)
        fromPath0Future = future.exceptionally(Channel.&handlerException)
    }

}
