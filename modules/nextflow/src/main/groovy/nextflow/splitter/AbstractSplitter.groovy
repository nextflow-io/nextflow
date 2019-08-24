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

package nextflow.splitter

import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPInputStream

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.exception.StopSplitIterationException
import nextflow.extension.CH
import nextflow.util.CheckHelper
/**
 * Generic data splitter, provide main methods/interfaces
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class AbstractSplitter<T> implements SplitterStrategy {

    protected Map fOptionsMap

    protected def into

    protected Closure closure

    protected boolean recordMode

    protected Map recordFields

    protected boolean autoClose = true

    protected Path sourceFile

    protected decompress

    protected String operatorName

    protected long limit

    protected Integer elem

    private targetObj

    private CollectorStrategy collector

    protected boolean multiSplit

    protected EntryCounter counter = new EntryCounter(1)

    AbstractSplitter() { }

    /**
     * Create a splitter object for the specified operator name
     *
     * @param name The name of an operator invoking the splitter. This value
     * is meant to be used only for reporting a meaningful error message
     */
    AbstractSplitter( String name ) {
        this.operatorName = name
    }

    /**
     * Create a splitter object with the specified option parameters
     *
     * See {@link #options(java.util.Map)}
     *
     * @param opt A map of named parameters
     */
    protected AbstractSplitter( Map opt ) {
        options(opt)
    }

    /**
     * @return A string representing the operator invoking the splitter
     */
    String getOperatorName() { operatorName ?: this.class.simpleName }

    /**
     * @return The splitter raw target object
     */
    protected Object getTargetObj() { targetObj }

    /**
     * @return The target object that receives the splitted chunks. It can be a {@link groovyx.gpars.dataflow.DataflowChannel} or a {@code List}
     */
    def getInto() { into }

    /**
     * @return Whenever each split is parsed to a record object or a chunk in the native format i.e. text line(s) or bytes
     */
    boolean getRecordMode() { recordMode }

    /**
     * @return The fields to be included in each parsed record
     */
    Map getRecordFields() { recordFields }

    AbstractSplitter setRecordFields( Map fields ) {
        recordMode = true
        recordFields = fields
        return this
    }

    AbstractSplitter setMultiSplit(boolean value) {
        this.multiSplit = value
        return this
    }

    /**
     * Apply the splitting operation on the given object
     *
     * @param index the current split count
     * @return Either {@link groovyx.gpars.dataflow.DataflowChannel} or a {@code List} which holds the splitted chunks
     */
    final apply() {

        def result = null
        def source = targetObj instanceof List ? findSource((List)targetObj) : targetObj

        setSource(source)


        final chunks = collector = createCollector()
        if( chunks instanceof CacheableCollector && chunks.checkCached() ) {
            log.debug "Operator `$operatorName` reusing cached chunks at path: ${chunks.getBaseFile()}"
            result = resumeFromCache(chunks)
        }

        else {
            try {
                def stream = normalizeSource(source)
                result = process(stream)
            }
            catch ( StopSplitIterationException e ) {
                log.trace 'Split iteration interrupted'
            }
        }


        /*
         * now close and return the result
         * - when the target it's a channel, send stop message
         * - when it's a list return it
         * - otherwise return the last value
         */
        if( into instanceof DataflowWriteChannel && autoClose ) {
            append(into, Channel.STOP)
            return into
        }
        if( into != null )
            return into

        return result
    }

    /**
     * @param tuple
     *      A non-empty list of objects
     * @return
     *      Returns the item at position defined by the attribute {@code #elem} in the list specified as parameter.
     *      When {@code #elem} is equal to -1 find find out the first occourence of a file object.
     *      If no file is available the first item in the list is returned
     */
    @PackageScope
    def findSource( List tuple ) {

        if( elem >= 0 )
            return tuple.get(elem)

        // find the elem-th item having Path or File type
        int pos = elem != null ? -elem : 1
        int count = 0
        for( int i=0; i<tuple.size(); i++ ) {
            def it = tuple[i]
            if(  it instanceof Path || it instanceof File ) {
                if( ++count == pos ) {
                    elem = i
                    return tuple.get(i)
                }
            }
        }

        // -- not found, if default was null just return the first item
        if( elem == null ) {
            elem = 0
            return tuple.get(0)
        }

        // -- otherwise if a `elem` value was specified but not found, raise an exception
        throw new IllegalArgumentException("Cannot find splittable file (elem=$elem)")
    }

    /**
     * Set the file object to be split (if any)
     * @param obj An instance of {@link Path} of {@code File}, otherwise do nothing
     */
    private void setSource( obj ) {
        if( obj instanceof Path )
            sourceFile = (Path)obj

        else if( obj instanceof File )
            sourceFile = (obj as File).toPath()
    }

    /**
     * Emits the cache file chunks
     *
     * @param collector
     * @param index
     * @return
     */
    protected resumeFromCache(CacheableCollector collector) {
        def result = null
        for( Path file : collector.allChunks ) {
            result = invokeEachClosure(closure, file)
        }
        return result
    }

    /**
     * Apply the splitting operation on the given object
     *
     * @param targetObject The actual object to be splitted
     * @param index the current split count
     * @return Either {@link groovyx.gpars.dataflow.DataflowChannel} or a {@code List} which holds the splitted chunks
     */
    protected abstract process( T targetObject )

    /**
     * Normalise the source object to be splitted
     *
     * @param object The object to be splitted
     * @return The normalised version of of the object to be splitted
     */
    abstract protected T normalizeSource( object )

    /**
     * Defines the splitter parameter. Subclass can override to provide format dependent options
     * <p>
     * Supporter parameters are:
     * <li>{@code by}: Defines the splitting interval e.g. how many lines are in each chunk when splitting a text file
     * <li>{@code into}: The receiving object, it can a {@link List} instance of a {@link DataflowQueue} instance
     * <li>{@code each}: The transforming closure invoke by each splitting chunk
     * <li>{@code record}:
     *          When {@code true} the splitting chunk is parsed into a record object, alternatively use to specify
     *          the field names required with a map of booleans
     * <li>{@code autoClose}:
     *          Then {@code into} parameter is {@link DataflowQueue} use this params to enable/disable the splitter
     *          to close the channel by sending a {@link nextflow.Channel#STOP} message when complete (default: {@code true})
     *
     * @param options The map holding the named parameters
     * @return The object itself
     */
    AbstractSplitter options( Map options ) {
        CheckHelper.checkParams(getOperatorName(), options, validOptions())

        fOptionsMap = options

        closure = (Closure)options.each

        if( options.by )
            counter = new EntryCounter(options.by as Integer)

        into = options.into

        recordMode = isTrueOrMap(options.record)

        if( options.record instanceof Map )
            recordFields = (Map)options.record

        if( options.autoClose instanceof Boolean )
            autoClose = options.autoClose as boolean

        if( options.decompress != null )
            decompress = options.decompress

        if( options.limit )
            limit = options.limit as long

        if( options.elem )
            elem = options.elem as int

        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    protected Map<String,Object> validOptions() {
        [
                each: Closure,
                by: Integer,
                into: [Collection, DataflowQueue, DataflowBroadcast],
                autoClose: Boolean,
                limit: Integer,
                elem: Integer,
                decompress: Boolean
        ]
    }

    /**
     * Set the target object to be splitter. This method invokes {@link #normalizeSource(java.lang.Object)}
     *
     * @param object The object to be splitted
     * @return The object itself
     */
    AbstractSplitter target( obj ) {
        targetObj = obj
        return this
    }

    /**
     * Start the slitting
     */
    def split() {
        apply()
    }

    /**
     * Apply the specified closure to each chunk in the target object
     *
     * @param closure A closure object
     */
    void each( Closure closure ) {
        this.closure = closure
        apply()
    }

    /**
     * @return The number of chunks in the target object
     */
    long count() {
        long result = 0
        closure = { result++ }
        apply()
        return result
    }

    /**
     * @return Split the target objects and return a list containing all chunks
     */
    List list() {
        into = []
        (List) apply()
    }

    /**
     * @return Split the target object and return a channel emitting the produced chunks
     */
    DataflowWriteChannel channel() {
        into = CH.create()
        (DataflowWriteChannel) apply()
    }

    /**
     * Invoke the each closure
     *
     * @param closure
     * @param chunk
     * @param index
     * @return
     */
    @PackageScope
    final invokeEachClosure( Closure closure, Object chunk ) {

        def result
        if( targetObj instanceof List ) {
            result = new ArrayList((List)targetObj)
            result.set(elem, chunk)
        }
        else {
            result = chunk
        }

        if( closure ) {
            result = closure.call(result)
        }

        if( into != null )
            append(into,result)

        return result
    }


    /**
     * Add a generic value to a target container, that can be either a {@code Collection}
     * or a {@code DataflowWriteChannel} instance
     *
     * @param into The target container, either a {@code Collection} or a {@code DataflowWriteChannel} instance
     * @param value Any value
     * @throws {@code IllegalArgumentException} whenever parameter {@code into} is not a valid object
     */

    private int debugCount = 0

    protected void append( into, value ) {
        log.trace "Splitter value: ${debugCount++}"

        if( into instanceof Collection )
            into.add(value)

        else if( into instanceof DataflowWriteChannel )
            into.bind(value)

        else
            throw new IllegalArgumentException("Not a valid 'into' target object: ${into?.class?.name}")
    }

    /**
     * @param value An object to check
     * @return {@code true} if the value is an instanceof {@link Map} of a boolean value equals to {@code true}
     */
    static protected boolean isTrueOrMap( value ) {
        if( value instanceof Map )
            return true

        return value instanceof Boolean && (value as Boolean)
    }

    /**
     * Given a {@link Path} return a new {@link InputStream} associated to it.
     * When the file name ends with {@code .gz} the stream is filtered by a {@link GZIPInputStream}
     *
     * @param path An path for an existing file
     * @return The {@link InputStream} object for the given file
     */
    protected InputStream newInputStream( Path path ) {

        def result = Files.newInputStream(path)

        if( decompress == null && path.name.endsWith('.gz') )
            decompress = true

        if( decompress ) {
            log.debug "Creating gzip splitter for: $path"
            return new GZIPInputStream(result)
        }

        return result
    }

    /**
     * @return The current {@link CollectorStrategy} object
     */
    final protected CollectorStrategy getCollector() {
        collector
    }

    /**
     * @return create a new {@link CollectorStrategy} object. Subclass must implement a valid
     * strategy
     */
    abstract protected CollectorStrategy createCollector()

}
