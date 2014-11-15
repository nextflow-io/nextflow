package nextflow.splitter
import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPInputStream

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.exception.StopSplitIterationException
import nextflow.util.CheckHelper
/**
 * Generic data splitter, provide main methods/interfaces
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class AbstractSplitter<T> implements SplitterStrategy {

    protected int count = 1

    protected def into

    protected Closure closure

    protected boolean recordMode

    protected Map recordFields

    protected boolean autoClose = true

    protected Path sourceFile

    protected meta = 'file'

    protected decompress

    protected String operatorName

    protected long limit

    private targetObj

    /**
     * Create a splitter object for the specified operator name
     *
     * @param name The name of an operator invoking the splitter. This value
     * is meant to be used only for reporting a meaningful error message
     */
    AbstractSplitter( String name = null ) {
        this.operatorName = name
    }

    /**
     * Create a splitter object with the specified option parameters
     *
     * See {@link #options(java.util.Map)}
     *
     * @param opt A map of named parameters
     */
    AbstractSplitter( Map opt ) {
        options(opt)
    }

    /**
     * @return A string representing the operator invoking the splitter
     */
    String getOperatorName() { operatorName ?: this.class.simpleName }

    /**
     * @return The number of entry of which each chunk is made up
     */
    int getCount() { count }

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

    /**
     * Apply the splitting operation on the given object
     *
     * @param index the current split count
     * @return Either {@link groovyx.gpars.dataflow.DataflowChannel} or a {@code List} which holds the splitted chunks
     */
    final apply( int index = 0 ) {
        def obj = normalizeType(targetObj)
        def result = null

        try {
            result = process(obj, index)
        }
        catch ( StopSplitIterationException e ) {
            log.trace 'Split iteration interrupted'
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
     * Apply the splitting operation on the given object
     *
     * @param targetObject The actual object to be splitted
     * @param index the current split count
     * @return Either {@link groovyx.gpars.dataflow.DataflowChannel} or a {@code List} which holds the splitted chunks
     */
    protected abstract process( T targetObject, int index )

    /**
     * Normalise the source object to be splitted
     *
     * @param object The object to be splitted
     * @return The normalised version of of the object to be splitted
     */
    abstract protected T normalizeType( object )

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

        closure = (Closure)options.each

        if( options.by )
            count = options.by as Integer

        into = options.into

        recordMode = isTrueOrMap(options.record)

        if( options.record instanceof Map )
            recordFields = (Map)options.record

        if( recordMode && count>1 )
            throw new IllegalArgumentException("When using 'record' option 'count' cannot be greater than 1")

        if( options.autoClose instanceof Boolean )
            autoClose = options.autoClose as boolean

        if( options.meta )
            meta = options.meta

        if( options.decompress != null )
            decompress = options.decompress

        if( options.limit )
            limit = options.limit as long

        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    protected Map<String,?> validOptions() {
        [
                each: Closure,
                by: Integer,
                into: [ Collection, DataflowQueue ],
                record: [ Boolean, Map ],
                autoClose: Boolean,
                meta: ['file','path','index'],
                limit: Integer,
                decompress: Boolean
        ]
    }

    /**
     * Set the target object to be splitter. This method invokes {@link #normalizeType(java.lang.Object)}
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
    DataflowQueue channel() {
        into = new DataflowQueue()
        (DataflowQueue) apply()
    }

    /**
     * Invoke the each closure
     *
     * @param closure
     * @param obj
     * @param index
     * @return
     */
    @PackageScope
    final invokeEachClosure( Closure closure, Object obj, int index ) {
        if( !closure )
            return obj

        def len = closure.getMaximumNumberOfParameters()
        if( len==1 )
            return closure.call(obj)

        return closure.call(obj, metaParam(index))
    }

    @PackageScope
    final metaParam( int index ) {

        if( meta == 'file' && sourceFile )
            return sourceFile.getName()

        if( meta == 'path' && sourceFile )
            return sourceFile

        if( meta == 'index' )
            return index

        return null
    }

    /**
     * Add a generic value to a target container, that can be either a {@code Collection}
     * or a {@code DataflowWriteChannel} instance
     *
     * @param into The target container, either a {@code Collection} or a {@code DataflowWriteChannel} instance
     * @param value Any value
     * @throws {@code IllegalArgumentException} whenever parameter {@code into} is not a valid object
     */
    protected void append( into, value ) {
        if( into instanceof Collection )
            into.add(value)

        else if( into instanceof DataflowWriteChannel )
            into.bind(value)

        else
            throw new IllegalArgumentException("Not a valid 'into' target object: ${into?.class?.name}")
    }

    static protected boolean isTrueOrMap( value ) {
        if( value instanceof Map )
            return true

        return value instanceof Boolean && (value as Boolean)
    }

    protected InputStream newInputStream( Path path ) {
        this.sourceFile = path

        def result = Files.newInputStream(path)

        if( decompress == null && path.name.endsWith('.gz') )
            decompress = true

        if( decompress ) {
            log.debug "Creating gzip splitter for: $path"
            return new GZIPInputStream(result)
        }

        return result
    }

}
