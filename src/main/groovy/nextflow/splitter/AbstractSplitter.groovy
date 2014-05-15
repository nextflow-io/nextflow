package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
/**
 * Generic data splitter, provide main methods/interfaces
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class AbstractSplitter<T> implements SplitterStrategy {

    protected int count

    protected def into

    protected Closure closure

    protected boolean recordMode

    protected Map recordCols

    protected boolean autoClose = true

    private T targetObj

    AbstractSplitter( Map opt = [:] ) {
        options(opt)
    }

    int getCount() { count }

    def getInto() { into }

    Map getRecordCols() { recordCols }

    boolean getRecordMode() { recordMode }

    abstract apply( T targetObject, int index )

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
        assert options != null
        closure = (Closure)options.each
        count = options.by as Integer ?: 1

        if( options.count ) {
            log.warn "The 'count' parameter has been deprecated -- please use 'by' instead"
            count = options.count as Integer ?: 1
        }

        //TODO add 'remainder' flag

        into = options.into
        if( into && !(into instanceof Collection) && !(into instanceof DataflowQueue) )
            throw new IllegalArgumentException("Argument 'into' can be a subclass of Collection or a DataflowQueue type -- Entered value type: ${into.class.name}")

        recordMode = isTrueOrMap(options.record)

        if( options.record instanceof Map )
            recordCols = (Map)options.record

        if( recordMode && count>1 )
            throw new IllegalArgumentException("When using 'record' option 'count' cannot be greater than 1")

        if( options.autoClose instanceof Boolean )
            autoClose = options.autoClose as boolean

        return this
    }


    /**
     * Set the target object to be splitter. This method invokes {@link #normalizeType(java.lang.Object)}
     *
     * @param object The object to be splitted
     * @return The object itself
     */
    AbstractSplitter target( obj ) {
        targetObj = normalizeType(obj)
        return this
    }

    /**
     * Start the slitting
     */
    def split() {
        apply(targetObj, 0)
    }

    /**
     * Apply the specified closure to each chunk in the target object
     *
     * @param closure A closure object
     */
    void each( Closure closure ) {
        this.closure = closure
        apply(targetObj, 0)
    }

    /**
     * @return The number of chunks in the target object
     */
    long count() {
        long result = 0
        closure = { result++ }
        apply(targetObj, 0)
        return result
    }

    /**
     * @return Split the target objects and return a list containing all chunks
     */
    List list() {
        into = []
        (List) apply(targetObj, 0)
    }

    /**
     * @return Split the target object and return a channel emitting the produced chunks
     */
    DataflowQueue channel() {
        into = new DataflowQueue()
        (DataflowQueue) apply(targetObj, 0)
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
    static invokeEachClosure( Closure closure, Object obj, int index ) {
        if( !closure ) return obj
        def types = closure.getParameterTypes()
        if( types.size()>1 ) {
            return closure.call(obj, index)
        }
        else {
            return closure.call(obj)
        }
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

}
