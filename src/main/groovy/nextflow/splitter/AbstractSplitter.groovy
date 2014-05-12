package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
/**
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

    AbstractSplitter options( Map options ) {
        assert options != null
        closure = (Closure)options.each
        count = options.by as Integer ?: 1

        if( options.count ) {
            log.warn "The 'count' parameter has been deprecated -- please use 'by' instead"
            count = options.count as Integer ?: 1
        }

        // add 'remainder' flag

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
     * Abstract splitter method
     *
     * @param object The object to be splitted
     * @param params The map holding the splitting named parameters
     * @param closure An option closure applied to each split entry
     * @return
     */
    AbstractSplitter target( obj ) {
        targetObj = normalizeType(obj)
        return this
    }

    def split() {
        apply(targetObj, 0)
    }

    void each( Closure closure ) {
        this.closure = closure
        apply(targetObj, 0)
    }

    long count() {
        long result = 0
        closure = { result++ }
        apply(targetObj, 0)
        return result
    }

    List list() {
        into = []
        (List) apply(targetObj, 0)
    }

    DataflowQueue channel() {
        into = new DataflowQueue()
        (DataflowQueue) apply(targetObj, 0)
    }



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
