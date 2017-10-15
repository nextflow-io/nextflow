package nextflow.extension
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.splitter.AbstractSplitter
import nextflow.splitter.FastqSplitter
import nextflow.splitter.SplitterFactory
/**
 * Implements splitter operators:
 * - splitText
 * - splitCsv
 * - splitFasta
 * - splitFastq
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SplitOp {

    /**
     * The channel to which this operator is applied
     */
    private DataflowReadChannel source

    /**
     * The channel emitting the splitter chunks
     */
    private DataflowQueue output

    /**
     * Index of the elements to which a split operation need to be applied
     */
    private List<Integer> indexes

    /**
     * The name of the operator eg. {code splitFasta}
     */
    private String methodName

    /**
     * Operator named parameters
     */
    private Map params

    /**
     * Whenever the splitter is applied to a paired-end read files (only valid for {@code splitFastaq} operator.
     */
    private boolean pairedEnd

    /**
     * Creates a splitter operator
     *
     * @param source The source channel to which apply to operator
     * @param methodName The operator method name eg. {@code splitFasta}, {@code splitCsv}, etc.
     * @param opts The operator named options
     */
    SplitOp( DataflowReadChannel source, String methodName, Map opts ) {

        this.source = source
        this.params = opts != null ? new HashMap(opts) : new HashMap<>()
        this.methodName = methodName

        if( params.pe && methodName != 'splitFastq' )
            throw new IllegalArgumentException("Unknown argument 'pe' for operator 'splitFastq'")

        if( params.pe==true && params.elem )
            throw new IllegalArgumentException("Parameter `pe` and `elem` conflicts")

        if( params.pe == true ) {
            params.remove('pe') // <-- remove to avoid entering in a loop
            indexes = [-1,-2]
        }
        if( params.elem instanceof List<Integer> )
            indexes = params.elem as List<Integer>

        // -- validate options
        if( params.autoClose == true )
            throw new IllegalArgumentException('Parameter `autoClose` do not supported')
        // turn off channel auto-close
        params.autoClose = false

        if( params.into && !(params.into instanceof DataflowQueue) )
            throw new IllegalArgumentException('Parameter `into` must reference a channel object')

    }

    SplitOp setSplitIndex( int value ) {
        params.elem = value
        return this
    }

    DataflowQueue apply() {
        indexes ? splitMultiEntries() : splitSingleEntry()
        return output
    }

    /**
     * Split more than one elements. Each split operation is handled
     * on a separate channel. All channels are then merged to a
     * single output result channel.
     */
    private void splitMultiEntries() {

        final cardinality = indexes.size()

        // creates a copy of `source` channel for each element to split
        def copies = new IntoOp(source, cardinality).apply().getOutputs()

        // applies the splitter the each channel copy
        def splitted = new ArrayList(cardinality)
        for( int i=0; i<cardinality; i++ ) {
            def channel = (DataflowQueue)copies.get(i)
            def op = new SplitOp(channel, methodName, params)
            op.splitIndex = indexes.get(i)
            op.pairedEnd = true
            splitted.add( op.apply() )
        }

        // now merge the result
        this.output = new DataflowQueue()
        DataflowHelper.newOperator(splitted, [output], new SplitterMergeClosure(indexes))
    }

    /**
     * Apply the split operation to a single element
     */
    private void splitSingleEntry() {

        // create a new DataflowChannel that will receive the splitter entries
        if( params.into instanceof DataflowQueue ) {
            output = (DataflowQueue)params.into
        }
        else {
            output = params.into = new DataflowQueue<>()
        }

        // create the splitter and set the options
        def splitter = SplitterFactory
                .create(methodName)
                .options(params) as AbstractSplitter

        if( pairedEnd ) {
            (splitter as FastqSplitter).emitSplitIndex = true
        }

        DataflowHelper.subscribeImpl ( source, [
                onNext: { entry -> splitter.target(entry).apply() },
                onComplete: { output << Channel.STOP }
        ] )

    }

}
