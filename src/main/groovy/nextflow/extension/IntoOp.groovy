package nextflow.extension

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
/**
 * Implements the {@link DataflowExtensions#into} operators logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IntoOp {

    private DataflowReadChannel source

    private List<DataflowWriteChannel> outputs

    private Session session = (Session)Global.session


    IntoOp( DataflowReadChannel source, List<DataflowWriteChannel> targets ) {
        assert source
        assert targets

        this.source = source
        this.outputs = targets
    }

    IntoOp( DataflowReadChannel source, int n ) {
        assert source
        assert n

        def targets = new ArrayList(n)
        for( int i=0; i<n; i++ )
            targets << new DataflowQueue()

        this.source = source
        this.outputs = targets
    }

    IntoOp( DataflowReadChannel source, Closure holder ) {
        assert source
        assert holder

        final names = CaptureProperties.capture(holder)
        final binding = Global.session.binding

        def targets = []
        names.each { identifier ->
            def channel = DataflowExtensions.newChannelBy(source)
            targets.add(channel)
            binding.setVariable(identifier, channel)
        }

        this.source = source
        this.outputs = targets
    }

    List<DataflowWriteChannel> getOutputs() { outputs }

    IntoOp apply() {

        final params = [:]
        params.inputs = [source]
        params.outputs = outputs
        params.listeners = createListener()

        DataflowExtensions.newOperator(params, new ChainWithClosure(new CopyChannelsClosure()))

        return this
    }

    private createListener() {

        final stopOnFirst = source instanceof DataflowExpression
        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( !stopOnFirst ) return
                // -- terminate the process
                processor.terminate()
                // -- close the output channels
                for( def it : outputs ) {
                    if( !(it instanceof DataflowExpression))
                        it.bind(Channel.STOP)

                    else if( !(it as DataflowExpression).isBound() )
                        it.bind(Channel.STOP)

                }
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        return [listener]
    }

}
