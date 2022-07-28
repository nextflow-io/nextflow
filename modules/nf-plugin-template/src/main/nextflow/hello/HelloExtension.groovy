package nextflow.hello

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.extension.ChannelExtensionPoint
import nextflow.extension.CH
import nextflow.NF
import nextflow.extension.DataflowHelper
import nextflow.plugin.Scoped

import java.util.concurrent.CompletableFuture

/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Slf4j
@Scoped('hello')
class HelloExtension extends ChannelExtensionPoint{

    /*
     * A session hold information about current execution of the script
     */
    private Session session

    /*
     * nf-core initializes the plugin once loaded and session is ready
     * @param session
     */
    @Override
    protected void init(Session session) {
        this.session = session
    }

    /*
     * reverse is a `producer` method and will be available to the script because:
     *
     * - it's public
     * - it returns a DataflowWriteChannel
     *
     * nf-core will inspect the extension class and allow the script to call all these kind of methods
     *
     * the method can require arguments but it's not mandatory, it depends of the business logic of the method
     *
     * business logic can write into the channel once ready and values will be consumed from it
     */
    DataflowWriteChannel reverse(String message) {
        createReverseChannel(message)
    }

    static String goodbyeMessage

    /*
    * goodbye is a `consumer` method as it receives values from a channel to perform some logic.
    *
    * Consumer methods are introspected by nextflow-core and include into the DSL if the method:
    *
    * - it's public
    * - it returns a DataflowWriteChannel
    * - it has only one arguments of DataflowReadChannel class
    *
    * a consumer method needs to proporcionate 2 closures:
    * - a closure to consume items (one by one)
    * - a finalizer closure
    *
    * in this case `goodbye` will consume a message and will store it as an upper case
    */
    DataflowWriteChannel goodbye(DataflowReadChannel source) {
        final target = CH.createBy(source)
        final next = {
            goodbyeMessage = "$it".toString().toUpperCase()
            target.bind(it)
        }
        final done = {
            target.bind(Channel.STOP)
        }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        target
    }

    protected DataflowWriteChannel createReverseChannel(final String message){
        final channel = CH.create()
        if( NF.isDsl2() ){
            session.addIgniter { ->
                businessLogicHere(channel, message)
            }
        }else{
            businessLogicHere(channel, message)
        }
        channel
    }

    /*
    * businessLogicHere will send, across the channel, the message reversed
    * and after will send an STOP signal to let know the channel it has been finished
    */
    protected static businessLogicHere(final DataflowWriteChannel channel, final String message){
        def future = CompletableFuture.runAsync({
            channel.bind(message.reverse())
            channel.bind(Channel.STOP)
        })
        future.exceptionally(this.&handlerException)
    }

    /*
    * an util class to trace exceptions
    */
    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }
}
