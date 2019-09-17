package nextflow.extension

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.stream.DataflowStreamReadAdapter
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import nextflow.Channel
import nextflow.Global
import nextflow.NF
import nextflow.Session
import static nextflow.Channel.STOP

/**
 * Helper class to handle channel internal api ops
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CH {

    static private Session session() {
        return (Session) Global.session
    }

    static private Map<DataflowQueue, DataflowBroadcast> bridges = new HashMap<>(10)

    static DataflowReadChannel getReadChannel(channel) {
        if (channel instanceof DataflowQueue)
            return getRead1(channel)

        if (channel instanceof DataflowBroadcast)
            return getRead2(channel)

        if (channel instanceof DataflowReadChannel)
            return channel

        throw new IllegalArgumentException("Illegal channel source type: ${channel?.getClass()?.getName()}")
    }

    static synchronized private DataflowReadChannel getRead1(DataflowQueue queue) {
        if( !NF.isDsl2() )
            return queue

        def broadcast = bridges.get(queue)
        if( broadcast == null ) {
            broadcast = new DataflowBroadcast()
            bridges.put(queue, broadcast)
        }
        return broadcast.createReadChannel()
    }

    static private DataflowReadChannel getRead2(DataflowBroadcast channel) {
        if( !NF.isDsl2() )
            throw new IllegalStateException("Broadcast channel are only allowed in a workflow definition scope")
        channel.createReadChannel()
    }

    static synchronized boolean isBridge(DataflowQueue queue) {
        bridges.get(queue) != null
    }

    static void broadcast() {
        // connect all dataflow queue variables to associated broadcast channel 
        for( DataflowQueue queue : bridges.keySet() ) {
            log.debug "Bridging dataflow queue=$queue"
            def broadcast = bridges.get(queue)
            queue.into(broadcast)
        }
    }

    static void init() { bridges.clear() }

    @PackageScope
    static DataflowWriteChannel close0(DataflowWriteChannel source) {
        if( source instanceof DataflowExpression ) {
            if( !source.isBound() )
                source.bind(Channel.STOP)
        }
        else {
            source.bind(Channel.STOP)
        }
        return source
    }

    static DataflowWriteChannel createBy(DataflowReadChannel channel) {
        create( channel instanceof DataflowExpression )
    }

    static DataflowWriteChannel create(boolean value=false) {
        if( value )
            return new DataflowVariable()

        if( NF.isDsl2() )
            return new DataflowBroadcast()

        return new DataflowQueue()
    }

    static boolean isChannel(obj) {
        obj instanceof DataflowReadChannel || obj instanceof DataflowWriteChannel
    }

    static boolean isChannelQueue(obj) {
        obj instanceof DataflowQueue || obj instanceof DataflowStreamReadAdapter || obj instanceof DataflowStreamWriteAdapter
    }

    static boolean allScalar(List args) {
        for( def el : args ) {
            if( isChannelQueue(el) ) {
                return false
            }
        }
        return true
    }

    static DataflowVariable value() {
        return new DataflowVariable()
    }

    static DataflowVariable value(obj) {
        final result = new DataflowVariable()
        emit(result, obj)
        return result
    }

    static DataflowQueue queue() {
        new DataflowQueue()
    }


    static DataflowWriteChannel emit(DataflowWriteChannel ch, Object value) {
        if(NF.isDsl2()) {
            session().addIgniter { ch.bind(value) }
        }
        else {
            ch.bind(value)
        }
        return ch
    }

    static <T extends DataflowWriteChannel> T emitValues(T ch, Collection items) {
        if(NF.dsl2) {
            session().addIgniter {-> for( def it : items ) ch.bind(it) }
        }
        else {
            for( def it : items ) ch.bind(it)
        }
        return ch
    }

    static <T extends DataflowWriteChannel> T emitAndClose(T ch, Collection items) {
        def values = items!=null ? new ArrayList(items) : new ArrayList<>(1)
        values.add(STOP)
        emitValues(ch, values)
    }
}
