package nextflow.extension

import static nextflow.Channel.*

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

    static class Topic {
        DataflowBroadcast broadcaster = new DataflowBroadcast()
        List<DataflowWriteChannel> sources = new ArrayList<>(10)
    }

    static final private Map<String, Topic> allTopics = new HashMap<>(10)

    static final private Map<DataflowQueue, DataflowBroadcast> bridges = new HashMap<>(10)

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
        channel.createReadChannel()
    }

    static synchronized boolean isBridge(DataflowQueue queue) {
        bridges.get(queue) != null
    }

    static void broadcast() {
        // connect all broadcast topics, note this must be before the following
        // "bridging" step because it can modify the final network topology
        connectTopics()
        // bridge together all broadcast channels
        bridgeChannels()
    }

    static private void bridgeChannels() {
        // connect all dataflow queue variables to associated broadcast channel
        for( DataflowQueue queue : bridges.keySet() ) {
            log.trace "Bridging dataflow queue=$queue"
            def broadcast = bridges.get(queue)
            queue.into(broadcast)
        }
    }

    static private void connectTopics() {
        for( Topic topic : allTopics.values() ) {
            if( topic.sources ) {
                // get the list of source channels for this topic
                final ch = new ArrayList(topic.sources)
                // the mix operator requires at least two sources, add an empty channel if needed
                if( ch.size()==1 )
                    ch.add(empty())
                // map write channels to read channels
                final sources = ch.collect(it -> getReadChannel(it))
                // mix all of them 
                new MixOp(sources).withTarget(topic.broadcaster).apply()
            }
            else {
                topic.broadcaster.bind(STOP)
            }
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

    static DataflowBroadcast topic(String name) {
        synchronized (allTopics) {
            def topic = allTopics[name]
            if( topic!=null )
                return topic.broadcaster
            // create a new topic
            topic = new Topic()
            allTopics[name] = topic
            return topic.broadcaster
        }
    }

    static DataflowWriteChannel createTopicSource(String name) {
        synchronized (allTopics) {
            def topic = allTopics.get(name)
            if( topic==null ) {
                topic = new Topic()
                allTopics.put(name, topic)
            }
            final result = CH.create()
            topic.sources.add(result)
            return result
        }
    }

    static boolean isChannel(obj) {
        obj instanceof DataflowReadChannel || obj instanceof DataflowWriteChannel
    }

    static boolean isValue(obj) {
        return obj instanceof DataflowExpression
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
