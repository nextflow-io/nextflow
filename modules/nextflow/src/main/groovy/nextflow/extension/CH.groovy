/*
 * Copyright 2013-2026, Seqera Labs
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
            return readChannelFromQueue(channel)

        if (channel instanceof DataflowBroadcast)
            return readChannelFromBroadcast(channel)

        if (channel instanceof DataflowReadChannel)
            return channel

        throw new IllegalArgumentException("Illegal channel source type: ${channel?.getClass()?.getName()}")
    }

    static synchronized private DataflowReadChannel readChannelFromQueue(DataflowQueue queue) {
        def broadcast = bridges.get(queue)
        if( broadcast == null ) {
            broadcast = new DataflowBroadcast()
            bridges.put(queue, broadcast)
        }
        return broadcast.createReadChannel()
    }

    static private DataflowReadChannel readChannelFromBroadcast(DataflowBroadcast channel) {
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
                    ch.add(Channel.empty())
                // map write channels to read channels
                final sources = ch.collect(it -> getReadChannel(it))
                // mix all of them
                new MixOp(sources).withTarget(topic.broadcaster).apply()
            }
            else {
                topic.broadcaster.bind(Channel.STOP)
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

        return new DataflowBroadcast()
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
        session().addIgniter {
            ch.bind(value)
        }
        return ch
    }

    static <T extends DataflowWriteChannel> T emitValues(T ch, Collection items) {
        session().addIgniter {->
            for( final it : items ) ch.bind(it)
        }
        return ch
    }

    static <T extends DataflowWriteChannel> T emitAndClose(T ch, Collection items) {
        def values = items!=null ? new ArrayList(items) : new ArrayList<>(1)
        values.add(Channel.STOP)
        emitValues(ch, values)
    }
}
