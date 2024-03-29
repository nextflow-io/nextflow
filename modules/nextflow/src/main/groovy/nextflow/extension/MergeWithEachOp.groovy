/*
 * Copyright 2013-2024, Seqera Labs
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

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
/**
 * Operator for merging many source channels into a single channel,
 * with the option to combine channels that are marked as "iterators".
 *
 * @see ProcessDef#collectInputs(Object[])
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class MergeWithEachOp {

    private List<DataflowReadChannel> sources

    private List<Integer> iterators

    /**
     * List of queues to receive values from source channels.
     */
    private List<List> queues = []

    /**
     * Mask of source channels that are singletons.
     */
    private List<Boolean> singletons

    /**
     * True when all source channels are singletons and therefore
     * the operator should emit a singleton channel.
     */
    private boolean emitSingleton

    /**
     * True when all source channels are iterators and therefore
     * the operator should simply emit the combinations.
     */
    private boolean emitCombination

    private transient List<List> combinations

    MergeWithEachOp(List<DataflowReadChannel> sources, List<Integer> iterators) {
        this.sources = sources
        this.iterators = iterators
        this.queues = sources.collect( ch -> [] )
        this.singletons = sources.collect( ch -> !CH.isChannelQueue(ch) )
        this.emitSingleton = iterators.size() == 0 && singletons.every()
        this.emitCombination = iterators.size() > 0 && singletons.every()
    }

    DataflowWriteChannel apply() {
        final target = emitSingleton
            ? new DataflowVariable()
            : new DataflowQueue()
        final counter = new AtomicInteger(sources.size())
        for( int i = 0; i < sources.size(); i++ )
            DataflowHelper.subscribeImpl( sources[i], eventsMap(i, target, counter) )

        return target
    }

    private Map eventsMap(int index, DataflowWriteChannel target, AtomicInteger counter) {
        final opts = new LinkedHashMap(2)
        opts.onNext = this.&take.curry(target, index)
        opts.onComplete = {
            if( counter.decrementAndGet() == 0 && !emitSingleton && !emitCombination )
                target.bind(Channel.STOP)
        }
        return opts
    }

    private synchronized void take(DataflowWriteChannel target, int index, Object value) {
        queues[index].add(value)

        // wait until every source has a value
        if( queues.any(q -> q.size() == 0) )
            return

        // emit singleton value if every source is a singleton
        if( emitSingleton ) {
            final args = queues.collect(q -> q.first())
            target.bind(args)
            return
        }

        // emit combinations once if every source is an iterator
        if( emitCombination ) {
            emit(target)
            target.bind(Channel.STOP)
            return
        }

        // otherwise emit as many items as are available
        while( queues.every(q -> q.size() > 0) )
            emit(target)
    }

    private void emit(DataflowWriteChannel target) {
        // emit the next item if there are no iterators
        if( iterators.size() == 0 ) {
            final args = (0..<queues.size()).collect( i ->
                singletons[i] ? queues[i].first() : queues[i].pop()
            )
            target.bind(args)
            return
        }

        // otherwise emit an item for every iterator combination
        if( combinations == null )
            combinations = iterators.collect( i -> queues[i].first() ).combinations()

        final args = (0..<queues.size()).collect( i ->
            i in iterators
                ? null
                : singletons[i] ? queues[i].first() : queues[i].pop()
        )
        for( List entries : combinations ) {
            for( int k = 0; k < entries.size(); k++ )
                args[iterators[k]] = entries[k]

            target.bind(new ArrayList(args))
        }
    }
}
