/*
 * Copyright 2013-2023, Seqera Labs
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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
/**
 * Operator for combining many source channels into a single channel,
 * with the option to only merge channels that are not marked as "iterators".
 *
 * @see ProcessDef#collectInputs(Object[])
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class CombineManyOp {

    private List<DataflowReadChannel> sources

    private List<Integer> iterators

    private boolean singleton

    private List<List> queues = []

    private transient List<List> combinations

    CombineManyOp(List<DataflowReadChannel> sources, List<Integer> iterators) {
        this.sources = sources
        this.iterators = iterators
        this.singleton = iterators.size() == 0 && sources.every(ch -> !CH.isChannelQueue(ch))
        this.queues = sources.collect( it -> [] )
    }

    private Map handler(int index, DataflowWriteChannel target, AtomicInteger counter) {
        final opts = new LinkedHashMap(2)
        opts.onNext = {
            onNext(target, index, it)
        }
        opts.onComplete = {
            if( counter.decrementAndGet() == 0 && !singleton )
                target.bind(Channel.STOP)
        }
        return opts
    }

    private synchronized void onNext(DataflowWriteChannel target, int index, Object value) {
        queues[index].add(value)

        // wait until every source has a value
        if( queues.any(q -> q.size() == 0) )
            return

        // emit the next item if there are no iterators
        if( iterators.size() == 0 ) {
            final args = queues.collect(q -> q.pop())
            target.bind(args)
            return
        }

        // otherwise emit an item for every iterator combination
        if( combinations == null )
            combinations = iterators.collect( i -> queues[i].first() ).combinations()

        final args = (0..<queues.size()).collect( i -> i in iterators ? null : queues[i].pop() )
        for( List entries : combinations ) {
            for( int k = 0; k < entries.size(); k++ )
                args[iterators[k]] = entries[k]

            target.bind(args)
        }
    }

    DataflowWriteChannel apply() {
        final target = CH.create(singleton)
        final counter = new AtomicInteger(sources.size())
        for( int i = 0; i < sources.size(); i++ )
            DataflowHelper.subscribeImpl( sources[i], handler(i, target, counter) )

        return target
    }
}
