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
 *
 */

package nextflow.extension

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.ContextRunPerThread
import nextflow.extension.op.Op
import nextflow.extension.op.OpContext

/**
 * Implements Nextflow Mix operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MixOp {

    private DataflowReadChannel source
    private List<DataflowReadChannel> others
    private DataflowWriteChannel target
    private OpContext context = new ContextRunPerThread()

    MixOp(DataflowReadChannel source, DataflowReadChannel other) {
        this.source = source
        this.others = List.of(other)
    }

    MixOp(DataflowReadChannel source, DataflowReadChannel[] others) {
        this.source = source
        this.others = others.toList()
    }

    MixOp(List<DataflowReadChannel> channels) {
        if( channels.size()<2 )
            throw new IllegalArgumentException("Mix operator requires at least 2 source channels")
        this.source = channels.get(0)
        this.others = channels.subList(1, channels.size())
    }

    MixOp withTarget(DataflowWriteChannel target) {
        this.target = target
        return this
    }

    DataflowWriteChannel apply() {
        if( target == null )
            target = CH.create()
        final count = new AtomicInteger( others.size()+1 )
        final handlers = [
                onNext: { DataflowProcessor dp, it -> Op.bind(dp, target, it) },
                onComplete: { DataflowProcessor dp -> if(count.decrementAndGet()==0) { Op.bind(dp, target, Channel.STOP) } }
        ]

        subscribe0(source, handlers)
        for( def it : others ) {
            subscribe0(it, handlers)
        }

        final allSources = [source]
        allSources.addAll(others)
        return target
    }

    private void subscribe0(final DataflowReadChannel source, final Map<String,Closure> events ) {
        new SubscribeOp()
            .withInput(source)
            .withContext(context)
            .withEvents(events)
            .apply()
    }

}
