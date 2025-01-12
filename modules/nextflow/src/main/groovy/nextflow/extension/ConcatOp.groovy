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

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op

/**
 * Implements the {@link OperatorImpl#concat} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConcatOp {

    private DataflowReadChannel source

    private DataflowReadChannel[] target

    ConcatOp( DataflowReadChannel source, DataflowReadChannel... target ) {
        assert source != null
        assert target

        this.source = source
        this.target = target
    }

    DataflowWriteChannel apply() {
        final result = CH.create()
        final allChannels = [source]
        allChannels.addAll(target)

        append(result, allChannels, 0)
        return result
    }


    private static void append( DataflowWriteChannel result, List<DataflowReadChannel> channels, int index ) {
        final current = channels[index++]
        final next = index < channels.size() ? channels[index] : null

        final events = new HashMap<String,Closure>(2)
        events.onNext = { DataflowProcessor proc, it -> Op.bind(proc, result, it) }
        events.onComplete = { DataflowProcessor proc ->
            if(next) append(result, channels, index)
            else Op.bind(proc, result, Channel.STOP)
        }

        new SubscribeOp()
            .withSource(current)
            .withEvents(events)
            .apply()
    }
}
