/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel

/**
 * Implements the {@link DataflowExtensions#concat} operator
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


    DataflowQueue apply() {
        final result = new DataflowQueue()
        final allChannels = [source]
        allChannels.addAll(target)

        append(result, allChannels, 0)
        return result
    }


    private static void append( DataflowWriteChannel result, List<DataflowReadChannel> channels, int index ) {
        def current = channels[index++]
        def next = index < channels.size() ? channels[index] : null

        DataflowHelper.subscribeImpl(current, [
                onNext: { result.bind(it) },
                onComplete: {
                    if(next) append(result, channels, index)
                    else result.bind(Channel.STOP)
                }
        ])
    }
}
