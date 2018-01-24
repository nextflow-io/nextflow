/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
