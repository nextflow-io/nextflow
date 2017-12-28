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

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel

/**
 * Implements the {@link DataflowExtensions#cross} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CrossOp {

    private DataflowReadChannel source

    private DataflowReadChannel target

    private Closure mapper = DataflowExtensions.DEFAULT_MAPPING_CLOSURE

    CrossOp(DataflowReadChannel source, DataflowReadChannel target) {
        assert source
        assert target

        this.source = source
        this.target = target
    }

    CrossOp setMapper( Closure mapper ) {
        this.mapper = mapper ?: DataflowExtensions.DEFAULT_MAPPING_CLOSURE
        return this
    }

    DataflowQueue apply() {

        def result = new DataflowQueue()
        Map<Object,Map<Integer,List>> state = [:]

        final count = 2
        final stopCount = new AtomicInteger(count)

        DataflowHelper.subscribeImpl( source, crossHandlers(state, count, 0, result, mapper, stopCount ) )
        DataflowHelper.subscribeImpl( target, crossHandlers(state, count, 1, result, mapper, stopCount ) )

        return result
    }

    static private final Map crossHandlers( Map<Object,Map<Integer,List>> buffer, int size, int index, DataflowWriteChannel target, Closure mapper, AtomicInteger stopCount ) {

        [
                onNext: {
                    synchronized (buffer) {  // phaseImpl is NOT thread safe, synchronize it !
                        while( true ) {
                            def entries = PhaseOp.phaseImpl(buffer, size, index, it, mapper, true)
                            log.trace "Cross #${target.hashCode()} ($index) > item: $it; entries: $entries "

                            if( entries ) {
                                target.bind(entries)
                                // when it is invoked on the 'left' operator channel
                                // try to invoke it one more time to consume value eventually produced and accumulated by the 'right' channel
                                if( index == 0 )
                                    continue
                            }
                            break
                        }

                    }},

                onComplete: {
                    log.trace "Cross #${target.hashCode()} ($index) > Complete"
                    if( stopCount.decrementAndGet()==0) {
                        log.trace "Cross #${target.hashCode()} ($index) > STOP"
                        target << Channel.STOP
                    }}

        ]

    }

}
