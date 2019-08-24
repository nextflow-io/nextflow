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

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
/**
 * Implements the {@link OperatorEx#cross} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CrossOp {

    private DataflowReadChannel source

    private DataflowReadChannel target

    private Closure mapper = OperatorEx.DEFAULT_MAPPING_CLOSURE

    CrossOp(DataflowReadChannel source, DataflowReadChannel target) {
        assert source
        assert target

        this.source = source
        this.target = target
    }

    CrossOp setMapper( Closure mapper ) {
        this.mapper = mapper ?: OperatorEx.DEFAULT_MAPPING_CLOSURE
        return this
    }

    DataflowWriteChannel apply() {

        final result = CH.create()
        Map<Object,Map<Integer,List>> state = [:]

        final count = 2
        final stopCount = new AtomicInteger(count)

        DataflowHelper.subscribeImpl( source, crossHandlers(state, count, 0, result, mapper, stopCount ) )
        DataflowHelper.subscribeImpl( target, crossHandlers(state, count, 1, result, mapper, stopCount ) )

        return result
    }

    static private final Map<String,Closure> crossHandlers( Map<Object,Map<Integer,List>> buffer, int size, int index, DataflowWriteChannel target, Closure mapper, AtomicInteger stopCount ) {

        def result = new HashMap<String,Closure>(2)

        result.onNext = {
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

            }}

        result.onComplete = {
            log.trace "Cross #${target.hashCode()} ($index) > Complete"
            if( stopCount.decrementAndGet()==0) {
                log.trace "Cross #${target.hashCode()} ($index) > STOP"
                target << Channel.STOP
            }}

        return result
    }

}
