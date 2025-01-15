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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.ContextRunPerThread
import nextflow.extension.op.Op
import nextflow.extension.op.OpContext
import nextflow.extension.op.OpDatum
/**
 * Implement "collate" operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CollateOp {
    private DataflowReadChannel source
    private DataflowWriteChannel target
    private int size
    private int step
    private boolean remainder
    private OpContext context = new ContextRunPerThread()

    private static Session getSession() { Global.getSession() as Session }

    CollateOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    CollateOp withTarget( DataflowWriteChannel target ) {
        this.target = target
        return this
    }

    CollateOp withSize(int size) {
        this.size = size
        return this
    }

    CollateOp withStep(int step) {
        this.step = step
        return this
    }

    CollateOp withRemainder(boolean keepRemainder) {
        this.remainder = keepRemainder
        return this
    }

    DataflowWriteChannel apply() {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        if( step <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'step' for operator 'collate' -- it must be greater than zero: $step")
        }

        // the result queue
        final target = CH.create()

        // the list holding temporary collected elements
        List<List<OpDatum>> allBuffers = []

        // -- intercepts the PoisonPill and sent out the items remaining in the buffer when the 'remainder' flag is true
        final listener = new DataflowEventAdapter() {
            Object controlMessageArrived(final DataflowProcessor dp, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if( message instanceof PoisonPill && remainder && allBuffers.size() ) {
                    for(List<OpDatum> it : allBuffers) {
                        Op.bindRunValues(target, it, false)
                    }
                }
                return message
            }

            @Override
            boolean onException(DataflowProcessor dp, Throwable e) {
                CollateOp.log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }

        int index = 0
        new Op()
            .withInput(source)
            .withOutput(target)
            .withContext(context)
            .withListener(listener)
            .withCode {

                if( index++ % step == 0 ) {
                    allBuffers.add( [] )
                }

                final run = context.getOperatorRun()
                for( List<OpDatum> list : allBuffers ) {
                    list.add(OpDatum.of(it,run))
                }

                final buf = allBuffers.head()
                if( buf.size() == size )  {
                    Op.bindRunValues(target, buf, false)
                    allBuffers = allBuffers.tail()
                }
            }
            .apply()

        return target
    }

}
