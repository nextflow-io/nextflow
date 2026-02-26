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
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op
/**
 * Implements the "into" operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IntoOp {

    private DataflowReadChannel source

    private List<DataflowWriteChannel> outputs

    IntoOp( DataflowReadChannel source, List<DataflowWriteChannel> targets ) {
        assert source
        assert targets
        this.source = source
        this.outputs = targets
    }

    IntoOp( DataflowReadChannel source, int n ) {
        assert source
        assert n

        def targets = new ArrayList(n)
        for( int i=0; i<n; i++ )
            targets << new DataflowQueue()

        this.source = source
        this.outputs = targets
    }

    List<DataflowWriteChannel> getOutputs() { outputs }

    IntoOp apply() {
        new SubscribeOp()
            .withInput(source)
            .withOnNext(this.&doNext)
            .withOnComplete(this.&doComplete)
            .apply()
        return this
    }

    private void doNext(DataflowProcessor dp, Object it) {
        for( DataflowWriteChannel ch : outputs ) {
            Op.bind(dp, ch, it)
        }
    }

    private void doComplete(DataflowProcessor dp) {
        for( DataflowWriteChannel ch : outputs ) {
            if( ch instanceof DataflowExpression ) {
                if( !ch.isBound()) Op.bind(dp, ch, Channel.STOP)
            }
            else {
                Op.bind(dp, ch, Channel.STOP)
            }
        }
    }

}
