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

import static nextflow.extension.DataflowHelper.*

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TakeOp {

    private DataflowReadChannel source;
    private int length

    TakeOp(final DataflowReadChannel source, int n) {
        this.source = source
        this.length = n
    }

    DataflowWriteChannel apply() {

        def count = 0
        final target = CH.create()

        if( length==0 ) {
            target.bind(Channel.STOP)
            return target
        }

        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                if( ++count >= length ) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }
            }

            boolean onException(final DataflowProcessor processor, final Throwable e) {
                TakeOp.log.error("@unknown", e)
                (Global.session as Session).abort(e)
                return true;
            }
        }

        newOperator(
                inputs: [source],
                outputs: [target],
                listeners: (length > 0 ? [listener] : []),
                new ChainWithClosure(new CopyChannelsClosure()))

        return target
    }
}
