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
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.Op
/**
 * Implement "take" operator
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
            void afterRun(final DataflowProcessor dp, final List<Object> messages) {
                if( ++count >= length ) {
                    dp.bindOutput( Channel.STOP )
                    dp.terminate()
                }
            }

            boolean onException(final DataflowProcessor dp, final Throwable e) {
                TakeOp.log.error("@unknown", e)
                (Global.session as Session).abort(e)
                return true;
            }
        }

        new Op()
            .withInput(source)
            .withOutput(target)
            .withListener(length>0 ? listener : null)
            .withCode {
                final proc = getDelegate() as DataflowProcessor
                Op.bind(proc, target, it)
            }
            .apply()

        return target
    }
}
