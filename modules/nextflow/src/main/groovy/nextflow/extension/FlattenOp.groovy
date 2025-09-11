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
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.Op
/**
 * Implements "flatten" operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FlattenOp {

    private static Session getSession() { Global.getSession() as Session }

    private DataflowReadChannel source
    private DataflowWriteChannel target

    FlattenOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    FlattenOp setTarget( DataflowWriteChannel target ) {
        this.target = target
        return this
    }

    DataflowWriteChannel apply() {
        final target = CH.create()
        final stopOnFirst = source instanceof DataflowExpression

        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(final DataflowProcessor dp, final List<Object> messages) {
                if( stopOnFirst )
                    dp.terminate()
            }

            @Override
            void afterStop(final DataflowProcessor dp) {
                Op.bind(dp, target, Channel.STOP)
            }

            boolean onException(final DataflowProcessor dp, final Throwable e) {
                FlattenOp.log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        new Op()
            .withInput(source)
            .withListener(listener)
            .withCode { Object item ->
                final dp = getDelegate() as DataflowProcessor
                switch( item ) {
                    case Collection:
                        ((Collection)item).flatten().each { value -> Op.bind(dp, target, value) }
                        break

                    case (Object[]):
                        ((Collection)item).flatten().each { value -> Op.bind(dp, target, value) }
                        break

                    case Channel.VOID:
                        break

                    default:
                        Op.bind(dp, target, item)
                }
            }
            .apply()

        return target
    }

}
