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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.extension.DataflowHelper
/**
 * Implements the {@link OperatorImpl#topic} operator
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class IntoTopicOp {

    private DataflowReadChannel source

    private String name

    private Session session = Global.session as Session

    IntoTopicOp( DataflowReadChannel source, String name ) {
        this.source = source
        this.name = name
    }

    void apply() {
        final target = CH.createTopicSource(name)
        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( source !instanceof DataflowExpression )
                    return
                // -- terminate the process
                processor.terminate()
                // -- send a poison pill if needed
                if( target !instanceof DataflowExpression )
                    target.bind(Channel.STOP)
                else if( !(target as DataflowExpression).isBound() )
                    target.bind(Channel.STOP)
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }

        final params = [
            inputs: List.of(source),
            outputs: List.of(target),
            listeners: List.of(listener)
        ]
        DataflowHelper.newOperator(params, new ChainWithClosure(new CopyChannelsClosure()))
    }

}
