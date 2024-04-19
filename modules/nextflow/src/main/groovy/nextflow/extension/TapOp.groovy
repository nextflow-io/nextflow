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
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.NF
import nextflow.Session
/**
 * Implements the {@link OperatorImpl#tap} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TapOp {

    private DataflowReadChannel source

    private List<DataflowWriteChannel> outputs

    private DataflowWriteChannel result

    private Session session = (Session)Global.session

    /**
     * Create the operator instance
     *
     * @param source An instance of {@link DataflowReadChannel} used to feed the operator
     * @param outputs A list of {@link DataflowWriteChannel}'s in which to 'tap' the source channel
     */
    TapOp( DataflowReadChannel source, List<DataflowWriteChannel> outputs ) {
        this.source = source
        this.result = CH.createBy(source)
        this.outputs = [result, *outputs]
    }

    /**
     * Create the operator instance
     *
     * @param source An instance of {@link DataflowReadChannel} used to feed the operator
     * @param holder A closure used to declare the target channel name e.g. {@code { ch_foo ; ch_bar } }
     */
    TapOp( DataflowReadChannel source, Closure holder ) {
        this(source, [])

        // -- set the target variable in the script binding context
        final names = CaptureProperties.capture(holder)
        if( !names )
            throw new IllegalArgumentException("Missing target channel on `tap` operator")

        final binding = NF.binding
        names.each { name ->
            final channel = CH.createBy(source)
            if( binding.hasVariable(name) )
                log.warn "A variable named '${name}' already exists in the script global context -- Consider renaming it "

            binding.setVariable(name, channel)
            outputs << channel
        }

    }

    List<DataflowWriteChannel> getOutputs() { outputs }

    DataflowWriteChannel getResult() { result }

    TapOp apply() {

        final params = [:]
        params.inputs = [source]
        params.outputs = outputs
        params.listeners = createListener()

        DataflowHelper.newOperator(params, new ChainWithClosure(new CopyChannelsClosure()))
        return this
    }

    private createListener() {

        final stopOnFirst = source instanceof DataflowExpression
        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( !stopOnFirst )
                    return
                // -- terminate the process
                processor.terminate()
                // -- close the output channels
                for( final output : outputs ) {
                    if( output !instanceof DataflowExpression )
                        output.bind(Channel.STOP)

                    else if( !(output as DataflowExpression).isBound() )
                        output.bind(Channel.STOP)
                }
            }

            @Override
            boolean onException(DataflowProcessor processor, Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }

        return [listener]
    }

}
