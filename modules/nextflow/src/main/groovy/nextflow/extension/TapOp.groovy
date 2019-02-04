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
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import nextflow.Global
/**
 * Implements the {@link DataflowExtensions#tap} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TapOp {

    /**
     * Operator input channel
     */
    DataflowReadChannel source

    /**
     * Operator `tapped` channel
     */
    List<DataflowWriteChannel> outputs

    /**
     * Operator output channel
     */
    DataflowWriteChannel result

    /**
     * Create the operator instance
     *
     * @param source An instance of {@link DataflowReadChannel} used to feed the operator
     * @param holder A closure used to declare the target channel name e.g. {@code { targetChannelName } }
     */
    TapOp( DataflowReadChannel source, Closure holder ) {
        assert source != null
        assert holder != null

        this.source = source
        this.result = DataflowExtensions.newChannelBy(source)
        this.outputs = [result]

        // -- set the target variable in the script binding context
        final names = CaptureProperties.capture(holder)
        if( !names )
            throw new IllegalArgumentException("Missing target channel on `tap` operator")

        final binding = Global.session.binding
        names.each { item ->
            def channel = DataflowExtensions.newChannelBy(source)
            if( binding.hasVariable(item) )
                log.warn "A variable named '${item}' already exists in script global context -- Consider renaming it "

            binding.setVariable(item, channel)
            outputs << channel
        }

    }

    /**
     * Create the operator instance
     *
     * @param source An instance of {@link DataflowReadChannel} used to feed the operator
     * @param target An instance of {@link DataflowWriteChannel} that will receive the items emitted by the source
     */
    TapOp( DataflowReadChannel source, DataflowWriteChannel target ) {
        assert source != null
        assert target != null
        if( source.class != target.class ) {
                throw new IllegalArgumentException("Operator `tap` source and target channel types must match -- source: ${source.class.name}, target: ${target.class.name} ")
        }

        this.source = source
        this.result = DataflowExtensions.newChannelBy(source)
        this.outputs = [result, target]
    }

    /**
     * @return A list holding the output channels of the `tap` operator
     */
    List<DataflowWriteChannel> getOutputs() { outputs }

    /**
     * Apply the operator
     * @return An instance of {@link TapOp} itself
     */
    TapOp apply() {
        DataflowHelper.newOperator([source], outputs, new ChainWithClosure(new CopyChannelsClosure()));
        return this
    }


}
