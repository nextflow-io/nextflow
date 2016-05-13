/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
    DataflowWriteChannel target

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
        this.target = DataflowExtensions.newChannelBy(source)

        // -- set the target variable in the script binding context
        final name = CaptureProperties.capture(holder)
        if( !name )
            throw new IllegalArgumentException("Missing target channel on `tap` operator")

        if( name.size()>1 )
            throw new IllegalArgumentException("Operator `tap` does not allow more than one target channel")

        final binding = Global.session.binding
        if( binding.hasVariable(name[0]) ) {
            log.warn "A variable named '${name[0]}' already exists in script global context -- Consider renaming it "
        }
        binding.setVariable(name[0], target)

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
        this.target = target
        this.result = DataflowExtensions.newChannelBy(source)
    }

    /**
     * @return A list holding the output channels of the `tap` operator
     */
    List<DataflowWriteChannel> getOutputs() { [result, target] }

    /**
     * Apply the operator
     * @return An instance of {@link TapOp} itself
     */
    TapOp apply() {
        DataflowExtensions.newOperator([source], [result, target], new ChainWithClosure(new CopyChannelsClosure()));
        return this
    }


}
