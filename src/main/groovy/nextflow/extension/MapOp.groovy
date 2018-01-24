/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
/**
 * Implements {@link DataflowExtensions#map(groovyx.gpars.dataflow.DataflowReadChannel, groovy.lang.Closure)} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MapOp {

    private DataflowReadChannel source

    private Closure closure

    MapOp( final DataflowReadChannel<?> source, final Closure closure ) {
        this.source = source
        this.closure = closure
    }

    DataflowChannel apply() {

        final stopOnFirst = source instanceof DataflowExpression
        final target = DataflowExtensions.newChannelBy(source)
        DataflowHelper.newOperator(source, target) { it ->

            def result = closure.call(it)
            def proc = (DataflowProcessor) getDelegate()

            // bind the result value
            if (result != Channel.VOID)
                proc.bindOutput(result)

            // when the `map` operator is applied to a dataflow flow variable
            // terminate the processor after the first emission -- Issue #44
            if( result == Channel.STOP || stopOnFirst )
                proc.terminate()

        }

        return target
    }
}
