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

import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.extension.op.Op
/**
 * Implements the "distinct" operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DistinctOp {

    private DataflowReadChannel source
    private DataflowWriteChannel target
    private Closure comparator

    DistinctOp() {}

    DistinctOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    DistinctOp withTarget(DataflowWriteChannel target) {
        assert target!=null
        this.target = target
        return this
    }

    DistinctOp withComparator(Closure comparator) {
        assert comparator!=null
        this.comparator = comparator
        return this
    }

    DataflowWriteChannel apply() {
        assert source != null
        assert comparator != null

        if( !target )
            target = CH.createBy(source)

        final params = new DataflowHelper.OpParams()
            .withInput(source)
            .withOutput(target)

        def previous = null
        DataflowHelper.newOperator(params) {
            final proc = ((DataflowProcessor) getDelegate())
            def key = comparator.call(it)
            if( key != previous ) {
                previous = key
                Op.bind(proc, target, it)
            }
        }

        return target
    }

}
