/*
 * Copyright 2013-2023, Seqera Labs
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

import static nextflow.extension.DataflowHelper.createOpParams
import static nextflow.extension.DataflowHelper.newOperator
import static nextflow.extension.DataflowHelper.stopErrorListener

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ChainWithClosure

/**
 * Implements merge operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MergeOp {
    private DataflowReadChannel source
    private List<DataflowReadChannel> others
    private boolean flat
    private Closure closure

    MergeOp(DataflowReadChannel source, List<DataflowReadChannel> others, Map opts=null, Closure closure=null) {
        this.source = source
        this.others = others
        this.flat = opts?.flat!=null ? opts?.flat : true
        this.closure = closure
    }

    MergeOp(DataflowReadChannel source, DataflowReadChannel other, Map opts=null, Closure closure=null) {
        this(source, Collections.singletonList(other), opts, closure)
    }

    DataflowWriteChannel apply() {
        final result = CH.createBy(source)
        final List<DataflowReadChannel> inputs = new ArrayList<DataflowReadChannel>(1 + others.size())
        final action = closure ? new ChainWithClosure<>(closure) : new DefaultMergeClosure(1 + others.size(), flat)
        inputs.add(source)
        inputs.addAll(others)
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, action)
        return result
    }
}
