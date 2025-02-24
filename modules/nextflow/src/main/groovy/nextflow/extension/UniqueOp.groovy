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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.extension.op.Op

/**
 * Implements the "unique" operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class UniqueOp {

    private DataflowReadChannel source
    private DataflowWriteChannel target
    private Closure comparator

    UniqueOp() {}

    UniqueOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    UniqueOp withTarget(DataflowWriteChannel target) {
        assert target!=null
        this.target = target
        return this
    }

    UniqueOp withComparator(Closure comparator) {
        assert comparator!=null
        this.comparator = comparator
        return this
    }

    DataflowWriteChannel apply() {
        assert source != null
        assert comparator != null

        if( !target )
            target = CH.createBy(source)

        final history = new LinkedHashMap()
        final stopOnFirst = source instanceof DataflowExpression

        // when the operator stop clear the history map
        final listener = new DataflowEventAdapter() {
            void afterStop(final DataflowProcessor dp) {
                history.clear()
            }
        }

        final code = {
            final proc = getDelegate() as DataflowProcessor
            try {
                final key = comparator.call(it)
                if( !history.containsKey(key) ) {
                    history.put(key,true)
                    Op.bind(proc, target, it)
                }
            }
            finally {
                if( stopOnFirst ) {
                    proc.terminate()
                }
            }
        }

        new Op()
            .withInput(source)
            .withOutput(target)
            .withListener(listener)
            .withCode(code)
            .apply()

        return target
    }

}
