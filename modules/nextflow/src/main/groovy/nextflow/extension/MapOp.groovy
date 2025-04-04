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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op

/**
 * Implements {@link OperatorImpl#map(groovyx.gpars.dataflow.DataflowReadChannel, groovy.lang.Closure)} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MapOp {

    private DataflowReadChannel source

    private Closure mapper

    private DataflowWriteChannel target

    MapOp() {}

    MapOp( final DataflowReadChannel<?> source, final Closure mapper ) {
        this.source = source
        this.mapper = mapper
    }

    MapOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    MapOp withMapper(Closure code) {
        assert code!=null
        this.mapper = code
        return this
    }

    MapOp setTarget( DataflowWriteChannel target ) {
        this.target = target
        return this
    }

    DataflowWriteChannel apply() {

        if( target == null )
            target = CH.createBy(source)

        final stopOnFirst = source instanceof DataflowExpression
        new Op()
            .withInput(source)
            .withOutput(target)
            .withCode { it ->
                final result = mapper.call(it)
                final proc = getDelegate() as DataflowProcessor

                // bind the result value
                if (result != Channel.VOID)
                    Op.bind(proc, target, result)

                // when the `map` operator is applied to a dataflow flow variable
                // terminate the processor after the first emission -- Issue #44
                if( stopOnFirst )
                    proc.terminate()
            }
            .apply()

        return target
    }
}
