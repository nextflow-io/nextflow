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

import static nextflow.extension.DataflowHelper.*

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op
/**
 * Implement "flatMap" operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FlatMapOp {

    private DataflowReadChannel source
    private DataflowWriteChannel target
    private Closure mapper

    FlatMapOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    FlatMapOp withMapper(Closure code) {
        this.mapper = code
        return this
    }

    FlatMapOp setTarget( DataflowWriteChannel target ) {
        this.target = target
        return this
    }

    DataflowWriteChannel apply() {
        assert source!=null

        if( target==null )
            target = CH.create()

        new Op()
            .withInput(source)
            .withOutput(target)
            .withListener(stopErrorListener(source,target))
            .withCode { Object item ->
                    final result = mapper != null ? mapper.call(item) : item
                    final dp = getDelegate() as DataflowProcessor

                    switch( result ) {
                        case Collection:
                            result.each { it -> Op.bind(dp, target,it) }
                            break

                        case (Object[]):
                            result.each { it -> Op.bind(dp, target,it) }
                            break

                        case Map:
                            result.each { it -> Op.bind(dp, target,it) }
                            break

                        case Map.Entry:
                            Op.bind(dp, target, (result as Map.Entry).key )
                            Op.bind(dp, target, (result as Map.Entry).value )
                            break

                        case Channel.VOID:
                            break

                        default:
                            Op.bind(dp, target, result)
                    }
                }
            .apply()
        return target
    }

}
