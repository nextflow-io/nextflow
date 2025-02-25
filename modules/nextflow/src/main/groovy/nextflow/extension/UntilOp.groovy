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

import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
/**
 * Implements Nextflow `until` operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class UntilOp {

    private DataflowReadChannel source
    private final Closure<Boolean> closure

    UntilOp(DataflowReadChannel source, final Closure<Boolean> closure) {
        this.source = source
        this.closure = closure
    }

    DataflowWriteChannel apply() {
        final target = CH.createBy(source)
        
        newOperator(source, target, {
            final result = DefaultTypeTransformation.castToBoolean(closure.call(it))
            final proc = ((DataflowProcessor) getDelegate())

            if( result ) {
                proc.bindOutput(Channel.STOP)
                proc.terminate()
            }
            else {
                proc.bindOutput(it)
            }
        })

        return target
    }

}
