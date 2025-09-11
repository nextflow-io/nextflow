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
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
/**
 * Implements Nextflow `until` operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class UntilOp {

    private DataflowReadChannel source
    private final Closure<Boolean> closure

    UntilOp(DataflowReadChannel source, final Closure<Boolean> closure) {
        this.source = source
        this.closure = closure
    }

    DataflowWriteChannel apply() {
        final target = CH.createBy(source)

        new Op()
            .withInput(source)
            .withOutput(target)
            .withCode {
                final result = DefaultTypeTransformation.castToBoolean(closure.call(it))
                final dp = getDelegate() as DataflowProcessor
                if( result ) {
                    Op.bind(dp, target, Channel.STOP)
                    dp.terminate()
                }
                else {
                    Op.bind(dp, target, it)
                }
            }
        .apply()

        return target
    }

}
