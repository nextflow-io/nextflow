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
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
/**
 * Implements the "filter" operator logic
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FilterOp {

    private DataflowReadChannel source

    private Object criteria

    FilterOp withSource(DataflowReadChannel source) {
        this.source = source
        return this
    }

    FilterOp withCriteria(Closure<Boolean> criteria) {
        this.criteria = criteria
        return this
    }

    FilterOp withCriteria(Object criteria) {
        this.criteria = criteria
        return this
    }

    DataflowWriteChannel apply() {
        assert source!=null
        assert criteria!=null

        final discriminator = criteria !instanceof Closure
                            ?  new BooleanReturningMethodInvoker("isCase")
                            : null
        final target = CH.createBy(source)
        final stopOnFirst = source instanceof DataflowExpression
        new Op()
            .withInput(source)
            .withOutput(target)
            .withCode {
                final result = criteria instanceof Closure<Boolean>
                    ? DefaultTypeTransformation.castToBoolean(criteria.call(it))
                    : discriminator.invoke(criteria, (Object)it)
                final dp = getDelegate() as DataflowProcessor
                if( result ) {
                    Op.bind(dp, target, it)
                }
                if( stopOnFirst ) {
                    Op.bind(dp, target, Channel.STOP)
                }
            }
            .apply()

        return target
    }
}
