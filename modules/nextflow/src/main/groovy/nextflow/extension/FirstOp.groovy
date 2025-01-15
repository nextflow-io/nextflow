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
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.Op
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FirstOp {

    private DataflowReadChannel source
    private DataflowVariable target
    private Object criteria

    FirstOp() {}

    FirstOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    FirstOp withTarget(DataflowVariable target) {
        assert target!=null
        this.target = target
        return this
    }

    FirstOp withCriteria(Object criteria) {
        assert criteria!=null
        this.criteria = criteria
        return this
    }

    DataflowWriteChannel apply() {
        assert source!=null
        assert criteria!=null

        if( target==null )
            target = new DataflowVariable()

        final stopOnFirst = source instanceof DataflowExpression
        final discriminator = new BooleanReturningMethodInvoker("isCase");

        final listener = new DataflowEventAdapter() {
            @Override
            void afterStop(DataflowProcessor dp) {
                if( stopOnFirst && !target.isBound() )
                    Op.bind(dp, target, Channel.STOP)
            }
        }

        final code = {
            final dp = getDelegate() as DataflowProcessor
            final accept = discriminator.invoke(criteria, it)
            if( accept )
                Op.bind(dp, target, it)
            if( accept || stopOnFirst )
                dp.terminate()
        }

        new Op()
            .withInput(source)
            .withListener(listener)
            .withCode(code)
            .apply()

        return target
    }
}
