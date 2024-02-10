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
import org.codehaus.groovy.runtime.InvokerHelper
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
/**
 * Implements Nextflow `until` operator supporting more than one source channel
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class UntilManyOp {

    private List<DataflowReadChannel> sources
    private final Closure<Boolean> closure

    UntilManyOp(List<DataflowReadChannel> sources, final Closure<Boolean> closure) {
        this.sources = sources
        this.closure = closure
    }

    List<DataflowWriteChannel> apply() {
        final len = sources.size()
        final targets = new ArrayList(len)
        for( DataflowReadChannel it : sources )
            targets << CH.createBy(it)

        newOperator(sources, targets, new UntilCondition(len, closure))

        return targets
    }

    static class UntilCondition extends Closure {
        private int numOfChannels
        private Closure<Boolean> condition

        UntilCondition(int len, Closure<Boolean> condition) {
            super(null, null)
            this.numOfChannels = len
            this.condition = condition
            final numOfParams = condition.maximumNumberOfParameters
            if( (numOfChannels==1 && numOfParams>1) || ( numOfChannels>1 && numOfChannels!=numOfParams && numOfParams>1 ) )
                throw new IllegalArgumentException("Number of arguments mismatch")
        }

        @Override
        int getMaximumNumberOfParameters() {
            return numOfChannels
        }

        @Override
        Class[] getParameterTypes() {
            def result = new Class[numOfChannels]
            for( int i=0; i<result.size(); i++ )
                result[i] = Object
            return result
        }

        def spread(Object xyz) {
            return xyz
        }

        @Override
        Object call(final Object... args) {
            final single = args.size()>1 && condition.maximumNumberOfParameters==1
            final value = single
                        ? condition.call(args.toList())
                        : condition.call(*args) // not clear why it requires the use of the spread operator to invoke the `call(Object[])` method signature
            final result = DefaultTypeTransformation.castToBoolean(value)
            final proc = ((DataflowProcessor) getDelegate())

            if( result ) {
                proc.bindOutput(Channel.STOP)
                proc.terminate()
            }
            else {
                proc.bindAllOutputValues(InvokerHelper.asArray(args))
            }

            return null
        }

        @Override
        Object call(final Object arguments) {
            throw new UnsupportedOperationException()
        }

        @Override
        Object call() {
            throw new UnsupportedOperationException()
        }
    }
}
