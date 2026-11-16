/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.dataflow.ops

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.dataflow.ChannelImpl
import nextflow.dataflow.ValueImpl
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.util.ArrayTuple
import nextflow.util.RecordMap
/**
 * Implements the `combine` operator for typed workflows.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class CombineOpV2 {

    private DataflowReadChannel left

    private DataflowReadChannel right

    private Map<String,Object> namedArgs

    CombineOpV2(DataflowReadChannel left, Object other) {
        this.left = left
        if( other instanceof DataflowReadChannel )
            this.right = other
        else if( other instanceof Map )
            this.namedArgs = (Map<String,Object>) other
        else
            throw new IllegalArgumentException()
    }

    DataflowWriteChannel apply() {
        if( namedArgs )
            return applyNamedArgs()
        final target = CH.createBy(left)
        DataflowHelper.subscribeImpl( left,  handler(target, 0) )
        DataflowHelper.subscribeImpl( right, handler(target, 1) )
        return target
    }

    private DataflowWriteChannel applyNamedArgs() {
        for( final name : namedArgs.keySet() ) {
            if( namedArgs[name] instanceof ChannelImpl )
                throw new ScriptRuntimeException("Operator `combine` named argument '${name}' cannot be a channel")
        }

        // the named dataflow values are operator inputs, read on every firing, so that the
        // operator never blocks a dataflow thread waiting for them
        final inputs = new ArrayList<DataflowReadChannel>()
        inputs.add(left)
        for( final value : namedArgs.values() ) {
            if( value instanceof ValueImpl )
                inputs.add(value.getSource())
        }

        final target = CH.createBy(left)
        final params = new HashMap<String,Object>()
        params.inputs = inputs
        params.outputs = [target]
        params.listeners = [DataflowHelper.stopErrorListener(left, target)]
        DataflowHelper.newOperator(params, new NamedArgsClosure(namedArgs, inputs.size()))
        return target
    }

    /**
     * Adds the named arguments to each record, taking the dataflow values from the operator inputs
     */
    @CompileStatic
    private static class NamedArgsClosure extends Closure {

        private Map<String,Object> namedArgs

        private int numOfParams

        NamedArgsClosure(Map<String,Object> namedArgs, int numOfParams) {
            super(null, null)
            this.namedArgs = namedArgs
            this.numOfParams = numOfParams
        }

        @Override
        int getMaximumNumberOfParameters() {
            return numOfParams
        }

        @Override
        Class[] getParameterTypes() {
            return Collections.nCopies(numOfParams, Object) as Class[]
        }

        @Override
        Object call(Object... args) {
            final r = args[0]
            if( r !instanceof RecordMap )
                throw new ScriptRuntimeException("Operator `combine` with named arguments expected a channel of records but received: ${r} [${r.class.simpleName}]")

            final fields = new HashMap<String,Object>(namedArgs.size())
            int i = 1
            for( final name : namedArgs.keySet() ) {
                final value = namedArgs[name]
                fields.put(name, value instanceof ValueImpl ? args[i++] : value)
            }

            final result = ((RecordMap) r).plus(new RecordMap(fields))
            ((DataflowProcessor) getDelegate()).bindOutput(result)
            return result
        }
    }

    private Map<String,Closure> handler(DataflowWriteChannel target, int index) {
        final result = new HashMap<String,Closure>(2)
        result.onNext = { value ->
            onNext(target, value, index)
        }
        result.onComplete = {
            onComplete(target, index)
        }
        return result
    }

    private Collection<Object> leftValues = []
    private Collection<Object> rightValues = []
    private int count = 2

    private synchronized void onNext(DataflowWriteChannel target, Object value, int index) {
        if( index == 0 ) {
            final leftValue = value
            for( final rightValue : rightValues )
                target << combine0(leftValue, rightValue)
            leftValues.add(leftValue)
        }
        else if( index == 1 ) {
            final rightValue = value
            for( final leftValue : leftValues )
                target << combine0(leftValue, rightValue)
            rightValues.add(rightValue)
        }
    }

    private static ArrayTuple combine0(Object left, Object right) {
        final values = []
        append0(values, left)
        append0(values, right)
        return new ArrayTuple<>(values)
    }

    private static void append0(List values, Object value) {
        if( value instanceof ArrayTuple )
            values.addAll(value)
        else
            values.add(value)
    }

    private synchronized void onComplete(DataflowWriteChannel target, int index) {
        count--
        if( count > 0 )
            return

        if( !CH.isValue(target) )
            target << CH.stop()
    }

}
