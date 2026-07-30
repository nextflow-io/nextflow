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

        final target = CH.createBy(left)
        final onNext = { r ->
            if( r !instanceof RecordMap )
                throw new ScriptRuntimeException("Operator `combine` with named arguments expected a channel of records but received: ${r} [${r.class.simpleName}]")

            final fields = new HashMap<String,Object>(namedArgs.size())
            for( final name : namedArgs.keySet() ) {
                final value = namedArgs[name]
                final rawValue = value instanceof ValueImpl
                    ? value.getSource().get()
                    : value
                fields.put(name, rawValue)
            }

            target << ((RecordMap) r).plus(new RecordMap(fields))
        }
        final onComplete = {
            if( !CH.isValue(target) )
                target << CH.stop()
        }
        DataflowHelper.subscribeImpl(left, [onNext: onNext, onComplete: onComplete])
        return target
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

    private Collection<?> leftValues = []
    private Collection<?> rightValues = []
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
