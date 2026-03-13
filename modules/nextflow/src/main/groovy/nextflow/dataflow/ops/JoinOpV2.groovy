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
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.types.Record
import nextflow.util.CheckHelper
import nextflow.util.RecordMap
/**
 * Implements the `join` operator for typed dataflow.
 *
 * This operator is specialized for records. It expects records
 * from the source channels and emits records to the output channel.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class JoinOpV2 {

    private static final Map JOIN_PARAMS = [by: String, remainder: Boolean]

    private DataflowReadChannel left

    private DataflowReadChannel right

    private String pivot

    private boolean remainder

    JoinOpV2(DataflowReadChannel left, DataflowReadChannel right, Map opts = [:]) {
        if( opts.by == null )
            throw new ScriptRuntimeException("Operator `join` requires the `by` option")
        CheckHelper.checkParams('join', opts, JOIN_PARAMS)
        this.left = left
        this.right = right
        this.pivot = opts.by as String
        this.remainder = opts.remainder ? opts.remainder as boolean : false
    }

    DataflowWriteChannel apply() {
        final target = CH.create()
        DataflowHelper.subscribeImpl( left,  handler(target, 0) )
        DataflowHelper.subscribeImpl( right, handler(target, 1) )
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

    private Map<?,Collection<RecordMap>> mappingsLeft = [:]
    private Map<?,Collection<RecordMap>> mappingsRight = [:]
    private Set<?> emitted = []
    private int count = 2

    private synchronized void onNext(DataflowWriteChannel target, Object value, int index) {
        if( value !instanceof RecordMap )
            throw new ScriptRuntimeException("Operator `join` expected a record but received: ${value} [${value.class.simpleName}]")
        if( index == 0 ) {
            final leftValue = value as RecordMap
            final key = leftValue[pivot]
            final leftValues = mappingsLeft.computeIfAbsent(key, (k) -> new ArrayList<>())
            final rightValues = mappingsRight.computeIfAbsent(key, (k) -> new ArrayList<>())
            for( final rightValue : rightValues )
                target << leftValue.plus(rightValue)
            if( !rightValues.isEmpty() )
                emitted.add(key)
            leftValues.add(leftValue)
        }
        else if( index == 1 ) {
            final rightValue = value as RecordMap
            final key = rightValue[pivot]
            final leftValues = mappingsLeft.computeIfAbsent(key, (k) -> new ArrayList<>())
            final rightValues = mappingsRight.computeIfAbsent(key, (k) -> new ArrayList<>())
            for( final leftValue : leftValues )
                target << leftValue.plus(rightValue)
            if( !leftValues.isEmpty() )
                emitted.add(key)
            rightValues.add(rightValue)
        }
    }

    private synchronized void onComplete(DataflowWriteChannel target, int index) {
        count--
        if( count > 0 )
            return

        if( remainder ) {
            // emit unmatched left-hand values
            final keysLeft = mappingsLeft.keySet()
            keysLeft.removeAll(emitted)
            for( final key : keysLeft ) {
                for( final value : mappingsLeft.get(key) )
                    target << value
            }

            // emit unmatched right-hand values
            final keysRight = mappingsRight.keySet()
            keysRight.removeAll(emitted)
            for( final key : keysRight ) {
                for( final value : mappingsRight.get(key) )
                    target << value
            }
        }

        target << CH.stop()
    }

}
