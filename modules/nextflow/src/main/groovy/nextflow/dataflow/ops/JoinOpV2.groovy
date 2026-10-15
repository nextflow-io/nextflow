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
 * Implements the `join` operator for typed workflows.
 *
 * This operator is specialized for records. It expects records
 * from the source channels and emits records to the output channel.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class JoinOpV2 {

    private static final Map JOIN_PARAMS = [by: [String, Closure], remainder: Boolean]

    private DataflowReadChannel left

    private DataflowReadChannel right

    private String pivot

    private Closure keyFn

    private Closure predicate

    private boolean remainder

    JoinOpV2(DataflowReadChannel left, DataflowReadChannel right, Map opts = [:]) {
        if( opts.by == null )
            throw new ScriptRuntimeException("Operator `join` requires the `by` option")
        CheckHelper.checkParams('join', opts, JOIN_PARAMS)
        this.left = left
        this.right = right
        if( opts.by instanceof Closure ) {
            final closure = opts.by as Closure
            final arity = closure.maximumNumberOfParameters
            if( arity == 1 )
                this.keyFn = closure
            else if( arity == 2 )
                this.predicate = closure
            else
                throw new ScriptRuntimeException("Operator `join` expected a closure with one parameter (key function) or two parameters (join predicate) for the `by` option but received a closure with ${arity}")
        }
        else
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
            onNextSafe(target, value, index)
        }
        result.onComplete = {
            onComplete(target, index)
        }
        return result
    }

    private Map<?,Collection<RecordMap>> mappingsLeft = [:]
    private Map<?,Collection<RecordMap>> mappingsRight = [:]
    private Set<?> emitted = []

    // buffers for the predicate (closure) join
    private List<RecordMap> leftBuf = []
    private List<RecordMap> rightBuf = []
    private BitSet leftMatched = new BitSet()
    private BitSet rightMatched = new BitSet()

    private volatile int count = 2
    private volatile boolean failed

    private synchronized void onNextSafe(DataflowWriteChannel target, Object value, int index) {
        if( failed )
            return
        try {
            onNext(target, value, index)
        } catch( Throwable e ) {
            this.failed = true
            target << CH.stop()
            throw e
        }
    }

    private void onNext(DataflowWriteChannel target, Object value, int index) {
        if( value !instanceof RecordMap )
            throw new ScriptRuntimeException("Operator `join` expected a record but received: ${value} [${value.class.simpleName}]")
        final record = value as RecordMap
        if( predicate != null )
            onNextPredicate(target, record, index)
        else
            onNextKey(target, record, index)
    }

    private void onNextKey(DataflowWriteChannel target, RecordMap record, int index) {
        if( index == 0 ) {
            final key = keyOf(record)
            final leftValues = mappingsLeft.computeIfAbsent(key, (k) -> new ArrayList<>())
            final rightValues = mappingsRight.computeIfAbsent(key, (k) -> new ArrayList<>())
            for( final rightValue : rightValues )
                target << record.plus(rightValue)
            if( !rightValues.isEmpty() )
                emitted.add(key)
            leftValues.add(record)
        }
        else if( index == 1 ) {
            final key = keyOf(record)
            final leftValues = mappingsLeft.computeIfAbsent(key, (k) -> new ArrayList<>())
            final rightValues = mappingsRight.computeIfAbsent(key, (k) -> new ArrayList<>())
            for( final leftValue : leftValues )
                target << leftValue.plus(record)
            if( !leftValues.isEmpty() )
                emitted.add(key)
            rightValues.add(record)
        }
    }

    private Object keyOf(RecordMap record) {
        return pivot != null ? record[pivot] : keyFn.call(record)
    }

    // nested-loop join, O(n*m) predicate calls
    private void onNextPredicate(DataflowWriteChannel target, RecordMap record, int index) {
        if( index == 0 ) {
            final i = leftBuf.size()
            for( int j = 0; j < rightBuf.size(); j++ ) {
                if( predicate.call(record, rightBuf.get(j)) ) {
                    target << record.plus(rightBuf.get(j))
                    leftMatched.set(i)
                    rightMatched.set(j)
                }
            }
            leftBuf.add(record)
        }
        else if( index == 1 ) {
            final j = rightBuf.size()
            for( int i = 0; i < leftBuf.size(); i++ ) {
                if( predicate.call(leftBuf.get(i), record) ) {
                    target << leftBuf.get(i).plus(record)
                    leftMatched.set(i)
                    rightMatched.set(j)
                }
            }
            rightBuf.add(record)
        }
    }

    private synchronized void onComplete(DataflowWriteChannel target, int index) {
        if( failed )
            return

        if( --count > 0 )
            return

        if( remainder && predicate != null ) {
            // emit unmatched values from the predicate (closure) join
            for( int i = 0; i < leftBuf.size(); i++ )
                if( !leftMatched.get(i) )
                    target << leftBuf.get(i)
            for( int j = 0; j < rightBuf.size(); j++ )
                if( !rightMatched.get(j) )
                    target << rightBuf.get(j)
        }
        else if( remainder ) {
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
