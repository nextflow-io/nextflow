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
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.types.Tuple
import nextflow.util.ArrayTuple
/**
 * Implements the `cross` operator for typed dataflow.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class CrossOpV2 {

    private DataflowReadChannel left

    private DataflowReadChannel right

    CrossOpV2(DataflowReadChannel left, DataflowReadChannel right) {
        this.left = left
        this.right = right
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

    private Collection<?> leftValues = []
    private Collection<?> rightValues = []
    private int count = 2

    private synchronized void onNext(DataflowWriteChannel target, Object value, int index) {
        if( index == 0 ) {
            final leftValue = value
            for( final rightValue : rightValues )
                target << cross0(leftValue, rightValue)
            leftValues.add(leftValue)
        }
        else if( index == 1 ) {
            final rightValue = value
            for( final leftValue : leftValues )
                target << cross0(leftValue, rightValue)
            rightValues.add(rightValue)
        }
    }

    private static Tuple cross0(Object left, Object right) {
        return new ArrayTuple<>([left, right])
    }

    private synchronized void onComplete(DataflowWriteChannel target, int index) {
        count--
        if( count > 0 )
            return

        target << CH.stop()
    }

}
