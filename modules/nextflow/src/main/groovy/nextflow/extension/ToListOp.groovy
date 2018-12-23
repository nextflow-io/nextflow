/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.extension

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.expression.DataflowExpression

/**
 * Implements {@link DataflowExtensions#toList(groovyx.gpars.dataflow.DataflowReadChannel)}  operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ToListOp {

    private DataflowReadChannel source

    private ordering

    /**
     * Implements both `toList` and `toSortedList`
     *
     * @param source The source channel
     * @param order Order the resulting channel, either boolean {@code true} or a comparing {@link Closure}
     */
    ToListOp( DataflowReadChannel source, order=null ) {
        this.source = source
        this.ordering = order
    }

    DataflowVariable apply() {
        assert source != null

        final target = new DataflowVariable()
        if( source instanceof DataflowExpression ) {
            final result = new ArrayList(1)
            Map<String,Closure> events = [:]
            events.onNext = { result.add(it) }
            events.onComplete = { target.bind(result) }
            DataflowHelper.subscribeImpl(source, events )
            return target
        }

        DataflowHelper.reduceImpl(source, target, []) { List list, item -> list << item }
        if( ordering ) {
            final sort = { List list -> ordering instanceof Closure ? list.sort((Closure) ordering) : list.sort() }
            return (DataflowVariable)target.then(sort)
        }
        else {
            return target
        }
    }

    static DataflowVariable apply( DataflowReadChannel source ) {
        new ToListOp(source).apply()
    }

}
