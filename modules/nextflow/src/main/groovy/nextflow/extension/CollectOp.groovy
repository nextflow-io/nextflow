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
 */

package nextflow.extension

import static nextflow.util.CheckHelper.*

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.ContextGrouping
import nextflow.extension.op.Op
import nextflow.util.ArrayBag

/**
 * Implements {@link OperatorImpl#collect(groovyx.gpars.dataflow.DataflowReadChannel)}  operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CollectOp {

    private static final Map COLLECT_PARAMS = [ flat: Boolean, sort: [Boolean, Comparator, Closure] ]

    private DataflowReadChannel source

    private sort

    private flat

    private Closure action

    CollectOp( DataflowReadChannel source, Closure action, Map params=null ) {
        checkParams('collect', params, COLLECT_PARAMS)
        this.source = source
        this.action = action
        this.sort = params?.sort
        this.flat = params?.flat
    }

    DataflowVariable apply( ) {
        final result = []
        final target = new DataflowVariable()

        new SubscribeOp()
            .withInput(source)
            .withContext(new ContextGrouping())
            .withOnNext { append(result, it) }
            .withOnComplete { DataflowProcessor dp ->
                final msg = result ? new ArrayBag(normalise(result)) : Channel.STOP
                Op.bind(dp, target, msg)
            }
            .apply()

        return target
    }

    private void append(List list, item) {

        final obj = action ? action.call(item) : item

        if( obj instanceof List && (flat==null || flat==true) ) {
            list.addAll( (List)obj )
        }
        else
            list.add(obj)
    }

    private List normalise(List list) {
        if( !sort )
            return list

        if( sort == true )
            return list.sort()

        if( sort instanceof Closure )
            return list.sort( true, (Closure)sort )

        if( sort instanceof Comparator )
            return list.sort( true, (Comparator)sort )

        throw new IllegalArgumentException("Not a valid `collect` sort parameter: $sort [${sort.getClass().getName()}]")
    }

}
