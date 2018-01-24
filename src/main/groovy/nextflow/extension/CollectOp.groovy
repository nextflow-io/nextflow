/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.extension

import static nextflow.util.CheckHelper.checkParams

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.util.ArrayBag
/**
 * Implements {@link DataflowExtensions#collect(groovyx.gpars.dataflow.DataflowReadChannel)}  operator
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

        def next = { append(result, it) }
        def done = { target << ( result ? new ArrayBag(normalise(result)) : Channel.STOP )  }

        DataflowHelper.subscribeImpl(source, [onNext:next, onComplete:done])
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
