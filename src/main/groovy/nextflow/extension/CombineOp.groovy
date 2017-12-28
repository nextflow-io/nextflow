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
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Channel
import nextflow.Nextflow

import static nextflow.extension.DataflowHelper.addToList
import static nextflow.extension.DataflowHelper.split

/**
 * Implements the {@link DataflowExtensions#spread(groovyx.gpars.dataflow.DataflowReadChannel, java.lang.Object)} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CombineOp {

    static final List<Integer> NONE = []

    private DataflowReadChannel leftChannel

    private DataflowReadChannel rightChannel

    private DataflowQueue target

    private Map<Object,List> leftValues = [:]

    private Map<Object,List> rightValues = [:]

    private static final int LEFT = 0

    private static final int RIGHT = 1

    private List<Integer> pivot = NONE

    CombineOp(DataflowReadChannel left, Object right) {

        leftChannel = left

        switch(right) {
            case DataflowQueue:
            case DataflowExpression:
                rightChannel = (DataflowReadChannel)right;
                break

            case Collection:
            case (Object[]):
                rightChannel = Nextflow.channel(right as List)
                break

            default:
                throw new IllegalArgumentException("Not a valid argument for 'combine' operator [${right?.class?.simpleName}]: ${right} -- Use a List or a channel instead. ")
        }

    }

    CombineOp setPivot( pivot ) {
        this.pivot = pivot instanceof List ? (List)pivot : [pivot]
        return this
    }


    List<DataflowReadChannel> getInputs() {
        def result = [leftChannel]
        if( rightChannel )
            result << rightChannel
        return result
    }

    private Map handler(int index, DataflowWriteChannel target, AtomicInteger stopCount) {

        def opts = [:]
        opts.onNext = {
            if( pivot ) {
                def pair = split(pivot, it)
                emit(target, index, pair.keys, pair.values)
            }
            else {
                emit(target, index, NONE, it)
            }
        }

        opts.onComplete = {
            if( stopCount.decrementAndGet()==0) {
                target << Channel.STOP
            }}

        return opts
    }

    @PackageScope
    @CompileDynamic
    def tuple( List p, a, b ) {
        List result = new LinkedList()
        result.addAll(p)
        addToList(result, a)
        addToList(result, b)

        result.size()==1 ? result[0] : result
    }

    @PackageScope
    synchronized void emit( DataflowWriteChannel target, int index, List p, v ) {

        if( leftValues[p] == null ) leftValues[p] = []
        if( rightValues[p] == null ) rightValues[p] = []

        if( index == LEFT ) {
            log.trace "combine >> left >> by=$p; val=$v; right-values: ${rightValues[p]}"
            for ( Object x : rightValues[p] ) {
                target.bind( tuple(p, v, x) )
            }
            leftValues[p].add(v)
            return
        }

        if( index == RIGHT ) {
            log.trace "combine >> right >> by=$p; val=$v; right-values: ${leftValues[p]}"
            for ( Object x : leftValues[p] ) {
                target.bind( tuple(p, x, v) )
            }
            rightValues[p].add(v)
            return
        }

        throw new IllegalArgumentException("Not a valid spread operator index: $index")
    }

    public DataflowReadChannel apply() {

        target = new DataflowQueue()

        if( rightChannel ) {
            final stopCount = new AtomicInteger(2)
            DataflowHelper.subscribeImpl( leftChannel, handler(LEFT, target, stopCount) )
            DataflowHelper.subscribeImpl( rightChannel, handler(RIGHT, target, stopCount) )
        }

        else if( rightValues != null ) {
            final stopCount = new AtomicInteger(1)
            DataflowHelper.subscribeImpl( leftChannel, handler(LEFT, target, stopCount) )
        }

        else
            throw new IllegalArgumentException("Not a valid spread operator state -- Missing right operand")

        return target
    }
}
