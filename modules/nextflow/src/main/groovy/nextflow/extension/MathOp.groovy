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


import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.ContextGrouping
import nextflow.extension.op.Op
/**
 * Implements the logic for "sum" and "mean" operators
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MathOp {

    private DataflowReadChannel source
    private DataflowVariable target
    private Closure action
    private Aggregate aggregate

    private MathOp(Aggregate aggregate) {
        this.aggregate = aggregate
    }

    static MathOp sum() {
        new MathOp(new Aggregate(name: 'sum'))
    }

    static MathOp mean() {
        new MathOp(new Aggregate(name: 'mean', mean: true))
    }

    MathOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    MathOp withTarget(DataflowVariable target) {
        assert target!=null
        this.target = target
        return this
    }

    MathOp withAction(Closure action) {
        this.action = action
        this.aggregate.withAction(action)
        return this
    }

    DataflowWriteChannel apply() {
        assert source!=null
        if( target==null )
            target = new DataflowVariable()
        new SubscribeOp()
            .withInput(source)
            .withContext(new ContextGrouping())
            .withOnNext(aggregate.&process)
            .withOnComplete(this.&completion)
            .apply()
        return target
    }

    private void completion(DataflowProcessor dp) {
        Op.bind(dp, target, aggregate.result )
    }

    /**
     * Implements the logic for sum and mean operators
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @CompileDynamic
    static class Aggregate {
        def accum
        long count = 0
        boolean mean
        Closure action
        String name

        def process(it) {
            if( it == null || it == Channel.VOID )
                return

            count++

            def item = action != null ? action.call(it) : it
            if( accum == null )
                accum = item

            else if( accum instanceof Number )
                accum += item

            else if( accum instanceof List && item instanceof List)
                for( int i=0; i<accum.size() && i<item.size(); i++ )
                    accum[i] += item.get(i)

            else
                throw new IllegalArgumentException("Invalid `$name` item: $item [${item.class.simpleName}]")
        }

        def getResult() {
            if( !mean || count == 0 )
                return accum

            if( accum instanceof List )
                return accum.collect { it / count }
            else
                return accum / count
        }

        Aggregate withAction(Closure action) {
            this.action = action
            return this
        }
    }

}
