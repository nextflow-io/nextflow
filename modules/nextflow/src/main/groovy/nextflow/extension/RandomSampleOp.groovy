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

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.ContextSequential
import nextflow.extension.op.Op
import nextflow.extension.op.OpContext
import nextflow.extension.op.OpDatum
/**
 * Implements Reservoir sampling of channel content
 *
 * See https://en.wikipedia.org/wiki/Reservoir_sampling
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class RandomSampleOp {

    private DataflowReadChannel source

    private DataflowWriteChannel result

    private int N

    private Random rng

    private List<OpDatum> reservoir = []

    private int counter

    private OpContext context = new ContextSequential()

    RandomSampleOp( DataflowReadChannel source, int N, Long seed = null) {
        this.source = source
        this.N = N
        this.rng = seed != null ? new Random(seed) : new Random()
    }

    private void sampling(Object it) {
        counter++
        //Fill reservoir
        if (counter <= N){
            reservoir.add(OpDatum.of(it, context.getOperatorRun()))
        }
        else {
            //Pick a random number
            int i = rng.nextInt(counter)
            if (i < N)
                reservoir[i] = OpDatum.of(it, context.getOperatorRun())
        }
    }

    private void emit(DataflowProcessor dp) {
        if( counter <= N )
            Collections.shuffle(reservoir, rng)
        for( OpDatum it : reservoir ) {
            if( it!=null )
            Op.bind(it.run, result, it.value)
        }
        Op.bind(dp, result, Channel.STOP)
    }

    DataflowWriteChannel apply() {
        result = CH.create()
        new SubscribeOp()
            .withInput(source)
            .withContext(context)
            .withOnNext(this.&sampling)
            .withOnComplete(this.&emit)
            .apply()
        return result
    }
}
