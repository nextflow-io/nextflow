/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import nextflow.Channel

import static DataflowHelper.eventsMap
import static DataflowHelper.subscribeImpl

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

    private List reservoir = []

    private int counter

    RandomSampleOp( DataflowReadChannel source, int N, Long seed = null) {
        this.source = source
        this.N = N
        this.rng = seed != null ? new Random(seed) : new Random()
    }


    private void sampling(it) {

        counter++
        //Fill reservoir
        if (counter <= N){
            reservoir << it
        }
        else {
            //Pick a random number
            int i = rng.nextInt(counter)
            if (i < N)
                reservoir[i] = it
        }

    }

    private void emit(nop) {

        if( counter <= N )
            Collections.shuffle(reservoir, rng)

        reservoir.each { it!=null ? result.bind(it) : null }
        result.bind(Channel.STOP)
    }

    DataflowWriteChannel apply() {
        result = CH.create()
        subscribeImpl(source, eventsMap(this.&sampling, this.&emit))
        return result
    }
}
