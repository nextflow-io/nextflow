package nextflow.extension

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel

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

    private DataflowQueue result

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

    DataflowQueue apply() {
        result = new DataflowQueue()
        DataflowHelper.subscribeImpl(source, [onNext: this.&sampling, onComplete: this.&emit])
        return result
    }
}
