/*
 * Copyright 2020-2021, Seqera Labs
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

package nextflow.script


import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.extension.CH
import nextflow.extension.MixOp
import nextflow.extension.TakeOp
/**
 * Implements common logic for recurse execution of nextflow processes and workflows
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait IterableDef {

    abstract Object invoke_a(Object[] args)

    private Object[] args

    private int times
    private int[] counters
    private List inputChannels
    private List outputChannels

    List getFeedbackChannels() { outputChannels }

    IterableDef loop(Object[] args) {
        this.args = args
        return this
    }

    Object times(int n) {
        // expand the args
        this.times = n
        createFeedbackChannels(ChannelOut.spread(args))
        return invoke_a(inputChannels.toArray())
    }

    private void createFeedbackChannels(List args) {
        this.inputChannels = new ArrayList(args.size())
        this.outputChannels = new ArrayList(args.size())
        this.counters = new int[ args.size() ]

        for( int i=0; i<args.size(); i++ ) {
            final it = args[i]
            // reset counters
            counters[i] = 0
            // normalise input i-th to a dataflow channel
            final source = it instanceof DataflowWriteChannel ? it : CH.value(it)
            // create the output channel
            final output = CH.create()
            // termination condition
            final feedback = new TakeOp(CH.getReadChannel(output), times-1).apply()
            // the input is made as the original source + the emission from the feedback
            final input = new MixOp( CH.getReadChannel(source), CH.getReadChannel(feedback) ).apply()
            inputChannels.add(input)
            outputChannels.add(output)
        }
    }

}
