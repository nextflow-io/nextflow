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

package nextflow.script

import com.google.common.annotations.Beta
import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.NF
import nextflow.extension.CH
import nextflow.extension.MapOp
import nextflow.extension.MergeOp
import nextflow.extension.MixOp
import nextflow.extension.TakeOp
import nextflow.extension.UntilManyOp
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Implements common logic for recurse execution of nextflow processes and workflows
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Beta
@CompileStatic
trait IterableDef {

    abstract Object invoke_a(Object[] args)

    private List args

    private List inputChannels
    private List outputChannels

    private List<List<Object>> accumulators

    /**
     * To support component recursion, the feedback channels must be used as the channel instances
     * for the outputs instead of creating new ones
     * 
     * @return A list of output channels bringing the process feedback
     */
    List<DataflowWriteChannel> getFeedbackChannels() { outputChannels }

    /**
     * Implements a `recurse` idiom for a process or workflow
     *
     * @param args The process or workflow initial input
     * @return The object itself
     */
    IterableDef recurse(Object[] args) {
        if( !NF.recurseEnabled ) throw new MissingMethodException('recurse', this.getClass(), InvokerHelper.EMPTY_ARGS)
        this.args = ChannelOut.spread(args)
        checkRecurseArgs(this.args)
        return this
    }

    Object times(int n) {
        if( !NF.recurseEnabled ) throw new MissingMethodException('times', this.getClass(), InvokerHelper.EMPTY_ARGS)
        // expand the args
        createRecurseFeedbackWithTimes(n, args)
        return invoke_a(inputChannels.toArray())
    }

    Object until(Closure<Boolean> condition) {
        if( !NF.recurseEnabled ) throw new MissingMethodException('until', this.getClass(), InvokerHelper.EMPTY_ARGS)
        createRecurseFeedbackWithUntil(condition, args)
        return invoke_a(inputChannels.toArray())
    }

    private void createRecurseFeedbackWithUntil(Closure<Boolean> condition, List args) {
        this.inputChannels = new ArrayList(args.size())
        this.outputChannels = new ArrayList(args.size())
        final sources = new ArrayList(args.size())

        // create the output channel
        for( int i=0; i<args.size(); i++ ) {
            final it = args[i]
            // create the output channel
            final output = CH.create()
            outputChannels.add(output)
            // normalise input i-th to a dataflow channel
            sources << (it instanceof DataflowWriteChannel ? it : CH.value(it))
        }

        // create the feedback as list of many channel
        final readers = outputChannels.collect( it->CH.getReadChannel(it) )
        final feedbacks = new UntilManyOp(readers, condition).apply()

        // another iteration to create the new input mixing the original source + the feedback
        for( int i=0; i<sources.size(); i++ ) {
            final input = new MixOp( CH.getReadChannel(sources[i]), CH.getReadChannel(feedbacks[i]) ).apply()
            inputChannels << input
        }
    }

    static private void checkRecurseArg0(Object input, int index) {
        if( input !instanceof DataflowWriteChannel ) {
            // any non channel value is fine, because it will be wrapped into
            // a channel value
            return 
        }
        if( !CH.isValue(input) )
            throw new IllegalArgumentException("Recurse operation only allows value inputs -- Check ${index+1}-th argument")
        if( input==null )
            throw new IllegalArgumentException("Recurse operation does not allow null input values -- Check ${index+1}-th argument")
    }

    private List checkRecurseArgs(List items) {
        for( int i=0; i<items.size(); i++ ) {
            checkRecurseArg0(items.get(i), i)
        }
        return items
    }

    private void createRecurseFeedbackWithTimes(int times, List args) {
        this.inputChannels = new ArrayList(args.size())
        this.outputChannels = new ArrayList(args.size())

        checkRecurseArgs(args)

        for( int i=0; i<args.size(); i++ ) {
            // `it` represents the input value for the channel i-th
            final it = args[i]
            // create the output channel
            final output = CH.create()
            outputChannels.add(output)
            // normalise input i-th to a dataflow channel
            final source = it instanceof DataflowWriteChannel ? it : CH.value(it)
            // termination condition
            final feedback = new TakeOp(CH.getReadChannel(output), times-1).apply()
            // the input is made as the original source + the emission from the feedback
            final input = new MixOp( CH.getReadChannel(source), CH.getReadChannel(feedback) ).apply()
            inputChannels.add(input)
        }
    }

    Object scan(Object[] args) {
        if( !NF.recurseEnabled ) throw new MissingMethodException('scan', this.getClass(), InvokerHelper.EMPTY_ARGS)
        this.args = ChannelOut.spread(args)
        createScanFeedback(this.args)
        return invoke_a(inputChannels.toArray())
    }

    private void createScanFeedback(List args) {
        this.inputChannels = new ArrayList(args.size())
        this.outputChannels = new ArrayList(args.size())
        this.accumulators = new ArrayList<List<Object>>(args.size())

        for( int i=0; i<args.size(); i++ ) {
            // `it` represents the input value for the channel i-th
            final it = args[i]
            // create the output channel
            final output = CH.create()
            outputChannels.add(output)
            // normalise input i-th to a dataflow channel
            final source = it instanceof DataflowWriteChannel ? it : CH.value(it)
            // create the i-th channel accumulator
            accumulators[i] = new ArrayList<>()
            // create the feedback channel
            // this takes care to inject previous values back in the input
            final feedback = accumulator( CH.getReadChannel(output), i )
            // create the input composition
            final target = new MixOp(CH.getReadChannel(empty()), CH.getReadChannel(feedback)).apply()
            final input = new MergeOp(CH.getReadChannel(source), CH.getReadChannel(target)).apply()
            inputChannels.add(input)
        }
    }

    private DataflowWriteChannel accumulator(DataflowReadChannel source, int i) {
        final mapper = { value -> accumulators[i].add(value); return accumulators[i].clone() }
        return new MapOp(source, mapper).apply()
    }

    private DataflowWriteChannel empty() {
        CH.emitValues(CH.create(), [new ArrayList(0), Channel.STOP])
    }
}
