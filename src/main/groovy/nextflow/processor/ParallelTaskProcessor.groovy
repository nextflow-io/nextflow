/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow.processor
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.script.EachInParam
/**
 * Defines the parallel tasks execution logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@InheritConstructors
class ParallelTaskProcessor extends TaskProcessor {

    /**
     * Keeps track of the task instance executed by the current thread
     */
    protected final ThreadLocal<TaskRun> currentTask = new ThreadLocal<>()

    @Override
    protected void createOperator() {

        def opInputs = new ArrayList(config.getInputs().getChannels())

        /*
         * check if there are some iterators declaration
         * the list holds the index in the list of all *inputs* for the {@code each} declaration
         */
        def iteratorIndexes = []
        config.getInputs().eachWithIndex { param, index ->
            if( param instanceof EachInParam ) {
                log.trace "Process ${name} > got each param: ${param.name} at index: ${index} -- ${param.dump()}"
                iteratorIndexes << index
            }
        }

        /*
         * When one (or more) {@code each} are declared as input, it is created an extra
         * operator which will receive the inputs from the channel (excepts the values over iterate)
         *
         * The operator will *expand* the received inputs, iterating over the user provided value and
         * forwarding the final values the the second *parallel* processor executing the user specified task
         */
        if( iteratorIndexes ) {
            log.debug "Creating *combiner* operator for each param(s) at index(es): ${iteratorIndexes}"

            // don't care about the last channel, being the control channel it doesn't bring real values
            final size = opInputs.size()-1

            // the iterator operator needs to executed just one time
            // thus add a dataflow queue binding a single value and then a stop signal
            def termination = new DataflowQueue<>()
            termination << Boolean.TRUE
            opInputs[size] = termination

            // the channel forwarding the data from the *iterator* process to the target task
            final linkingChannels = new ArrayList(size)
            size.times { linkingChannels[it] = new DataflowQueue() }

            // the script implementing the iterating process
            final forwarder = createForwardWrapper(size, iteratorIndexes)

            // instantiate the iteration process
            def stopAfterFirstRun = allScalarValues
            def interceptor = new BaseProcessInterceptor(opInputs, stopAfterFirstRun)
            def params = [inputs: opInputs, outputs: linkingChannels, maxForks: 1, listeners: [interceptor]]
            session.allProcessors << (processor = new DataflowOperator(group, params, forwarder))
            // fix issue #41
            processor.start()

            // set as next inputs the result channels of the iteration process
            // adding the 'control' channel removed previously
            opInputs = new ArrayList(size+1)
            opInputs.addAll( linkingChannels )
            opInputs.add( config.getInputs().getChannels().last() )
        }


        /*
         * create a mock closure to trigger the operator
         */
        final wrapper = createCallbackWrapper( opInputs.size(), this.&invokeTask )

        /*
         * define the max forks attribute:
         * - by default the process execution is parallel using the poolSize value
         * - otherwise use the value defined by the user via 'taskConfig'
         */
        def maxForks = session.poolSize
        if( config.maxForks ) {
            maxForks = config.maxForks
            blocking = true
        }
        log.debug "Creating operator > $name -- maxForks: $maxForks"

        /*
         * finally create the operator
         */
        // note: do not specify the output channels in the operator declaration
        // this allows us to manage them independently from the operator life-cycle
        def stopAfterFirstRun = allScalarValues && !hasEachParams
        def interceptor = new TaskProcessorInterceptor(opInputs, stopAfterFirstRun)
        def params = [inputs: opInputs, maxForks: maxForks, listeners: [interceptor] ]
        session.allProcessors << (processor = new DataflowOperator(group, params, wrapper))

        // notify the creation of a new vertex the execution DAG
        session.dag.addProcessNode(name, config.getInputs(), config.getOutputs())

        // fix issue #41
        processor.start()

    }

    /**
     * Implements the closure which *combines* all the iteration
     *
     * @param numOfInputs Number of in/out channel
     * @param indexes The list of indexes which identify the position of iterators in the input channels
     * @return The closure implementing the iteration/forwarding logic
     */
    protected createForwardWrapper( int len, List indexes ) {

        final args = []
        (len+1).times { args << "x$it" }

        final outs = []
        len.times { outs << "x$it" }

        /*
         * Explaining the following closure:
         *
         * - it has to be evaluated as a string since the number of input must much the number input channel
         *   that is known only at runtime
         *
         * - 'out' holds the list of all input values which need to be forwarded (bound) to output as many
         *   times are the items in the iteration list
         *
         * - the iteration list(s) is (are) passed like in the closure inputs like the other values,
         *   the *indexes* argument defines the which of them are the iteration lists
         *
         * - 'itr' holds the list of all iteration lists
         *
         * - using the groovy method a combination of all values is create (cartesian product)
         *   see http://groovy.codehaus.org/groovy-jdk/java/util/Collection.html#combinations()
         *
         * - the resulting values are replaced in the 'out' array of values and forwarded out
         *
         */

        final str =
            """
            { ${args.join(',')} ->
                def out = [ ${outs.join(',')} ]
                def itr = [ ${indexes.collect { 'x'+it }.join(',')} ]
                def cmb = itr.combinations()
                for( entries in cmb ) {
                    def count = 0
                    ${len}.times { i->
                        if( i in indexes ) { out[i] = entries[count++] }
                    }
                    bindAllOutputValues( *out )
                }
            }
            """


        final Binding binding = new Binding( indexes: indexes )
        final result = (Closure)grengine.run(str, binding)
        return result

    }


    /**
     * The processor execution body
     *
     * @param processor
     * @param values
     */
    final protected void invokeTask( def args ) {

        // create and initialize the task instance to be executed
        final List values = args instanceof List ? args : [args]
        if( log.isTraceEnabled() )
            log.trace "Setup new process > $name"

        // -- create the task run instance
        final task = createTaskRun()
        // -- set the task instance as the current in this thread
        currentTask.set(task)

        // -- map the inputs to a map and use to delegate closure values interpolation
        final secondPass = [:]
        int count = makeTaskContextStage1(task, secondPass, values)
        makeTaskContextStage2(task, secondPass, count)

        // verify that `when` guard, when specified, is satisfied
        if( !checkWhenGuard(task) )
            return

        // -- verify if exists a stored result for this case,
        //    if true skip the execution and return the stored data
        if( checkStoredOutput(task) )
            return

        def hash = createTaskHashKey(task)
        checkCachedOrLaunchTask(task, hash, resumable)
    }


    /**
     *  Intercept dataflow process events
     */
    class TaskProcessorInterceptor extends BaseProcessInterceptor {

        TaskProcessorInterceptor(List<DataflowChannel> inputs, boolean stop) {
            super(inputs, stop)
        }

        @Override
        public List<Object> beforeRun(final DataflowProcessor processor, final List<Object> messages) {
            log.trace "<${name}> Before run -- messages: ${messages}"
            // the counter must be incremented here, otherwise it won't be consistent
            state.update { StateObj it -> it.incSubmitted() }
            return messages;
        }


        @Override
        void afterRun(DataflowProcessor processor, List<Object> messages) {
            log.trace "<${name}> After run"
            currentTask.remove()
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def channelName = config.getInputs()?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Message arrived -- ${channelName} => ${message}"
            }

            super.messageArrived(processor, channel, index, message)
        }

        @Override
        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def channelName = config.getInputs()?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Control message arrived ${channelName} => ${message}"
            }

            super.controlMessageArrived(processor, channel, index, message)

            if( message == PoisonPill.instance ) {
                log.debug "<${name}> Poison pill arrived"
                state.update { StateObj it -> it.poison() }
            }

            return message
        }

        @Override
        public void afterStop(final DataflowProcessor processor) {
            log.debug "<${name}> After stop"
        }

        /**
         * Invoked if an exception occurs. Unless overridden by subclasses this implementation returns true to terminate the operator.
         * If any of the listeners returns true, the operator will terminate.
         * Exceptions outside of the operator's body or listeners' messageSentOut() handlers will terminate the operator irrespective of the listeners' votes.
         * When using maxForks, the method may be invoked from threads running the forks.
         * @param processor
         * @param error
         * @return
         */
        public boolean onException(final DataflowProcessor processor, final Throwable error) {
            // return `true` to terminate the dataflow processor
            handleException( error, currentTask.get() )
        }

    }

    /*
     * logger class for the *iterator* processor
     */
    class BaseProcessInterceptor extends DataflowEventAdapter {

        final List<DataflowChannel> inputs

        final boolean stopAfterFirstRun

        final int len

        final DataflowQueue control

        final int first

        BaseProcessInterceptor( List<DataflowChannel> inputs, boolean stop ) {
            this.inputs = new ArrayList<>(inputs)
            this.stopAfterFirstRun = stop
            this.len = inputs.size()
            this.control = (DataflowQueue)inputs.get(len-1)
            this.first = inputs.findIndexOf { it instanceof DataflowQueue }
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( len == 1 || stopAfterFirstRun ) {
                // -- kill itself
                control.bind(PoisonPill.instance)
            }
            else if( index == first ) {
                // -- keep things rolling
                control.bind(Boolean.TRUE)
            }

            return message;
        }

    }


}
