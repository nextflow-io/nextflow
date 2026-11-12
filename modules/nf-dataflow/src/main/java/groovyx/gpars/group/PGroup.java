// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2012  The original author or authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Modified by Seqera Labs for Nextflow (2026): removed the actor, select and selector factory methods, tasks, splitters and whenAllBound; no pool shutdown on finalization, since groups share the scheduler pool.

package groovyx.gpars.group;

import groovy.lang.Closure;
import groovy.lang.DelegatesTo;
import groovyx.gpars.agent.Agent;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import groovyx.gpars.dataflow.operator.DataflowOperator;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import groovyx.gpars.scheduler.Pool;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static java.util.Arrays.asList;

/**
 * Provides a common super class of pooled parallel groups.
 *
 * @author Vaclav Pech, Alex Tkachman
 *         Date: May 8, 2009
 */
@SuppressWarnings({"RawUseOfParameterizedType"})
public abstract class PGroup {

    /**
     * Stores the group actors' thread pool
     */
    private final Pool threadPool;

    public Pool getThreadPool() {
        return threadPool;
    }

    /**
     * Creates a group for actors, agents, tasks and operators. The actors will share a common daemon thread pool.
     *
     * @param threadPool The thread pool to use by the group
     */
    protected PGroup(final Pool threadPool) {
        this.threadPool = threadPool;
    }

    /**
     * Creates an agent instance initialized with the given state
     *
     * @param state The initial internal state of the new Agent instance
     * @return The created instance
     */
    public final <T> Agent<T> agent(final T state) {
        final Agent<T> safe = new Agent<T>(state);
        safe.attachToThreadPool(threadPool);
        return safe;
    }

    /**
     * Creates an agent instance initialized with the given state
     *
     * @param state The initial internal state of the new Agent instance
     * @param copy  A closure to use to create a copy of the internal state when sending the internal state out
     * @return The created instance
     */
    public final <T> Agent<T> agent(final T state, final Closure copy) {
        final Agent<T> safe = new Agent<T>(state, copy);
        safe.attachToThreadPool(threadPool);
        return safe;
    }

    /**
     * Creates an agent instance initialized with the given state, which will cooperate in thread sharing with other agents and actors in a fair manner.
     *
     * @param state The initial internal state of the new Agent instance
     * @return The created instance
     */
    public final <T> Agent<T> fairAgent(final T state) {
        final Agent<T> safe = agent(state);
        safe.makeFair();
        return safe;
    }

    /**
     * Creates an agent instance initialized with the given state, which will cooperate in thread sharing with other agents and actors in a fair manner.
     *
     * @param copy  A closure to use to create a copy of the internal state when sending the internal state out
     * @param state The initial internal state of the new Agent instance
     * @return The created instance
     */
    public final <T> Agent<T> fairAgent(final T state, final Closure copy) {
        final Agent<T> safe = agent(state, copy);
        safe.makeFair();
        return safe;
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param channels A map specifying "inputs" and "outputs" - dataflow channels (instances of the DataflowQueue or DataflowVariable classes) to use for inputs and outputs
     * @param code     The operator's body to run each time all inputs have a value to read
     * @return A new started operator instance with all the channels set
     */
    public final DataflowProcessor operator(final Map channels, @DelegatesTo(DataflowOperator.class) final Closure code) {
        return new DataflowOperator(this, channels, code).start();
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param inputChannels  dataflow channels to use for input
     * @param outputChannels dataflow channels to use for output
     * @param code           The operator's body to run each time all inputs have a value to read
     * @return A new started operator instance with all the channels set
     */
    public final DataflowProcessor operator(final List inputChannels, final List outputChannels, @DelegatesTo(DataflowOperator.class) final Closure code) {
        final HashMap<String, List> params = new HashMap<String, List>(5);
        params.put(DataflowProcessor.INPUTS, inputChannels);
        params.put(DataflowProcessor.OUTPUTS, outputChannels);
        return new DataflowOperator(this, params, code).start();
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param inputChannels  dataflow channels to use for input
     * @param outputChannels dataflow channels to use for output
     * @param maxForks       Number of parallel threads running operator's body, defaults to 1
     * @param code           The operator's body to run each time all inputs have a value to read
     * @return A new started operator instance with all the channels set
     */
    public final DataflowProcessor operator(final List inputChannels, final List outputChannels, final int maxForks, @DelegatesTo(DataflowOperator.class) final Closure code) {
        final HashMap<String, Object> params = new HashMap<String, Object>(5);
        params.put(DataflowProcessor.INPUTS, inputChannels);
        params.put(DataflowProcessor.OUTPUTS, outputChannels);
        params.put(DataflowProcessor.MAX_FORKS, maxForks);
        return new DataflowOperator(this, params, code).start();
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param input  a dataflow channel to use for input
     * @param output a dataflow channel to use for output
     * @param code   The operator's body to run each time all inputs have a value to read
     * @return A new started operator instance with all the channels set
     */
    public final DataflowProcessor operator(final DataflowReadChannel input, final DataflowWriteChannel output, @DelegatesTo(DataflowOperator.class) final Closure code) {
        final HashMap<String, List> params = new HashMap<String, List>(5);
        params.put(DataflowProcessor.INPUTS, asList(input));
        params.put(DataflowProcessor.OUTPUTS, asList(output));
        return new DataflowOperator(this, params, code).start();
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param input    a dataflow channel to use for input
     * @param output   a dataflow channel to use for output
     * @param maxForks Number of parallel threads running operator's body, defaults to 1
     * @param code     The operator's body to run each time all inputs have a value to read
     * @return A new started operator instance with all the channels set
     */
    public final DataflowProcessor operator(final DataflowReadChannel input, final DataflowWriteChannel output, final int maxForks, @DelegatesTo(DataflowOperator.class) final Closure code) {
        final HashMap<String, Object> params = new HashMap<String, Object>(5);
        params.put(DataflowProcessor.INPUTS, asList(input));
        params.put(DataflowProcessor.OUTPUTS, asList(output));
        params.put(DataflowProcessor.MAX_FORKS, maxForks);
        return new DataflowOperator(this, params, code).start();
    }

    /**
     * Resizes the thread pool to the specified value
     *
     * @param poolSize The new pool size
     */
    public void resize(final int poolSize) {
        threadPool.resize(poolSize);
    }

    /**
     * Sets the pool size to the default
     */
    public void resetDefaultSize() {
        threadPool.resetDefaultSize();
    }

    /**
     * Retrieves the current thread pool size
     *
     * @return The pool size
     */
    public int getPoolSize() {
        return threadPool.getPoolSize();
    }

    /**
     * schedules a new task for processing with the pool
     *
     * @param task The task to schedule
     */
    public void execute(final Runnable task) {
        threadPool.execute(task);
    }

    /**
     * Gently stops the pool
     */
    public void shutdown() {
        threadPool.shutdown();
    }
}
