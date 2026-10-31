// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2013  The original author or authors
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
// Modified by Seqera Labs for Nextflow (2026): merged DataflowOperator into this class, which takes a typed constructor instead of a parameter map;
// run by the mailbox engine instead of a GPars actor; added the backpressure gate (setGate, resume, isPaused, getPauseCount) and isActive();
// removed registerChannelListenersToAllInputs, start(PGroup), terminateAfterNextRun, bindAllOutputValuesAtomically, the output and state getters,
// the state object, adding and removing listeners after construction and the registered, afterStart, messageSentOut and customEvent events.

package groovyx.gpars.dataflow.operator;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import groovyx.gpars.scheduler.Scheduler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Executor;
import java.util.function.BooleanSupplier;

/**
 * Dataflow operators (processors) form the basic units in dataflow networks. They are typically combined into oriented graphs that transform data.
 * Once a value is available in each of the input channels, the operator's body is triggered on the values, potentially generating values
 * to be written into the output channels.
 * The output channels at the same time are suitable to be used as input channels by some other dataflow processors.
 * The channels allow processors to communicate.
 * <p>
 * Operators run on a shared pool of threads, so the actual number of threads used by the calculation can be kept much lower
 * than the actual number of processors used in the network.
 * </p>
 *
 * @author Vaclav Pech
 *         Date: Sep 9, 2009
 */
@SuppressWarnings({"RawUseOfParameterizedType", "unchecked"})
public final class DataflowProcessor {

    /**
     * The engine running the processor
     */
    private final MailboxEngine engine;

    private final List<DataflowWriteChannel<?>> outputs;

    private final List<DataflowEventListener> listeners;

    private volatile BooleanSupplier gate;

    /**
     * Creates a processor on the dataflow pool.
     * After creation the processor needs to be started using the start() method.
     *
     * @param inputs    The input channels, one per parameter of the body
     * @param outputs   The output channels, or null for none
     * @param listeners The event listeners, or null for none
     * @param maxForks  The maximum number of bodies running at once
     * @param code      The processor's body to run each time all inputs have a value to read
     */
    public DataflowProcessor(final List<? extends DataflowReadChannel<?>> inputs, final List<? extends DataflowWriteChannel<?>> outputs, final List<? extends DataflowEventListener> listeners, final int maxForks, final Closure<?> code) {
        this(Scheduler.pool(), inputs, outputs, listeners, maxForks, code);
    }

    /**
     * Creates a processor that runs on the given pool instead of the dataflow pool
     */
    DataflowProcessor(final Executor pool, final List<? extends DataflowReadChannel<?>> inputs, final List<? extends DataflowWriteChannel<?>> outputs, final List<? extends DataflowEventListener> listeners, final int maxForks, final Closure<?> code) {
        if (inputs == null || inputs.isEmpty())
            throw new IllegalArgumentException("The processor body must take some inputs. The provided list of input channels is empty.");
        final int parameters = code.getMaximumNumberOfParameters();
        if (parameters != inputs.size())
            throw new IllegalArgumentException("The operator's body accepts " + parameters + " parameters while it is given " + inputs.size() + " input streams. The numbers must match.");
        if (maxForks < 1)
            throw new IllegalArgumentException("The maxForks argument must be a positive value. " + maxForks + " was provided.");
        this.outputs = copy(outputs);
        this.listeners = copy(listeners);
        //noinspection ThisEscapedInObjectConstruction
        code.setDelegate(this);
        this.engine = new MailboxEngine(this, pool, copy(inputs), (Closure<?>) code.clone(), maxForks);
    }

    private static <T> List<T> copy(final List<? extends T> list) {
        return list != null ? Collections.unmodifiableList(new ArrayList<T>(list)) : Collections.<T>emptyList();
    }

    /**
     * Starts the processor
     *
     * @return This operator instance
     */
    public DataflowProcessor start() {
        engine.start();
        return this;
    }

    /**
     * Stops the processor immediately, potentially loosing unhandled messages
     */
    public void terminate() {
        engine.terminate();
    }

    /**
     * Sets a gate that pauses the operator: while the gate returns true, the operator does not
     * read further inputs after a firing. Call {@link #resume()} once the gate may have opened.
     */
    public void setGate(final BooleanSupplier gate) {
        this.gate = gate;
    }

    boolean isGateClosed() {
        final BooleanSupplier g = gate;
        return g != null && g.getAsBoolean();
    }

    /**
     * Wakes the operator if it is paused by its gate
     */
    public void resume() {
        engine.resume();
    }

    public boolean isPaused() {
        return engine.isPaused();
    }

    /**
     * The number of times the operator was paused by its gate
     */
    public long getPauseCount() {
        return engine.getPauseCount();
    }

    /**
     * Joins the processor waiting for it to finish
     *
     * @throws InterruptedException If the thread gets interrupted
     */
    public void join() throws InterruptedException {
        engine.join();
    }

    /**
     * Whether the processor was started and has not stopped yet
     */
    public boolean isActive() {
        return engine.isActive();
    }

    /**
     * Used by the processor's body to send a value to the given output channel
     *
     * @param idx   The index of the channel to bind
     * @param value The value to bind
     */
    public void bindOutput(final int idx, final Object value) {
        ((DataflowWriteChannel<Object>) outputs.get(idx)).bind(value);
    }

    /**
     * Used by the processor's body to send a value to the first / only output channel
     *
     * @param value The value to bind
     */
    public void bindOutput(final Object value) {
        bindOutput(0, value);
    }

    /**
     * Used by the processor's body to send a value to all output channels.
     * If the maxForks value is set to a value greater than 1, calls to bindAllOutputs may result in values written to different
     * channels to be in different order. If this is a problem for the application logic, the bindAllOutputsAtomically
     * method should be considered instead.
     *
     * @param value The value to bind
     */
    public void bindAllOutputs(final Object value) {
        for (final DataflowWriteChannel<?> channel : outputs)
            ((DataflowWriteChannel<Object>) channel).bind(value);
    }

    /**
     * Used by the processor's body to send a value to all output channels. The values passed as arguments will each be sent
     * to an output channel with identical position index.
     * <p>
     * If the maxForks value is set to a value greater than 1, calls to bindAllOutputs may result in values written to different
     * channels to be in different order. If this is a problem for the application logic, the bindAllOutputsAtomically
     * method should be considered instead.
     * </p>
     *
     * @param values Values to send to output channels of the same position index
     */
    public void bindAllOutputValues(final Object... values) {
        for (int i = 0; i < outputs.size(); i++)
            ((DataflowWriteChannel<Object>) outputs.get(i)).bind(values[i]);
    }

    /**
     * Used by the processor's body to send a value to all output channels, while guaranteeing atomicity of the operation
     * and preventing other calls to bindAllOutputsAtomically() from interfering with one another.
     *
     * @param value The value to bind
     */
    public synchronized void bindAllOutputsAtomically(final Object value) {
        bindAllOutputs(value);
    }

    /**
     * Is invoked in case the body throws an exception.
     *
     * @param e The reported exception
     */
    @SuppressWarnings({"UseOfSystemOutOrSystemErr"})
    synchronized void reportError(final Throwable e) {
        if (listeners.isEmpty()) {
            System.err.println("The dataflow processor experienced an unhandled exception and is about to terminate. " + e);
            terminate();
        }
        if (fireOnException(e)) {
            terminate();
        }
    }

    void fireAfterStop() {
        for (final DataflowEventListener listener : listeners) {
            listener.afterStop(this);
        }
    }

    boolean fireOnException(final Throwable e) {
        boolean terminate = false;
        for (final DataflowEventListener listener : listeners) {
            terminate = terminate || listener.onException(this, e);
        }
        return terminate;
    }

    Object fireMessageArrived(final DataflowReadChannel channel, final int index, final Object message) {
        Object result = message;
        for (final DataflowEventListener listener : listeners) {
            result = listener.messageArrived(this, channel, index, result);
        }
        return result;
    }

    Object fireControlMessageArrived(final DataflowReadChannel channel, final int index, final Object message) {
        Object result = message;
        for (final DataflowEventListener listener : listeners) {
            result = listener.controlMessageArrived(this, channel, index, result);
        }
        return result;
    }

    List<Object> fireBeforeRun(final List<Object> messages) {
        List<Object> result = messages;
        for (final DataflowEventListener listener : listeners) {
            result = listener.beforeRun(this, result);
        }
        return result;
    }

    void fireAfterRun(final List<Object> messages) {
        for (final DataflowEventListener listener : listeners) {
            listener.afterRun(this, messages);
        }
    }
}
