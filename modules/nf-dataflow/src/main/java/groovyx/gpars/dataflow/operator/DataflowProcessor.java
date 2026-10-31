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
// Modified by Seqera Labs for Nextflow (2026): added the backpressure gate (setGate, resume, isPaused, getPauseCount); run by the mailbox engine instead of a GPars actor; added isActive(); removed registerChannelListenersToAllInputs,
// start(PGroup), terminateAfterNextRun, bindAllOutputValuesAtomically, the output and state getters, the state object, adding and removing listeners after construction and the
// registered, afterStart, messageSentOut and customEvent events.

package groovyx.gpars.dataflow.operator;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.function.BooleanSupplier;

/**
 * Dataflow selectors and operators (processors) form the basic units in dataflow networks. They are typically combined into oriented graphs that transform data.
 * They accept a set of input and output dataflow channels and following specific strategies they transform input values from the input channels
 * into new values written to the output channels.
 * The output channels at the same time are suitable to be used as input channels by some other dataflow processors.
 * The channels allow processors to communicate.
 * <p>
 * Dataflow selectors and operators enable creation of highly concurrent applications yet the abstraction hides the low-level concurrency primitives
 * and exposes much friendlier API.
 * Operators run on a shared pool of threads, so the actual number of threads used by the calculation can be kept much lower
 * than the actual number of processors used in the network.
 * </p>
 *
 * @author Vaclav Pech
 *         Date: Sep 9, 2009
 */
@SuppressWarnings({"RawUseOfParameterizedType", "AccessingNonPublicFieldOfAnotherObject", "unchecked", "AbstractClassWithoutAbstractMethods", "ConstantDeclaredInAbstractClass"})
public abstract class DataflowProcessor {

    public static final String INPUTS = "inputs";
    public static final String OUTPUTS = "outputs";
    public static final String MAX_FORKS = "maxForks";
    private static final String LISTENERS = "listeners";

    /**
     * The engine running the processor
     */
    MailboxEngine engine;

    final List<DataflowReadChannel<?>> inputs;

    final List<DataflowWriteChannel<?>> outputs;

    protected final Collection<DataflowEventListener> listeners = new CopyOnWriteArrayList<DataflowEventListener>();

    private volatile BooleanSupplier gate;

    /**
     * Creates a processor
     * After creation the processor needs to be started using the start() method.
     *
     * @param channels A map specifying "inputs" and "outputs" - dataflow channels (instances of the DataflowQueue or DataflowVariable classes) to use for inputs and outputs
     * @param code     The processor's body to run each time all inputs have a value to read
     */
    protected DataflowProcessor(final Map<String, Object> channels, final Closure<?> code) {
        //noinspection ThisEscapedInObjectConstruction
        code.setDelegate(this);

        listeners.addAll(extractListeners(channels));
        inputs = channels != null ? extractInputs(channels) : Collections.<DataflowReadChannel<?>>emptyList();
        outputs = channels != null ? extractOutputs(channels) : Collections.<DataflowWriteChannel<?>>emptyList();
        if (channels == null) return;
        final Collection<?> declaredInputs = (Collection<?>) channels.get(INPUTS);
        if (declaredInputs == null || declaredInputs.isEmpty()) {
            throw new IllegalArgumentException("The processor body must take some inputs. The provided list of input channels is empty.");
        }
    }

    static boolean shouldBeMultiThreaded(final Map<String, Object> channels) {
        final Integer maxForks = (Integer) channels.get(MAX_FORKS);
        return maxForks != null && maxForks != 1;
    }

    static List<DataflowReadChannel<?>> extractInputs(final Map<String, Object> channels) {
        final List<DataflowReadChannel<?>> inputs = (List<DataflowReadChannel<?>>) channels.get(INPUTS);
        if (inputs == null) return Collections.emptyList();
        return Collections.unmodifiableList(inputs);
    }

    static List<DataflowWriteChannel<?>> extractOutputs(final Map<String, Object> channels) {
        final List<DataflowWriteChannel<?>> outputs = (List<DataflowWriteChannel<?>>) channels.get(OUTPUTS);
        if (outputs == null) return Collections.emptyList();
        return Collections.unmodifiableList(outputs);
    }

    private static Collection<DataflowEventListener> extractListeners(final Map<String, Object> channels) {
        if (channels == null) return null;
        final Collection<DataflowEventListener> listeners = (Collection<DataflowEventListener>) channels.get(LISTENERS);
        return listeners != null ? listeners : Collections.<DataflowEventListener>emptyList();
    }

    protected static void checkMaxForks(final Map<?, ?> channels) {
        if ((Integer) channels.get(MAX_FORKS) < 1)
            throw new IllegalArgumentException("The maxForks argument must be a positive value. " + channels.get(MAX_FORKS) + " was provided.");
    }

    /**
     * Starts the processor
     *
     * @return This operator instance
     */
    public final DataflowProcessor start() {
        engine.start();
        return this;
    }

    /**
     * Stops the processor immediately, potentially loosing unhandled messages
     */
    public final void terminate() {
        engine.terminate();
    }

    /**
     * Sets a gate that pauses the operator: while the gate returns true, the operator does not
     * read further inputs after a firing. Call {@link #resume()} once the gate may have opened.
     */
    public final void setGate(final BooleanSupplier gate) {
        this.gate = gate;
    }

    final boolean isGateClosed() {
        final BooleanSupplier g = gate;
        return g != null && g.getAsBoolean();
    }

    /**
     * Wakes the operator if it is paused by its gate
     */
    public final void resume() {
        engine.resume();
    }

    public final boolean isPaused() {
        return engine.isPaused();
    }

    /**
     * The number of times the operator was paused by its gate
     */
    public final long getPauseCount() {
        return engine.getPauseCount();
    }

    /**
     * Joins the processor waiting for it to finish
     *
     * @throws InterruptedException If the thread gets interrupted
     */
    public final void join() throws InterruptedException {
        engine.join();
    }

    /**
     * Whether the processor was started and has not stopped yet
     */
    public final boolean isActive() {
        return engine.isActive();
    }

    /**
     * Used by the processor's body to send a value to the given output channel
     *
     * @param idx   The index of the channel to bind
     * @param value The value to bind
     */
    public final void bindOutput(final int idx, final Object value) {
        ((DataflowWriteChannel<Object>) outputs.get(idx)).bind(value);
    }

    /**
     * Used by the processor's body to send a value to the first / only output channel
     *
     * @param value The value to bind
     */
    public final void bindOutput(final Object value) {
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
    public final void bindAllOutputs(final Object value) {
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
    public final void bindAllOutputValues(final Object... values) {
        for (int i = 0; i < outputs.size(); i++)
            ((DataflowWriteChannel<Object>) outputs.get(i)).bind(values[i]);
    }

    /**
     * Used by the processor's body to send a value to all output channels, while guaranteeing atomicity of the operation
     * and preventing other calls to bindAllOutputsAtomically() from interfering with one another.
     *
     * @param value The value to bind
     */
    public final synchronized void bindAllOutputsAtomically(final Object value) {
        bindAllOutputs(value);
    }

    /**
     * Is invoked in case the body throws an exception.
     *
     * @param e The reported exception
     */
    @SuppressWarnings({"UseOfSystemOutOrSystemErr"})
    final synchronized void reportError(final Throwable e) {
        if (listeners.isEmpty()) {
            System.err.println("The dataflow processor experienced an unhandled exception and is about to terminate. " + e);
            terminate();
        }
        if (fireOnException(e)) {
            terminate();
        }
    }

    protected final void fireAfterStop() {
        for (final DataflowEventListener listener : listeners) {
            listener.afterStop(this);
        }
    }

    protected final boolean fireOnException(final Throwable e) {
        boolean terminate = false;
        for (final DataflowEventListener listener : listeners) {
            terminate = terminate || listener.onException(this, e);
        }
        return terminate;
    }

    protected final Object fireMessageArrived(final DataflowReadChannel channel, final int index, final Object message) {
        Object result = message;
        for (final DataflowEventListener listener : listeners) {
            result = listener.messageArrived(this, channel, index, result);
        }
        return result;
    }

    protected final Object fireControlMessageArrived(final DataflowReadChannel channel, final int index, final Object message) {
        Object result = message;
        for (final DataflowEventListener listener : listeners) {
            result = listener.controlMessageArrived(this, channel, index, result);
        }
        return result;
    }

    protected final List<Object> fireBeforeRun(final List<Object> messages) {
        List<Object> result = messages;
        for (final DataflowEventListener listener : listeners) {
            result = listener.beforeRun(this, result);
        }
        return result;
    }

    protected final void fireAfterRun(final List<Object> messages) {
        for (final DataflowEventListener listener : listeners) {
            listener.afterRun(this, messages);
        }
    }
}
