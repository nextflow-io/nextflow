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
// Modified by Seqera Labs for Nextflow (2026): merged AgentCore and AgentBase into this class; handles the messages
// on the dataflow pool instead of the GPars actor messaging core, and in place when the loop thread waits for them;
// removed fair agents, validators, the copy closure, the static factories, getInstantVal, valAsync, the public
// sendAndWait, setPGroup, attachToThreadPool and the leftShift and call aliases of send.

package groovyx.gpars.agent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.Executor;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.locks.ReentrantLock;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.scheduler.Scheduler;
import org.codehaus.groovy.runtime.InvokerInvocationException;
import org.codehaus.groovy.runtime.NullObject;

/**
 * A thread-safe reference to mutable state, inspired by Agents in Clojure.
 * <p>
 * The state is not accessed directly: closures are sent to the agent, which runs them one at a time, in send
 * order, on the dataflow pool, with the state as argument and the agent as delegate. A closure may replace the
 * state with {@link #updateValue}. A message that is not a closure replaces the state. An error thrown by a
 * closure is collected, see {@link #getErrors()}, and the agent goes on with the next message.
 *
 * @author Vaclav Pech
 *         Date: Jul 2, 2009
 */
@SuppressWarnings({"rawtypes", "unchecked"})
public class Agent<T> {

    private static final int BATCH = 128;

    private final Queue<Object> queue = new ConcurrentLinkedQueue<Object>();

    private final AtomicBoolean scheduled = new AtomicBoolean();

    // held while handling messages, so that they are handled one at a time and in order
    private final ReentrantLock drainLock = new ReentrantLock();

    private final Runnable drain = this::drain;

    private final Executor pool;

    private final Collection<Closure> listeners = new CopyOnWriteArrayList<Closure>();

    private volatile Thread currentThread;

    private T data;

    private List<Throwable> errors;

    /**
     * Creates a new Agent with the internal state set to null
     */
    public Agent() {
        this(null);
    }

    /**
     * Creates a new Agent around the supplied modifiable object
     *
     * @param data The object to use for storing the internal state of the variable
     */
    public Agent(final T data) {
        this(Scheduler.pool(), data);
    }

    /**
     * Creates an agent that runs on the given pool instead of the dataflow pool
     */
    Agent(final Executor pool, final T data) {
        this.pool = pool;
        this.data = data;
    }

    /**
     * Adds the message to the agent's message queue
     *
     * @param message A closure, or a new value for the state
     */
    public final void send(final Object message) {
        queue.add(message != null ? message : NullObject.getNullObject());
        schedule();
    }

    private void schedule() {
        if (!scheduled.compareAndSet(false, true))
            return;
        try {
            pool.execute(drain);
        }
        catch (RuntimeException | Error e) {
            scheduled.set(false);
            throw e;
        }
    }

    /**
     * Handles a batch of messages, then gives up the thread
     */
    private void drain() {
        drainLock.lock();
        try {
            handleMessages(BATCH, null);
        }
        finally {
            drainLock.unlock();
            // nor to the next pool task, e.g. an interrupt that came after the last message
            Thread.interrupted();
            scheduled.set(false);
            if (!queue.isEmpty())
                schedule();
        }
    }

    /**
     * Handles up to {@code max} messages, or until the given variable is bound
     */
    private void handleMessages(final int max, final DataflowVariable<?> until) {
        currentThread = Thread.currentThread();
        try {
            for (int n = 0; n < max && (until == null || !until.isBound()); n++) {
                final Object message = queue.poll();
                if (message == null)
                    break;
                try {
                    handleMessage(message);
                }
                catch (InvokerInvocationException e) {
                    registerError(e.getCause());
                }
                catch (Throwable e) {
                    registerError(e);
                }
                // like GPars, an interrupt set by a closure is an error, cleared so that it does not leak to the next message
                if (Thread.interrupted()) {
                    registerError(new InterruptedException());
                    break;
                }
            }
        }
        finally {
            currentThread = null;
        }
    }

    private void handleMessage(final Object message) {
        if (message instanceof Closure) {
            final Closure code = (Closure) message;
            code.setDelegate(this);
            code.call(data);
        }
        else if (message instanceof NullObject)
            updateValue(null);
        else
            updateValue((T) message);
    }

    /**
     * Allows closures to set the new internal state as a whole.
     * Do not call this method directly from the outside. It is designed to be used from within the submitted closures.
     *
     * @param newValue The value to set the internal state to
     */
    public final void updateValue(final T newValue) {
        final T oldValue = data;
        data = newValue;
        for (final Closure listener : listeners)
            listener.call(new Object[]{oldValue, newValue});
    }

    /**
     * Retrieves the state once all the messages sent before are handled
     *
     * @return The value of the agent's state
     * @throws InterruptedException If the thread gets interrupted while waiting
     */
    public final T getVal() throws InterruptedException {
        return (T) sendAndWait();
    }

    /**
     * Blocks until all the messages sent before are handled
     *
     * @throws InterruptedException If the thread gets interrupted while waiting
     */
    public final void await() throws InterruptedException {
        sendAndWait();
    }

    /**
     * Sends a message that returns the state, and waits for it. On the loop thread of a single-thread pool,
     * which the agent also runs on, the messages are handled in place instead.
     */
    private Object sendAndWait() throws InterruptedException {
        if (Thread.currentThread() == currentThread)
            throw new IllegalStateException("Cannot submit messages to agents inside submitted commands");
        // like a blocking wait, so that the inline path does not turn the interrupt into an agent error
        if (Thread.interrupted())
            throw new InterruptedException();
        final DataflowVariable<Object> result = new DataflowVariable<Object>();
        send(new Closure(this) {
            @Override
            public Object call(final Object state) {
                result.bind(state);
                return null;
            }
        });
        if (Scheduler.isLoopThread() && drainLock.tryLock()) {
            try {
                // an interrupt set by a closure stops a batch, not the wait
                do {
                    handleMessages(Integer.MAX_VALUE, result);
                } while (!result.isBound() && !queue.isEmpty());
            }
            finally {
                drainLock.unlock();
            }
        }
        return result.getVal();
    }

    /**
     * Adds a listener of the state updates: a closure accepting the old and the new value in this order,
     * plus optionally the agent as the first argument
     *
     * @param listener The closure with two or three arguments
     */
    public void addListener(final Closure listener) {
        final int parameters = listener.getMaximumNumberOfParameters();
        if (parameters < 2 || parameters > 3)
            throw new IllegalArgumentException("Agent listeners can only take two arguments plus optionally the current agent instance as the first argument.");
        listeners.add(parameters == 3 ? listener.curry(new Object[]{this}) : listener);
    }

    private synchronized void registerError(final Throwable e) {
        if (errors == null)
            errors = new ArrayList<Throwable>();
        errors.add(e);
    }

    /**
     * Retrieves the errors thrown by the closures, and clears them
     *
     * @return A detached collection of the errors
     */
    public synchronized List<Throwable> getErrors() {
        if (errors == null)
            return Collections.emptyList();
        try {
            return errors;
        }
        finally {
            errors = null;
        }
    }

    /**
     * Indicates whether the closures have thrown errors
     */
    public synchronized boolean hasErrors() {
        return errors != null && !errors.isEmpty();
    }
}
