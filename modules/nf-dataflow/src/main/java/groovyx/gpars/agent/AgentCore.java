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
// Modified by Seqera Labs for Nextflow (2026): a serial executor on the dataflow scheduler pool, instead of the GPars actor messaging core;
// removed setPGroup and the leftShift and call aliases of send.

package groovyx.gpars.agent;

import groovyx.gpars.scheduler.Pool;
import groovyx.gpars.scheduler.Scheduler;
import org.codehaus.groovy.runtime.InvokerInvocationException;
import org.codehaus.groovy.runtime.NullObject;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Handles the messages sent to an agent one at a time, in send order, on a pool thread
 *
 * @author Vaclav Pech
 *         Date: 13.4.2010
 */
@SuppressWarnings({"UnqualifiedStaticUsage"})
public abstract class AgentCore {

    private static final int BATCH = 128;

    private final Queue<Object> queue = new ConcurrentLinkedQueue<Object>();

    private final AtomicBoolean scheduled = new AtomicBoolean();

    private final Runnable drain = this::drain;

    private volatile Pool threadPool = Scheduler.pool();

    private volatile boolean fair = false;

    protected volatile Thread currentThread;

    /**
     * Sets a new thread pool to be used by the agent
     *
     * @param threadPool The thread pool to use
     */
    public final void attachToThreadPool(final Pool threadPool) {
        this.threadPool = threadPool;
    }

    /**
     * Retrieves the agent's fairness flag
     * Fair agents give up the thread after processing each message, non-fair agents keep a thread until their message queue is empty.
     * Non-fair agents tend to perform better than fair ones.
     *
     * @return True for fair agents, false for non-fair ones. Agents are non-fair by default.
     */
    public boolean isFair() {
        return fair;
    }

    /**
     * Makes the agent fair. Agents are non-fair by default.
     * Fair agents give up the thread after processing each message, non-fair agents keep a thread until their message queue is empty.
     * Non-fair agents tend to perform better than fair ones.
     */
    public void makeFair() {
        this.fair = true;
    }

    /**
     * Holds agent errors
     */
    private List<Throwable> errors;

    /**
     * Adds the message to the agent\s message queue
     *
     * @param message A value or a closure
     */
    public final void send(final Object message) {
        queue.add(message != null ? message : NullObject.getNullObject());
        schedule();
    }

    private void schedule() {
        if (!scheduled.compareAndSet(false, true))
            return;
        try {
            threadPool.execute(drain);
        }
        catch (RuntimeException | Error e) {
            scheduled.set(false);
            throw e;
        }
    }

    /**
     * Handles a batch of messages, then gives up the thread. A fair agent handles one message at a time.
     */
    private void drain() {
        currentThread = Thread.currentThread();
        try {
            final int batch = fair ? 1 : BATCH;
            for (int n = 0; n < batch; n++) {
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
            // nor to the next pool task, e.g. an interrupt that came after the last message
            Thread.interrupted();
            scheduled.set(false);
            if (!queue.isEmpty())
                schedule();
        }
    }

    /**
     * Dynamically dispatches the method call
     *
     * @param message A value or a closure
     */
    abstract void handleMessage(final Object message);

    /**
     * Adds the exception to the list of thrown exceptions
     *
     * @param e The exception to store
     */
    @SuppressWarnings({"MethodOnlyUsedFromInnerClass", "SynchronizedMethod"})
    synchronized void registerError(final Throwable e) {
        if (errors == null) errors = new ArrayList<Throwable>();
        errors.add(e);
    }

    /**
     * Retrieves a list of exception thrown within the agent's body.
     * Clears the exception history
     *
     * @return A detached collection of exception that have occurred in the agent's body
     */
    @SuppressWarnings({"SynchronizedMethod", "ReturnOfCollectionOrArrayField"})
    public synchronized List<Throwable> getErrors() {
        if (errors == null) return Collections.emptyList();
        try {
            return errors;
        } finally {
            errors = null;
        }
    }

    /**
     * Indicates whether there have been exception thrown within the agent's body.
     *
     * @return True, if any exceptions have occurred in the agent's body
     */
    @SuppressWarnings({"SynchronizedMethod"})
    public synchronized boolean hasErrors() {
        if (errors == null) return false;
        return !errors.isEmpty();
    }
}
