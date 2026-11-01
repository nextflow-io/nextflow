/*
 * Copyright 2013-2026, Seqera Labs
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
package groovyx.gpars.dataflow.operator;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executor;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.impl.ReadCallback;

import static groovyx.gpars.dataflow.impl.ReadCallback.NONE;

/**
 * Runs an operator as a mailbox drained serially on a shared pool.
 * <p>
 * Control events are posted to the mailbox, input values and completed forks are signalled with
 * flags. Either schedules a drain on the pool unless one is already scheduled, so at most one thread
 * handles the operator at a time, like a GPars actor. The drain handles a bounded batch of work, then
 * yields the thread.
 * <p>
 * Like the GPars operator actor, the engine reads one value from each input, fires the body once all
 * have arrived, then reads the next ones. A read takes a buffered value at once, otherwise it
 * subscribes to the input, which hands over its next value. Value inputs ({@link DataflowVariable})
 * are re-read on every firing instead. With {@code maxForks > 1}, the body runs on the pool and up to
 * {@code maxForks} bodies run at once: when all are in flight, a complete input set is parked, and no
 * further input is read until a body completes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
final class MailboxEngine implements Runnable {

    private static final int BATCH = 128;

    private static final int NEW = 0;
    private static final int RUNNING = 1;
    private static final int STOPPING = 2;
    private static final int STOPPED = 3;

    private static final Object START = new Object();
    private static final Object RESUME = new Object();
    private static final Object TERMINATE = new Object();

    private final DataflowProcessor owner;

    private final Closure<?> code;

    private final int maxForks;

    private final Port[] ports;

    private final Executor pool;

    private final Queue<Object> mailbox = new ConcurrentLinkedQueue<>();

    // the number of forked bodies that completed since the last drain
    private final AtomicInteger forksDone = new AtomicInteger();

    // whether an input handed over a value since the last drain
    private volatile boolean delivered;

    private final AtomicBoolean scheduled = new AtomicBoolean();

    private final AtomicInteger state = new AtomicInteger(NEW);

    private final CountDownLatch stopped = new CountDownLatch(1);

    // the drain state, only accessed by the thread running the drain

    private final Object[] values;

    private int arrived;

    // like GPars, a firing needs a value that arrived: re-reading the value inputs alone does not fire
    private boolean fresh;

    private boolean started;

    private int inFlight;

    private List<Object> parked;

    private PoisonPill pendingPill;

    // the backpressure gate

    private volatile boolean paused;

    private volatile long pauseCount;

    // the thread running the drain, so that an external terminate can interrupt a running body

    private final Object threadLock = new Object();

    private volatile Thread drainThread;

    private boolean interrupted;

    MailboxEngine(final DataflowProcessor owner, final Executor pool, final List<DataflowReadChannel<?>> inputs, final Closure<?> code, final int maxForks) {
        this.owner = owner;
        this.pool = pool;
        this.code = code;
        this.maxForks = maxForks;
        this.values = new Object[inputs.size()];
        this.ports = new Port[inputs.size()];
        for (int i = 0; i < ports.length; i++)
            ports[i] = new Port(i, inputs.get(i));
    }

    void start() {
        if (!state.compareAndSet(NEW, RUNNING))
            throw new IllegalStateException("The operator has already been started");
        post(START);
    }

    /**
     * Stops the operator. On the drain thread it stops at once, otherwise it interrupts a running
     * drain and stops on the next one. {@code afterStop} fires once, on the drain thread.
     */
    void terminate() {
        if (!state.compareAndSet(RUNNING, STOPPING))
            return;
        if (drainThread == Thread.currentThread()) {
            finishStop();
            return;
        }
        synchronized (threadLock) {
            if (drainThread != null) {
                interrupted = true;
                drainThread.interrupt();
            }
        }
        post(TERMINATE);
    }

    void join() throws InterruptedException {
        stopped.await();
    }

    boolean isActive() {
        final int s = state.get();
        return s == RUNNING || s == STOPPING;
    }

    void resume() {
        if (paused)
            post(RESUME);
    }

    boolean isPaused() {
        return paused;
    }

    long getPauseCount() {
        return pauseCount;
    }

    private void post(final Object message) {
        if (state.get() == STOPPED)
            return;
        mailbox.offer(message);
        // an unstarted operator only queues the message: start() schedules the drain
        if (state.get() != NEW)
            schedule();
    }

    /**
     * Signals work that bypasses the mailbox: inputs and forks exist only once started
     */
    private void signal() {
        if (state.get() != STOPPED)
            schedule();
    }

    private void schedule() {
        if (!scheduled.compareAndSet(false, true))
            return;
        try {
            pool.execute(this);
        }
        catch (RuntimeException | Error e) {
            scheduled.set(false);
            throw e;
        }
    }

    /**
     * Drains a batch of work
     */
    @Override
    public void run() {
        synchronized (threadLock) {
            drainThread = Thread.currentThread();
        }
        int budget = BATCH;
        try {
            while (budget > 0) {
                final int s = state.get();
                if (s == STOPPING) {
                    // afterStop must not see the interrupt meant for a body
                    clearInterrupt();
                    finishStop();
                    break;
                }
                if (s != RUNNING)
                    break;
                try {
                    final int work = step(budget);
                    if (work == 0)
                        break;
                    budget -= work;
                }
                catch (Throwable e) {
                    try {
                        owner.reportError(e);
                    }
                    finally {
                        terminate();
                    }
                }
            }
        }
        finally {
            synchronized (threadLock) {
                drainThread = null;
            }
            // so that an interrupt meant for a body does not leak to the next pool task
            clearInterrupt();
            scheduled.set(false);
            if (state.get() == STOPPED)
                mailbox.clear();
            else if (budget <= 0 || !mailbox.isEmpty() || forksDone.get() > 0 || delivered)
                schedule();
        }
    }

    /**
     * Handles a control message, the completed forks or the delivered values, then fires on the
     * available input sets
     *
     * @return the work done, 0 if there was nothing to do
     */
    private int step(final int budget) {
        final Object message = mailbox.poll();
        if (message != null) {
            handle(message);
            return 1;
        }
        if (forksDone.get() > 0) {
            onForksDone(forksDone.getAndSet(0));
            return 1;
        }
        int work = 0;
        if (delivered) {
            delivered = false;
            collect();
            work++;
        }
        if (started)
            work += pull(budget);
        return work;
    }

    /**
     * Clears the interrupt of an external terminate, if any
     */
    private void clearInterrupt() {
        final boolean clear;
        synchronized (threadLock) {
            clear = interrupted;
            interrupted = false;
        }
        if (clear)
            Thread.interrupted();
    }

    private void handle(final Object message) {
        if (message == RESUME)
            onResume();
        else if (message == START)
            started = true;
        // TERMINATE only wakes the drain, which stops on the state
    }

    /**
     * Takes the values handed over by the inputs
     */
    private void collect() {
        for (final Port port : ports) {
            final Object item = port.item;
            if (item == NONE)
                continue;
            port.item = NONE;
            port.requested = false;
            arrive(port, item);
        }
    }

    /**
     * Reads the inputs and fires as long as complete input sets are available, up to the budget
     *
     * @return the number of firings
     */
    private int pull(final int budget) {
        int fired = 0;
        while (fired < budget && state.get() == RUNNING && !paused && parked == null && pendingPill == null) {
            for (final Port port : ports) {
                if (port.has || port.requested)
                    continue;
                if (port.value && port.read) {
                    reread(port);
                    continue;
                }
                port.requested = true;
                final Object item = port.channel.takeOrSubscribe(port);
                if (item == NONE)
                    continue;
                port.requested = false;
                arrive(port, item);
                if (pendingPill != null || state.get() != RUNNING)
                    return fired;
            }
            if (arrived < values.length || !fresh)
                break;
            fire();
            fired++;
        }
        return fired;
    }

    private void arrive(final Port port, Object item) {
        // the operator stops once the bodies in flight complete
        if (pendingPill != null)
            return;
        if (item instanceof PoisonPill) {
            item = owner.fireControlMessageArrived(port.channel, port.index, item);
            if (item instanceof PoisonPill) {
                onPoison((PoisonPill) item);
                // like GPars, the input is not read again
                port.requested = true;
                return;
            }
        }
        values[port.index] = owner.fireMessageArrived(port.channel, port.index, item);
        port.has = true;
        port.read = true;
        arrived++;
        fresh = true;
    }

    /**
     * Reads a value input again, without firing the arrival listeners
     */
    private void reread(final Port port) {
        try {
            values[port.index] = port.channel.getVal();
            port.has = true;
            arrived++;
        }
        catch (InterruptedException e) {
            throw new IllegalStateException("couldn't read the value of a DataflowVariable inside an operator.", e);
        }
    }

    private void fire() {
        final List<Object> messages = new ArrayList<>(values.length);
        for (int i = 0; i < values.length; i++) {
            messages.add(values[i]);
            values[i] = null;
            ports[i].has = false;
        }
        arrived = 0;
        fresh = false;
        if (maxForks > 1 && inFlight >= maxForks) {
            // wait for a body to complete before taking more inputs. Like GPars, beforeRun fires
            // only once the body is about to run, so that a set dropped on stop has no afterRun pending
            parked = messages;
            return;
        }
        final List<Object> verified = owner.fireBeforeRun(messages);
        if (maxForks == 1)
            runBody(verified);
        else
            fork(verified);
        afterFiring();
    }

    private void afterFiring() {
        if (state.get() == RUNNING)
            pauseIfGateClosed();
    }

    private void runBody(final List<Object> messages) {
        try {
            code.call(messages.toArray(new Object[messages.size()]));
        }
        catch (Throwable e) {
            owner.reportError(e);
        }
        finally {
            owner.fireAfterRun(messages);
        }
    }

    private void fork(final List<Object> messages) {
        inFlight++;
        try {
            pool.execute(() -> {
                try {
                    runBody(messages);
                }
                finally {
                    forksDone.incrementAndGet();
                    signal();
                }
            });
        }
        catch (RuntimeException | Error e) {
            inFlight--;
            throw e;
        }
    }

    private void onForksDone(final int count) {
        inFlight -= count;
        if (pendingPill != null) {
            if (inFlight == 0)
                forwardPoisonPill();
            return;
        }
        if (parked != null) {
            final List<Object> messages = parked;
            parked = null;
            fork(owner.fireBeforeRun(messages));
            afterFiring();
        }
    }

    /**
     * Forwards the poison pill to the outputs and stops, once the bodies in flight complete
     */
    private void onPoison(final PoisonPill pill) {
        pendingPill = pill;
        if (inFlight == 0)
            forwardPoisonPill();
    }

    private void forwardPoisonPill() {
        final PoisonPill pill = pendingPill;
        pendingPill = null;
        owner.bindAllOutputsAtomically(pill);
        owner.terminate();
    }

    private void onResume() {
        if (paused && state.get() == RUNNING) {
            paused = false;
            pauseIfGateClosed();
        }
    }

    /**
     * Pauses the operator while the gate is closed: it reads no further inputs until resumed
     */
    private void pauseIfGateClosed() {
        if (!owner.isGateClosed())
            return;
        // count the pause before publishing it, so that a paused operator is always counted.
        // It counts the pause attempts, including one undone below because the gate reopened
        pauseCount++;
        paused = true;
        // check again after publishing the flag, so that a concurrent resume() is not lost
        if (!owner.isGateClosed())
            paused = false;
    }

    private void finishStop() {
        if (!state.compareAndSet(STOPPING, STOPPED))
            return;
        mailbox.clear();
        parked = null;
        // so that a queue input hands no further value over to a stopped operator. A value handed
        // over already is dropped, as GPars does when the operator is terminated externally
        for (final Port port : ports) {
            if (port.requested)
                port.channel.unsubscribe(port);
        }
        // like GPars, forward a poison pill that waits for the bodies in flight
        final PoisonPill pill = pendingPill;
        pendingPill = null;
        try {
            if (pill != null)
                owner.bindAllOutputsAtomically(pill);
            owner.fireAfterStop();
        }
        finally {
            stopped.countDown();
        }
    }

    /**
     * An operator input: subscribed to its channel while waiting for a value, which the channel hands over
     */
    private final class Port implements ReadCallback {

        final int index;

        final DataflowReadChannel channel;

        final boolean value;

        // the value handed over by the channel, published by `delivered`
        volatile Object item = NONE;

        // waiting for a value, or closed by a control message
        boolean requested;

        // holds a value of the current input set
        boolean has;

        // a value was read once
        boolean read;

        Port(final int index, final DataflowReadChannel channel) {
            this.index = index;
            this.channel = channel;
            this.value = channel instanceof DataflowVariable;
        }

        @Override
        public void onItem(final Object item) {
            this.item = item;
            delivered = true;
            signal();
        }
    }
}
