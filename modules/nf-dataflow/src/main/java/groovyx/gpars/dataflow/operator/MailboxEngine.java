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
import java.util.Arrays;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import groovy.lang.Closure;
import groovyx.gpars.actor.impl.MessageStream;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.group.PGroup;
import groovyx.gpars.scheduler.Pool;

/**
 * Runs an operator as a mailbox drained serially on a shared pool.
 * <p>
 * Input arrivals and control events are posted to the mailbox. A post schedules a drain on the
 * pool unless one is already scheduled, so at most one thread handles the operator at a time,
 * like a GPars actor. The drain handles a bounded batch of messages, then yields the thread.
 * <p>
 * Like the GPars operator actor, the engine requests one value from each input, fires the body
 * once all have arrived, then requests the next ones. Value inputs ({@link DataflowVariable})
 * are re-read on every firing instead of being requested. With {@code maxForks > 1}, the body runs
 * on the pool and up to {@code maxForks} bodies run at once: when all are in flight, a complete
 * input set is parked, and no further input is requested until a body completes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
final class MailboxEngine extends OperatorEngine implements Runnable {

    private static final int BATCH = 128;

    private static final int NEW = 0;
    private static final int RUNNING = 1;
    private static final int STOPPING = 2;
    private static final int STOPPED = 3;

    private static final Object START = new Object();
    private static final Object RESUME = new Object();
    private static final Object FORK_DONE = new Object();
    private static final Object TERMINATE = new Object();
    private static final Object STOP_GENTLY = new Object();

    private final DataflowProcessor owner;

    private final List<DataflowReadChannel<?>> inputs;

    private final Closure<?> code;

    private final int maxForks;

    private final Port[] ports;

    private volatile Pool pool;

    private final Queue<Object> mailbox = new ConcurrentLinkedQueue<>();

    private final AtomicBoolean scheduled = new AtomicBoolean();

    private final AtomicInteger state = new AtomicInteger(NEW);

    private final CountDownLatch stopped = new CountDownLatch(1);

    // the drain state, only accessed by the thread running the drain

    private final Object[] values;

    private int arrived;

    private int inFlight;

    private List<Object> parked;

    private PoisonPill pendingPill;

    private boolean stoppingGently;

    // the backpressure gate

    private volatile boolean paused;

    private volatile long pauseCount;

    // the thread running the drain, so that an external terminate can interrupt a running body

    private final Object threadLock = new Object();

    private volatile Thread drainThread;

    private boolean interrupted;

    MailboxEngine(final DataflowProcessor owner, final Pool pool, final List<DataflowReadChannel<?>> inputs, final Closure<?> code, final int maxForks) {
        this.owner = owner;
        this.pool = pool;
        this.inputs = inputs;
        this.code = code;
        this.maxForks = maxForks;
        this.values = new Object[inputs.size()];
        this.ports = new Port[inputs.size()];
        for (int i = 0; i < ports.length; i++)
            ports[i] = new Port(i, inputs.get(i) instanceof DataflowVariable);
    }

    @Override
    void start(final PGroup group) {
        if (group != null)
            pool = group.getThreadPool();
        if (!state.compareAndSet(NEW, RUNNING))
            throw new IllegalStateException("The operator has already been started");
        post(START);
    }

    /**
     * Stops the operator. On the drain thread it stops at once, otherwise it interrupts a running
     * drain and stops on the next one. {@code afterStop} fires once, on the drain thread.
     */
    @Override
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

    @Override
    void terminateAfterNextRun() {
        post(STOP_GENTLY);
    }

    @Override
    void join() throws InterruptedException {
        stopped.await();
    }

    @Override
    boolean isActive() {
        final int s = state.get();
        return s == RUNNING || s == STOPPING;
    }

    @Override
    void resume() {
        if (paused)
            post(RESUME);
    }

    @Override
    boolean isPaused() {
        return paused;
    }

    @Override
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
     * Drains a batch of messages
     */
    @Override
    public void run() {
        synchronized (threadLock) {
            drainThread = Thread.currentThread();
        }
        try {
            for (int n = 0; n < BATCH; n++) {
                final int s = state.get();
                if (s == STOPPING) {
                    // afterStop must not see the interrupt meant for a body
                    clearInterrupt();
                    finishStop();
                    break;
                }
                if (s != RUNNING)
                    break;
                final Object message = mailbox.poll();
                if (message == null)
                    break;
                try {
                    handle(message);
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
            else if (!mailbox.isEmpty())
                schedule();
        }
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
        if (message instanceof Port)
            onArrival((Port) message);
        else if (message == FORK_DONE)
            onForkDone();
        else if (message == RESUME)
            onResume();
        else if (message == START)
            onStart();
        else if (message == STOP_GENTLY)
            stoppingGently = true;
        // TERMINATE only wakes the drain, which stops on the state
    }

    private void onStart() {
        owner.fireAfterStart();
        requestInputs(true);
    }

    private void onArrival(final Port port) {
        final int index = port.index;
        Object result = port.take();
        // the operator stops once the bodies in flight complete
        if (pendingPill != null)
            return;
        final DataflowReadChannel channel = inputs.get(index);
        if (result instanceof ControlMessage) {
            result = owner.fireControlMessageArrived(channel, index, result);
            if (result instanceof PoisonPill)
                onPoison((PoisonPill) result);
            if (result instanceof ControlMessage)
                return;
        }
        values[index] = owner.fireMessageArrived(channel, index, result);
        if (++arrived < values.length)
            return;
        final List<Object> messages = new ArrayList<>(Arrays.asList(values));
        Arrays.fill(values, null);
        arrived = 0;
        fire(owner.fireBeforeRun(messages));
    }

    private void fire(final List<Object> messages) {
        if (maxForks == 1)
            runBody(messages);
        else if (inFlight < maxForks)
            fork(messages);
        else {
            // wait for a body to complete before taking more inputs
            parked = messages;
            return;
        }
        afterFiring();
    }

    private void afterFiring() {
        if (stoppingGently) {
            owner.terminate();
            return;
        }
        if (state.get() == RUNNING)
            requestInputsUnlessPaused();
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
                    post(FORK_DONE);
                }
            });
        }
        catch (RuntimeException | Error e) {
            inFlight--;
            throw e;
        }
    }

    private void onForkDone() {
        inFlight--;
        if (pendingPill != null) {
            if (inFlight == 0)
                forwardPoisonPill();
            return;
        }
        if (parked != null) {
            final List<Object> messages = parked;
            parked = null;
            fork(messages);
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
        pill.countDown();
    }

    private void onResume() {
        if (paused && state.get() == RUNNING) {
            paused = false;
            requestInputsUnlessPaused();
        }
    }

    /**
     * Requests the next input set, unless the gate is closed: then the operator is paused until resumed
     */
    private void requestInputsUnlessPaused() {
        if (owner.isGateClosed()) {
            // count the pause before publishing it, so that a paused operator is always counted.
            // It counts the pause attempts, including one undone below because the gate reopened
            pauseCount++;
            paused = true;
            // check again after publishing the flag, so that a concurrent resume() is not lost
            if (owner.isGateClosed())
                return;
            paused = false;
        }
        requestInputs(false);
    }

    private void requestInputs(final boolean initial) {
        for (final Port port : ports) {
            if (initial || !port.value)
                port.request();
            else
                port.reread();
        }
    }

    private void finishStop() {
        if (!state.compareAndSet(STOPPING, STOPPED))
            return;
        mailbox.clear();
        parked = null;
        // like the legacy engine, forward a poison pill that waits for the bodies in flight
        final PoisonPill pill = pendingPill;
        pendingPill = null;
        try {
            if (pill != null)
                owner.bindAllOutputsAtomically(pill);
            owner.fireAfterStop();
        }
        finally {
            if (pill != null)
                pill.countDown();
            stopped.countDown();
        }
    }

    /**
     * An operator input: receives one requested value at a time and posts itself to the mailbox
     */
    private final class Port extends MessageStream {

        private static final long serialVersionUID = 1L;

        final int index;

        final boolean value;

        // published to the drain by the mailbox
        private Object message;

        Port(final int index, final boolean value) {
            this.index = index;
            this.value = value;
        }

        void request() {
            inputs.get(index).getValAsync(this);
        }

        /**
         * Reads a value input again, without firing the arrival listeners
         */
        void reread() {
            try {
                values[index] = inputs.get(index).getVal();
                arrived++;
            }
            catch (InterruptedException e) {
                throw new IllegalStateException("couldn't read the value of a DataflowVariable inside an operator.", e);
            }
        }

        Object take() {
            final Object result = message;
            message = null;
            return result;
        }

        @Override
        public MessageStream send(final Object message) {
            this.message = message;
            post(this);
            return this;
        }
    }
}
