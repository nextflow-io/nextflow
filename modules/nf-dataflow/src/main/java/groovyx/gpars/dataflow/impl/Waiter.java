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
package groovyx.gpars.dataflow.impl;

import java.util.concurrent.locks.LockSupport;

import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.expression.DataflowExpression;
import groovyx.gpars.scheduler.Scheduler;

/**
 * A thread blocked in a channel read
 * <p>
 * It parks instead of waiting on a monitor, so that a blocked virtual thread does not pin its carrier.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public final class Waiter implements ReadCallback {

    private final Thread thread = Thread.currentThread();

    private volatile Object item = NONE;

    private Waiter() {}

    @Override
    public void onItem(final Object item) {
        this.item = item;
        LockSupport.unpark(thread);
    }

    /**
     * Reads the next item of a channel, blocking until there is one
     *
     * @param timeoutNanos the maximum time to wait, or a negative value to wait without a timeout
     * @return the item, or {@link ReadCallback#NONE} when the timeout expires
     */
    public static Object read(final DataflowReadChannel<?> channel, final long timeoutNanos) throws InterruptedException {
        final Waiter waiter = new Waiter();
        final Object result = channel.takeOrSubscribe(waiter);
        if (result != NONE)
            return result;
        if (Scheduler.isLoopThread() && channel.unsubscribe(waiter))
            throw new IllegalStateException("Cannot block on an unbound " + channel.getClass().getSimpleName() + " on the dataflow loop thread -- blocking reads are not allowed when the dataflow pool has a single thread (" + Scheduler.POOL_SIZE_ENV + "=1)");
        return waiter.await(channel, timeoutNanos);
    }

    private Object await(final DataflowReadChannel<?> channel, final long timeoutNanos) throws InterruptedException {
        final long deadline = System.nanoTime() + timeoutNanos;
        Object result;
        while ((result = item) == NONE) {
            if (Thread.interrupted()) {
                if (channel.unsubscribe(this))
                    throw new InterruptedException();
                // the item is being handed over: take it, so that it is not lost, and keep the interrupt
                return handover(channel, true);
            }
            if (timeoutNanos < 0) {
                LockSupport.park(this);
                continue;
            }
            final long remaining = deadline - System.nanoTime();
            if (remaining <= 0)
                return channel.unsubscribe(this) ? NONE : handover(channel, false);
            LockSupport.parkNanos(this, remaining);
        }
        return result;
    }

    /**
     * Waits for an item that the channel is handing over, then restores the interrupt status if requested or interrupted meanwhile
     */
    private Object handover(final DataflowReadChannel<?> channel, boolean interrupted) {
        Object result;
        while ((result = item) == NONE) {
            // a value channel is bound before its waiters are called
            if (channel instanceof DataflowExpression && channel.isBound()) {
                result = channel.takeOrSubscribe(this);
                break;
            }
            LockSupport.park(this);
            interrupted |= Thread.interrupted();
        }
        if (interrupted)
            Thread.currentThread().interrupt();
        return result;
    }
}
