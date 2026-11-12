// GPars - Groovy Parallel Systems
//
// Copyright © 2008-10  The original author or authors
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
// Modified by Seqera Labs for Nextflow (2026): removed remote message stream support (groovyx.gpars.remote, groovyx.gpars.serial).
package groovyx.gpars.actor.impl;

import groovy.time.Duration;
import groovyx.gpars.actor.Actor;
import groovyx.gpars.actor.ActorMessage;

import java.io.Serializable;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.LockSupport;

/**
 * Represents a stream of messages and forms the base class for actors
 *
 * @author Alex Tkachman, Vaclav Pech, Dierk Koenig
 */
public abstract class MessageStream implements Serializable {
    private static final long serialVersionUID = 7644465423857532477L;

    /**
     * Send message to stream and return immediately
     *
     * @param message message to send
     * @return always return message stream itself
     */
    public abstract MessageStream send(Object message);

    /**
     * Convenience method for send(new Object()).
     *
     * @return always return message stream itself
     */
    public final MessageStream send() {
        return send(new Object());
    }

    /**
     * Send message to stream and return immediately. Allows to specify an arbitrary actor to send replies to.
     * By default replies are sent to the originator (sender) of each message, however, when a different actor
     * is specified as the optional second argument to the send() method, this supplied actor will receive the replies instead.
     *
     * @param message message to send
     * @param replyTo where to send reply
     * @param <T>     type of message accepted by the stream
     * @return always return message stream itself
     */
    public final <T> MessageStream send(final T message, final MessageStream replyTo) {
        return send(new ActorMessage(message, replyTo));
    }

    /**
     * Same as send
     *
     * @param message to send
     * @return original stream
     */
    public final <T> MessageStream leftShift(final T message) {
        return send(message);
    }

    /**
     * Same as send
     *
     * @param message to send
     * @return original stream
     */
    public final <T> MessageStream call(final T message) {
        return send(message);
    }

    /**
     * Sends a message and waits for a reply.
     * Returns the reply or throws an IllegalStateException, if the target actor cannot reply.
     *
     * @param message message to send
     * @return The message that came in reply to the original send.
     * @throws InterruptedException if interrupted while waiting
     */
    public final <T, V> V sendAndWait(final T message) throws InterruptedException {
        final ResultWaiter<V> to = new ResultWaiter<V>();
        send(new ActorMessage(message, to));
        return to.getResult();
    }

    /**
     * Sends a message and waits for a reply. Timeouts after the specified timeout. In case of timeout returns null.
     * Returns the reply or throws an IllegalStateException, if the target actor cannot reply.
     *
     * @param message message to send
     * @param timeout timeout
     * @param units   units
     * @return The message that came in reply to the original send.
     * @throws InterruptedException if interrupted while waiting
     */
    public final <T> Object sendAndWait(final T message, final long timeout, final TimeUnit units) throws InterruptedException {
        final ResultWaiter<Object> to = new ResultWaiter<Object>();
        send(new ActorMessage(message, to));
        return to.getResult(timeout, units);
    }

    /**
     * Sends a message and waits for a reply. Timeouts after the specified timeout. In case of timeout returns null.
     * Returns the reply or throws an IllegalStateException, if the target actor cannot reply.
     *
     * @param message  message to send
     * @param duration timeout
     * @return The message that came in reply to the original send.
     * @throws InterruptedException if interrupted while waiting
     */
    public final <T> Object sendAndWait(final T message, final Duration duration) throws InterruptedException {
        return sendAndWait(message, duration.toMilliseconds(), TimeUnit.MILLISECONDS);
    }

    static void reInterrupt() throws InterruptedException {
        if (Thread.currentThread().isInterrupted()) {
            throw new InterruptedException();
        }
    }

    /**
     * Represents a pending request for a reply from an actor.
     *
     * @param <V> The type of expected reply message
     */
    private static class ResultWaiter<V> extends MessageStream {
        private static final long serialVersionUID = 6512046150496489148L;

        /**
         * Holds a reference to the calling thread, while waiting, and the received reply message, once it has arrived.
         */
        private volatile Object value;

        private volatile boolean isSet;

        private ResultWaiter() {
            value = Thread.currentThread();
        }

        /**
         * Accepts the message as a reply and wakes up the sleeping thread.
         *
         * @param message message to send
         * @return this
         */
        @Override
        public MessageStream send(final Object message) {
            if (isSet)
                throw new IllegalStateException("A reply has already been sent. The originator does not expect more than one reply.");
            final Thread thread = (Thread) this.value;
            if (message instanceof ActorMessage) {
                this.value = ((ActorMessage) message).getPayLoad();
            } else {
                this.value = message;
            }
            isSet = true;
            LockSupport.unpark(thread);
            return this;
        }

        /**
         * Retrieves the response blocking until a message arrives
         *
         * @return The received message
         * @throws InterruptedException If the thread gets interrupted
         */
        @SuppressWarnings("unchecked")
        public V getResult() throws InterruptedException {
            while (!isSet) {
                LockSupport.park();
                final Thread thread = Thread.currentThread();
                if (thread.isInterrupted()) {
                    throw new InterruptedException();
                }
            }
            rethrowException();
            return (V) value;
        }

        /**
         * Retrieves the response blocking until a message arrives
         *
         * @param timeout How long to wait
         * @param units   Unit for the timeout
         * @return The received message
         * @throws InterruptedException If the thread gets interrupted
         */
        public Object getResult(final long timeout, final TimeUnit units) throws InterruptedException {
            final long endNano = System.nanoTime() + units.toNanos(timeout);
            while (!isSet) {
                final long toWait = endNano - System.nanoTime();
                if (toWait <= 0L) {
                    return null;
                }
                LockSupport.parkNanos(toWait);
                MessageStream.reInterrupt();
            }
            rethrowException();
            return value;
        }

        private void rethrowException() {
            if (value instanceof Throwable) {
                if (value instanceof RuntimeException) {
                    throw (RuntimeException) value;
                } else {
                    throw new RuntimeException((Throwable) value);
                }
            }
        }

        /**
         * Handle cases when the message sent to the actor doesn't get delivered
         *
         * @param msg The message that failed to get delivered
         */
        @SuppressWarnings("unused") //  TODO:  Eclipse requires this to be tagged as unused.
        public void onDeliveryError(final Object msg) {
            send(new IllegalStateException("Delivery error. Maybe target actor is not active"));
        }
    }
}
