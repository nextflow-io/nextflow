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
package groovyx.gpars;

import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.function.BooleanSupplier;

import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import groovyx.gpars.dataflow.impl.ReadCallback;

/**
 * Thread helpers for concurrency tests.
 */
public final class TestThreads {

    private static final long DEADLINE_MILLIS = 5_000;

    private TestThreads() {}

    @FunctionalInterface
    public interface Task { void run() throws Throwable; }

    /**
     * Start a thread running the given task. Errors are rethrown on the thread, so they are
     * reported to its uncaught exception handler.
     */
    public static Thread start(Task task) {
        final Thread thread = new Thread(() -> {
            try {
                task.run();
            }
            catch( RuntimeException | Error e ) {
                throw e;
            }
            catch( Throwable e ) {
                throw new RuntimeException(e);
            }
        });
        thread.setDaemon(true);
        thread.start();
        return thread;
    }

    /**
     * Wait until a condition holds, failing after a deadline.
     */
    public static void await(BooleanSupplier condition) {
        final long deadline = System.currentTimeMillis() + DEADLINE_MILLIS;
        while( !condition.getAsBoolean() ) {
            if( System.currentTimeMillis() > deadline )
                throw new AssertionError("Condition not met in time");
            try {
                Thread.sleep(1);
            }
            catch( InterruptedException e ) {
                Thread.currentThread().interrupt();
                throw new AssertionError("Interrupted while waiting for a condition", e);
            }
        }
    }

    /**
     * Read the next value of a channel, like {@code getVal()}, waiting up to the timeout.
     *
     * @return the value, or null when the timeout expires
     */
    @SuppressWarnings("unchecked")
    public static <T> T poll(DataflowReadChannel<T> channel, long timeout, TimeUnit unit) throws InterruptedException {
        final CompletableFuture<Object> future = new CompletableFuture<>();
        final ReadCallback callback = future::complete;
        Object item = channel.takeOrSubscribe(callback);
        if( item == ReadCallback.NONE ) {
            try {
                item = future.get(Math.max(0, timeout), unit);
            }
            catch( TimeoutException e ) {
                if( channel.unsubscribe(callback) )
                    return null;
                // the value is being handed over
                item = future.join();
            }
            catch( ExecutionException e ) {
                throw new AssertionError(e);
            }
        }
        // a value channel reads as null when bound to an error
        return channel instanceof DataflowVariable ? ((DataflowVariable<T>) channel).getVal() : (T) item;
    }

    /**
     * Wait until a thread is parked, e.g. blocked reading an unbound channel.
     */
    public static void awaitBlocked(Thread thread) {
        await(() -> {
            final Thread.State state = thread.getState();
            return state == Thread.State.WAITING || state == Thread.State.TIMED_WAITING;
        });
    }
}
