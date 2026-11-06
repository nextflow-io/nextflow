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
package groovyx.gpars.scheduler;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;

/**
 * The dataflow runtime: the pool shared by operators, agents and callbacks.
 * <p>
 * {@code NXF_DATAFLOW_POOL_SIZE}, read once, sets a fixed number of dataflow threads. Otherwise the
 * pool comes from the executor factory, or is a default cached pool when no factory is set. With a
 * size of 1, every operator, agent and callback runs on a single loop thread, in a deterministic order.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public final class Scheduler {

    public static final String POOL_SIZE_ENV = "NXF_DATAFLOW_POOL_SIZE";

    private static final int POOL_SIZE = parsePoolSize(System.getenv(POOL_SIZE_ENV));

    private static final ThreadLocal<Boolean> LOOP_THREAD = new ThreadLocal<>();

    private static final AtomicLong THREAD_COUNT = new AtomicLong();

    private static final System.Logger log = System.getLogger(Scheduler.class.getName());

    private static volatile Supplier<ExecutorService> executorFactory;

    private static volatile boolean poolCreated;

    private Scheduler() {}

    private static final class Holder {
        static final ExecutorService POOL = createPool();
    }

    private static ExecutorService createPool() {
        poolCreated = true;
        if (POOL_SIZE > 0)
            return newPool(POOL_SIZE);
        final Supplier<ExecutorService> factory = executorFactory;
        return factory != null ? factory.get() : newDefaultPool();
    }

    /**
     * Sets the factory of the dataflow pool, used unless {@code NXF_DATAFLOW_POOL_SIZE} is set.
     * It has no effect once the pool has been created.
     */
    public static void setExecutorFactory(final Supplier<ExecutorService> factory) {
        // not an error: e.g. CmdRun installs its factory in a static initializer, which tests may load late
        if (poolCreated)
            log.log(System.Logger.Level.DEBUG, "The dataflow pool is already created -- ignoring executor factory: {0}", factory);
        executorFactory = factory;
    }

    /**
     * The shared dataflow pool
     */
    public static ExecutorService pool() {
        return Holder.POOL;
    }

    /**
     * The configured pool size, or 0 when the pool comes from the executor factory
     */
    public static int getPoolSize() {
        return POOL_SIZE;
    }

    /**
     * Whether the current thread is the loop thread of a single-thread pool: blocking on it
     * for a dataflow value would stall every operator, so it fails instead
     */
    public static boolean isLoopThread() {
        return LOOP_THREAD.get() != null;
    }

    /**
     * A JVM-unique name for a dataflow thread
     */
    public static String newThreadName() {
        return "Dataflow Thread " + THREAD_COUNT.incrementAndGet();
    }

    /**
     * Creates a fixed pool of dataflow threads. The thread of a single-thread pool is a loop thread.
     */
    static ExecutorService newPool(final int size) {
        if (size < 1)
            throw new IllegalArgumentException("The dataflow pool size must be a positive integer: " + size);
        final boolean loop = size == 1;
        return new ThreadPoolExecutor(size, size, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>(), task -> {
            final Runnable body = !loop ? task : () -> {
                LOOP_THREAD.set(Boolean.TRUE);
                task.run();
            };
            return newDaemonThread(body);
        });
    }

    /**
     * The pool used when no executor factory is set, e.g. in unit tests: a thread per task,
     * reused while warm, up to 1000 threads
     */
    static ExecutorService newDefaultPool() {
        return new ThreadPoolExecutor(1, 1000, 10L, TimeUnit.SECONDS, new SynchronousQueue<Runnable>(), Scheduler::newDaemonThread);
    }

    private static Thread newDaemonThread(final Runnable task) {
        final Thread thread = new Thread(task, newThreadName());
        thread.setDaemon(true);
        thread.setUncaughtExceptionHandler(Scheduler::logUncaught);
        return thread;
    }

    static void logUncaught(final Thread thread, final Throwable e) {
        log.log(System.Logger.Level.ERROR, "Unhandled exception in dataflow thread: " + thread.getName(), e);
    }

    static int parsePoolSize(final String value) {
        if (value == null || value.isEmpty())
            return 0;
        try {
            final int result = Integer.parseInt(value.trim());
            if (result >= 1)
                return result;
        }
        catch (NumberFormatException ignored) {
        }
        throw new IllegalArgumentException("Invalid " + POOL_SIZE_ENV + " value: '" + value + "' -- it must be a positive integer");
    }
}
