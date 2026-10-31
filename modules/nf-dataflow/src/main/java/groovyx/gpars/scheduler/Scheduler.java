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

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.GParsConfig;

/**
 * The dataflow runtime: the pool shared by operators, agents and callbacks, and the engine
 * running operators.
 * <p>
 * Configured with environment variables, read once:
 * <ul>
 * <li>{@code NXF_DATAFLOW_ENGINE}: {@code mailbox} (default) or {@code legacy} (GPars actors).
 * {@code legacy} only rolls back the operator engine: agents, the shared pool and the loop thread
 * check are the same on both engines.</li>
 * <li>{@code NXF_DATAFLOW_POOL_SIZE}: a fixed number of dataflow threads. By default the pool
 * comes from the {@link GParsConfig} pool factory. With a size of 1, every operator, agent and
 * callback runs on a single loop thread, in a deterministic order.</li>
 * </ul>
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public final class Scheduler {

    public static final String ENGINE_ENV = "NXF_DATAFLOW_ENGINE";

    public static final String POOL_SIZE_ENV = "NXF_DATAFLOW_POOL_SIZE";

    private static final boolean LEGACY = parseEngine(System.getenv(ENGINE_ENV));

    private static final int POOL_SIZE = parsePoolSize(System.getenv(POOL_SIZE_ENV));

    private static final ThreadLocal<Boolean> LOOP_THREAD = new ThreadLocal<>();

    private Scheduler() {}

    private static final class Holder {
        static final Pool POOL = POOL_SIZE > 0 ? newPool(POOL_SIZE) : GParsConfig.retrieveDefaultPool();
    }

    /**
     * The shared dataflow pool
     */
    public static Pool pool() {
        return Holder.POOL;
    }

    /**
     * Whether operators run on the GPars actor engine instead of the mailbox engine
     */
    public static boolean isLegacyEngine() {
        return LEGACY;
    }

    public static String getEngineName() {
        return LEGACY ? "legacy" : "mailbox";
    }

    /**
     * The configured pool size, or 0 when the pool comes from the pool factory
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
     * Creates a fixed pool of dataflow threads. The thread of a single-thread pool is a loop thread.
     */
    static Pool newPool(final int size) {
        if (size < 1)
            throw new IllegalArgumentException("The dataflow pool size must be a positive integer: " + size);
        final boolean loop = size == 1;
        final ThreadPoolExecutor executor = new ThreadPoolExecutor(size, size, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>(), task -> {
            final Runnable body = !loop ? task : () -> {
                LOOP_THREAD.set(Boolean.TRUE);
                task.run();
            };
            final Thread thread = new Thread(body, DefaultPool.createThreadName());
            thread.setDaemon(true);
            return thread;
        });
        return new DefaultPool(executor);
    }

    static boolean parseEngine(final String value) {
        if (value == null || value.isEmpty() || "mailbox".equals(value))
            return false;
        if ("legacy".equals(value))
            return true;
        throw new IllegalArgumentException("Invalid " + ENGINE_ENV + " value: '" + value + "' -- it must be 'mailbox' or 'legacy'");
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
