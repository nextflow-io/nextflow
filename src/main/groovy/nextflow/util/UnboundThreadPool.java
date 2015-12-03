package nextflow.util;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.scheduler.DefaultPool;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implements a resizeable unbounded Thread pool executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class UnboundThreadPool extends DefaultPool {

    private static final Logger log = LoggerFactory.getLogger(UnboundThreadPool.class);

    private static final long KEEP_ALIVE_TIME = 60L;

    /**
     * Creates the pool with specified number of threads.
     *
     * @param daemon   Sets the daemon flag of threads in the pool.
     * @param poolSize The required size of the pool
     */
    public UnboundThreadPool(final boolean daemon, final int poolSize) {
        super(createResizeablePool(daemon, poolSize));
    }

    /**
     * Creates a fixed-thread pool of given size. Each thread will have the uncaught exception handler set
     * to print the unhandled exception to standard error output.
     *
     * @param daemon   Sets the daemon flag of threads in the pool.
     * @param poolSize The required pool size  @return The created thread pool
     * @return The newly created thread pool
     */
    private static ThreadPoolExecutor createResizeablePool(final boolean daemon, final int poolSize) {
        assert poolSize > 0;
        log.debug("Creating gpars pool with size: "+poolSize);
        return new ThreadPoolExecutor(
                1,
                poolSize,
                KEEP_ALIVE_TIME,
                TimeUnit.SECONDS,
                new LinkedBlockingQueue<Runnable>(),
                newThreadFactory(daemon),
                new ThreadPoolExecutor.CallerRunsPolicy());
    }

    private static ThreadFactory newThreadFactory(final boolean daemon) {

        return new ThreadFactory() {
            @Override
            public Thread newThread(final Runnable task) {
                final Thread thread = new Thread(task, DefaultPool.createThreadName());
                thread.setDaemon(daemon);
                thread.setUncaughtExceptionHandler(new Thread.UncaughtExceptionHandler() {
                    @Override
                    public void uncaughtException(final Thread t, final Throwable e) {
                        log.error("Unexpected thread exception", e);
                    }
                });
                return thread;
            }
        };

    }
}
