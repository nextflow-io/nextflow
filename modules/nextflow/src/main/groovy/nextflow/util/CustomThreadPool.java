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

package nextflow.util;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.scheduler.Scheduler;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Factories of the dataflow thread pools
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public final class CustomThreadPool {

    private static final Logger log = LoggerFactory.getLogger(CustomThreadPool.class);

    private static final long KEEP_ALIVE_TIME = 60L;

    private CustomThreadPool() {}

    static ExecutorService boundedPool(int maxThreads, int queueSize) {
        log.debug(String.format("Creating bounded thread pool > maxThread: %s - queueSize: %s", maxThreads, queueSize));

        return new ThreadPoolExecutor(
                1,
                maxThreads,
                KEEP_ALIVE_TIME,
                TimeUnit.SECONDS,
                new LinkedBlockingQueue<Runnable>(queueSize),
                newDaemonThreadFactory(),
                new ThreadPoolExecutor.CallerRunsPolicy());
    }

    static ExecutorService unboundedPool(int maxThreads) {
        log.debug(String.format("Creating unbounded thread pool > maxThread: %s", maxThreads));

        return new ThreadPoolExecutor(
                1,
                maxThreads,
                KEEP_ALIVE_TIME,
                TimeUnit.SECONDS,
                new LinkedBlockingQueue<Runnable>(),
                newDaemonThreadFactory(),
                new ThreadPoolExecutor.CallerRunsPolicy());
    }

    static ExecutorService defaultPool(int poolSize, int maxThreads) {
        log.debug(String.format("Creating default thread pool > poolSize: %s; maxThreads: %s", poolSize, maxThreads));

        return new ThreadPoolExecutor(
                poolSize,
                maxThreads,
                KEEP_ALIVE_TIME,
                TimeUnit.SECONDS,
                new SynchronousQueue<Runnable>(),
                newDaemonThreadFactory(),
                new ThreadPoolExecutor.CallerRunsPolicy());
    }

    static ExecutorService synchronousPool(int maxThreads) {
        log.debug(String.format("Creating synchronous thread pool > maxThread: %s", maxThreads));

        return new ThreadPoolExecutor(
                1,
                maxThreads,
                KEEP_ALIVE_TIME,
                TimeUnit.SECONDS,
                new SynchronousQueue<Runnable>(),
                newDaemonThreadFactory(),
                new ThreadPoolExecutor.CallerRunsPolicy());
    }

    /**
     * The pool of GPars when no pool factory is set: one core thread, up to 1000 threads
     */
    static ExecutorService legacyPool() {
        log.debug("Creating legacy thread pool");

        return new ThreadPoolExecutor(
                1,
                1000,
                10L,
                TimeUnit.SECONDS,
                new SynchronousQueue<Runnable>(),
                newDaemonThreadFactory());
    }

    // a platform thread per task, not virtual threads, with the daemon factory as since #3871 (5c429046f)
    static ExecutorService virtualPool() {
        return Executors.newThreadPerTaskExecutor(newDaemonThreadFactory());
    }

    private static ThreadFactory newDaemonThreadFactory() {

        return new ThreadFactory() {
            @Override
            public Thread newThread(final Runnable task) {
                final Thread thread = new Thread(task, Scheduler.newThreadName());
                thread.setDaemon(true);
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
