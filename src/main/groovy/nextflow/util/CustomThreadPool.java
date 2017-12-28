/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.scheduler.DefaultPool;
import groovyx.gpars.scheduler.Pool;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implements a resizeable unbounded Thread pool executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CustomThreadPool extends DefaultPool {

    private static final Logger log = LoggerFactory.getLogger(CustomThreadPool.class);

    private static final long KEEP_ALIVE_TIME = 60L;

    static Pool boundedPool(int maxThreads, int queueSize) {
        log.debug(String.format("Creating bounded thread pool > maxThread: %s - queueSize: %s", maxThreads, queueSize));

        return new CustomThreadPool(

            new ThreadPoolExecutor(
                    1,
                    maxThreads,
                    KEEP_ALIVE_TIME,
                    TimeUnit.SECONDS,
                    new LinkedBlockingQueue<Runnable>(queueSize),
                    newDaemonThreadFactory(),
                    new ThreadPoolExecutor.CallerRunsPolicy())

        );

    }

    static Pool unboundedPool(int maxThreads) {
        log.debug(String.format("Creating unbounded thread pool > maxThread: %s", maxThreads));

        return new CustomThreadPool(

                new ThreadPoolExecutor(
                        1,
                        maxThreads,
                        KEEP_ALIVE_TIME,
                        TimeUnit.SECONDS,
                        new LinkedBlockingQueue<Runnable>(),
                        newDaemonThreadFactory(),
                        new ThreadPoolExecutor.CallerRunsPolicy())

        );

    }


    static Pool synchronousPool(int maxThreads) {
        log.debug(String.format("Creating synchronous thread pool > maxThread: %s", maxThreads));

        return new CustomThreadPool(

                new ThreadPoolExecutor(
                        1,
                        maxThreads,
                        KEEP_ALIVE_TIME,
                        TimeUnit.SECONDS,
                        new SynchronousQueue<Runnable>(),
                        newDaemonThreadFactory(),
                        new ThreadPoolExecutor.CallerRunsPolicy())

        );
    }

    /**
     * Creates the pool with specified executor.
     */
    private CustomThreadPool(ThreadPoolExecutor executor) {
        super(executor);
    }


    private static ThreadFactory newDaemonThreadFactory() {

        return new ThreadFactory() {
            @Override
            public Thread newThread(final Runnable task) {
                final Thread thread = new Thread(task, DefaultPool.createThreadName());
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
