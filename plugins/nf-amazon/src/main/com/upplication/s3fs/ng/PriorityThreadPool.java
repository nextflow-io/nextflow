/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package com.upplication.s3fs.ng;

import java.util.Comparator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.RunnableFuture;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implements a thread pool executing tasks based on a priority index
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Deprecated
public class PriorityThreadPool extends ThreadPoolExecutor {

    private static Logger log = LoggerFactory.getLogger(PriorityThreadPool.class);

    public PriorityThreadPool(int corePoolSize,
                              int maximumPoolSize,
                              long keepAliveTime,
                              TimeUnit unit,
                              BlockingQueue<Runnable> workQueue,
                              ThreadFactory threadFactory,
                              RejectedExecutionHandler handler) {
        super(corePoolSize, maximumPoolSize, keepAliveTime, unit, workQueue, threadFactory, handler);
    }

    @Override
    protected <T> RunnableFuture<T> newTaskFor(Callable<T> callable) {
        RunnableFuture<T> task = super.newTaskFor(callable);
        if( callable instanceof PriorityCallable )
            return new PriorityAwareFuture<T>(task, ((PriorityCallable<T>) callable).getPriority());
        throw new IllegalArgumentException("PriorityThreadPool task must subclass PriorityCallable or PriorityRunnable class - offending task: " + callable);
    }

    protected <T> RunnableFuture<T> newTaskFor(Runnable runnable, T value) {
        RunnableFuture<T> task = super.newTaskFor(runnable, value);
        if( runnable instanceof PriorityRunnable )
            return new PriorityAwareFuture<T>(task, ((PriorityRunnable) runnable).getPriority());
        throw new IllegalArgumentException("PriorityThreadPool task must subclass PriorityCallable or PriorityRunnable class - offending task: " + runnable);
    }

    static ThreadPoolExecutor create(String name, int maxThreads, int maxQueue) {
        PriorityComparator comparator = new PriorityComparator();
        BlockingQueue<Runnable> workQueue = new PriorityBlockingQueue<>(maxQueue, comparator);
        RejectedExecutionHandler rejectPolicy =  new ThreadPoolExecutor.CallerRunsPolicy();

        ThreadPoolExecutor pool = new PriorityThreadPool(
                maxThreads,
                maxThreads,
                60L, TimeUnit.SECONDS,
                workQueue,
                CustomThreadFactory.withName(name),
                rejectPolicy );

        pool.allowCoreThreadTimeOut(true);
        log.trace("Created priority thread pool -- max-treads: {}; max-queue={}", maxThreads, maxQueue);
        return pool;
    }

    /**
     * A callable holding a priority value
     */
    static abstract class PriorityCallable<T> implements Callable<T> {

        final private int priority;
        public int getPriority() { return priority; }

        PriorityCallable(int priority) {
            this.priority = priority;
        }
    }

    static abstract class PriorityRunnable implements Runnable {

        final private int priority;
        public int getPriority() { return priority; }

        PriorityRunnable(int priority) {
            this.priority = priority;
        }
    }

    /**
     * Model a {@link RunnableFuture} adding a priority information
     *
     * @param <T>
     */
    static class PriorityAwareFuture<T> implements RunnableFuture<T> {

        private final int priority;
        private final RunnableFuture<T> target;

        public PriorityAwareFuture(RunnableFuture<T> other, int priority) {
            this.target = other;
            this.priority = priority;
        }

        public int getPriority() {
            return priority;
        }

        public boolean cancel(boolean mayInterruptIfRunning) {
            return target.cancel(mayInterruptIfRunning);
        }

        public boolean isCancelled() {
            return target.isCancelled();
        }

        public boolean isDone() {
            return target.isDone();
        }

        public T get() throws InterruptedException, ExecutionException {
            return target.get();
        }

        public T get(long timeout, TimeUnit unit) throws InterruptedException, ExecutionException, TimeoutException {
            return target.get(timeout, unit);
        }

        public void run() {
            target.run();
        }
    }

    /**
     * Compare priority of future tasks
     */
    static class PriorityComparator implements Comparator<Runnable> {
        public int compare(Runnable o1, Runnable o2) {
            if (o1 == null && o2 == null)
                return 0;
            if (o1 == null)
                return -1;
            if (o2 == null)
                return 1;
            if( o1 instanceof PriorityAwareFuture && o2 instanceof PriorityAwareFuture) {
                int p1 = ((PriorityAwareFuture<?>) o1).getPriority();
                int p2 = ((PriorityAwareFuture<?>) o2).getPriority();
                return Integer.compare(p1, p2);
            }
            // no decision
            return 0;
        }
    }
}
