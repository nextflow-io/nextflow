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

package nextflow.util

import java.util.concurrent.Callable
import java.util.concurrent.FutureTask
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.RejectedExecutionException
import java.util.concurrent.RunnableFuture
import java.util.concurrent.Semaphore
import java.util.concurrent.ThreadFactory
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


/**
 * Implements the throttling and retry logic
 */
@Slf4j
@CompileStatic
class ThrottlingExecutor extends ThreadPoolExecutor {

    /**
     * A customised thread factory 
     */
    @CompileStatic
    static private class ThreadFactory0 implements ThreadFactory {

        private ThreadGroup group
        private AtomicInteger threadNumber = new AtomicInteger(1)
        private ExceptionHandler0 handler = new ExceptionHandler0();
        private prefix

        ThreadFactory0(String prefix) {
            this.prefix = prefix ?: ThrottlingExecutor.simpleName
            this.group = System.getSecurityManager()?.getThreadGroup() ?: Thread.currentThread().getThreadGroup()
        }

        Thread newThread(Runnable r) {
            def thread = new Thread(group, r, "${prefix}-${threadNumber.getAndIncrement()}", 0)
            if (thread.isDaemon())
                thread.setDaemon(false);
            if (thread.getPriority() != Thread.NORM_PRIORITY)
                thread.setPriority(Thread.NORM_PRIORITY);
            thread.setUncaughtExceptionHandler(handler)
            return thread
        }
    }

    /**
     * Redirect uncaught exceptionHandler to logging subsystem
     */
    @Slf4j
    @CompileStatic
    static private class ExceptionHandler0 implements Thread.UncaughtExceptionHandler {
        @Override
        void uncaughtException(Thread t, Throwable e) {
            log.error("Unexpected error task submission error", e)
        }
    }

    /**
     * Throttling options
     */
    @ToString(includeNames = true, excludes='successAction,failureAction,retryAction,rateLimitChangeAction,retryCondition')
    static class Options {

        String name

        RateLimiter limiter

        Closure successAction

        Closure failureAction

        Closure retryAction

        Closure rateLimitChangeAction

        Closure<Boolean> retryCondition

        int poolSize

        int maxPoolSize

        int queueSize

        int maxRetries = Integer.MAX_VALUE

        Duration keepAlive

        boolean autoThrottle

        long errorBurstDelayMillis = 1_000

        int rampUpInterval = 100

        float rampUpFactor = 1.2

        double rampUpMaxRate = Double.MAX_VALUE

        double backOffMinRate

        float backOffFactor = 2.0f

        Duration retryDelay = Duration.of('1 sec')

        Options normalise() {
            if( !poolSize )
                poolSize = Runtime.runtime.availableProcessors()
            if( !maxPoolSize || poolSize>maxPoolSize )
                maxPoolSize = poolSize
            if( keepAlive==null )
                keepAlive = Duration.of(0)
            if( !name )
                name = ThrottlingExecutor.simpleName

            return this
        }

        Options withName(String n) {
            this.name = n
            return this
        }

        Options withAutoThrottle(boolean value=true) {
            autoThrottle = value
            return this
        }

        Options withKeepAlive(Duration d)  {
            keepAlive = d
            return this
        }

        Options withPoolSize(int size) {
            poolSize = size
            return this
        }

        Options withMaxPoolSize(int size) {
            maxPoolSize = size
            return this
        }

        Options withQueueSize(int size ) {
            queueSize = size
            return this
        }

        Options withRateLimit(RateLimiter rate ) {
            limiter = rate
            return this
        }

        Options withRateLimit(String rate ) {
            limiter = new RateUnit(rate).getRateLimiter()
            return this
        }

        Options withMaxRetries(int max) {
            maxRetries = max
            return this
        }

        Options retryOn(Class<? extends Throwable> clazz) {
            retryCondition = { failure -> clazz.isAssignableFrom(clazz) }
            return this
        }

        Options retryOn(Closure<Boolean> condition) {
            retryCondition = condition
            return this
        }

        Options withRampUp(double maxRate, Float factor=null, Integer interval=null) {
            this.rampUpMaxRate = maxRate
            if( factor!=null )
                this.rampUpFactor = factor
            if( interval!=null )
                this.rampUpInterval = interval
            return this
        }

        Options withBackOff(double minRate, Float factor=null) {
            backOffMinRate = minRate
            if( factor != null )
                backOffFactor = factor
            return this
        }

        Options withBackOff(String minRate, float factor) {
            withBackOff(new RateUnit(minRate).rate, factor)
        }

        Options withErrorBurstDelay(Duration duration) {
            this.errorBurstDelayMillis = duration.millis
            return this
        }

        Options withRetryDelay(Duration d) {
            this.retryDelay = d
            return this
        }

        Options onSuccess(Closure action) {
            successAction = action
            return this
        }

        Options onFailure( @ClosureParams(value = SimpleType.class, options = "java.lang.Throwable") Closure action) {
            failureAction = action
            return this
        }

        Options onRetry( @ClosureParams(value = SimpleType.class, options = "java.lang.Throwable") Closure action) {
            retryAction = action
            return this
        }

        Options onRateLimitChange( @ClosureParams(value = SimpleType.class, options = "nextflow.util.RateUnit") Closure action) {
            rateLimitChangeAction = action
            return this
        }

        Options withOptions(Map opts) {
            if( opts == null )
                opts = Collections.emptyMap()

            if( opts.rateLimit )
                limiter = RateLimiter.create( (opts.rateLimit as RateUnit).rate  )

            if( opts.poolSize )
                poolSize = opts.poolSize as int

            if( opts.maxPoolSize )
                maxPoolSize = opts.maxPoolSize as int

            if( opts.queueSize )
                queueSize = opts.queueSize as int

            if( opts.maxRetries )
                maxRetries = opts.maxRetries as int

            if( opts.keepAlive )
                keepAlive = opts.keepAlive as Duration

            if( opts.autoThrottle )
                autoThrottle = opts.autoThrottle as boolean

            if( opts.errorBurstDelay )
                errorBurstDelayMillis = (opts.errorBurstDelay as Duration).millis

            if( opts.rampUpInterval )
                rampUpInterval = opts.rampUpInterval as int

            if( opts.rampUpFactor)
                rampUpFactor = opts.rampUpFactor as float

            if( opts.rampUpMaxRate)
                rampUpMaxRate = opts.rampUpMaxRate as double

            if( opts.backOffMinRate )
                backOffMinRate = (opts.backOffMinRate as RateUnit).rate

            if( opts.backOffFactor )
                backOffFactor = opts.backOffFactor as float

            if( opts.retryDelay )
                retryDelay = opts.retryDelay as Duration

            return this
        }
    }

    private abstract class Recoverable implements Callable, Runnable {

        abstract Object invoke()

        @Override
        final Object call() throws Exception {
            while(true)
                try {
                    final result = invoke()
                    handleSuccess0(result)
                    return result
                }
                catch (Throwable t) {
                    if( !handleRetry0(t) )
                        return null
                }
        }

        @Override
        final void run() { call() }
    }

    // ---- executor settings ----

    private Options opts

    private Semaphore semaphore

    private AtomicLong successCount = new AtomicLong()

    private AtomicInteger retryCount = new AtomicInteger()

    private long lastFailureMillis

    /**
     * Executor factory method
     *
     * @param opts Throttling and executor options
     * @return The {@link ThrottlingExecutor} object instance 
     */
    static ThrottlingExecutor create( Options opts ) {
        assert opts
        opts.normalise()
        log.debug "Creating throttling executor with opts: $opts"
        new ThrottlingExecutor(opts)
    }

    protected ThrottlingExecutor(Options opts) {
        super(
                opts.poolSize,
                opts.maxPoolSize,
                opts.keepAlive.millis, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new ThreadFactory0(opts.name))
        // the semaphore is bounding both the number of tasks currently executing
        // and those queued up
        this.opts = opts
        if( opts.queueSize )
            this.semaphore = new Semaphore(opts.maxPoolSize + opts.queueSize)
    }

    @Override
    void execute(final Runnable task) {
        def wrapper = task instanceof FutureTask ? task : wrap(task)

        acquireSemaphore()
        log.trace "Executing task"

        try {
            super.execute(wrapper)
        }
        catch (RejectedExecutionException e) {
            semaphore?.release()
            throw e
        }
    }

    protected <T> Callable<T> wrap(Callable<T> task) {
        if( task instanceof Recoverable )
            return task

        new Recoverable() {
            @Override Object invoke() { return task.call() }
        }
    }

    protected  Runnable wrap(Runnable task) {
        if( task instanceof Recoverable )
            return task

        new Recoverable() {
            @Override Object invoke() { return task.run() }
        }
    }

    @Override
    protected <T> RunnableFuture<T> newTaskFor(Callable<T> task) {
        return super.newTaskFor(wrap(task))
    }

    @Override
    protected <T> RunnableFuture<T> newTaskFor(Runnable task, T value) {
        return super.newTaskFor(wrap(task), value)
    }


    protected void beforeExecute(Thread thread, Runnable task) {
        log.trace "Before acquire limit=$opts.limiter"
        if( opts.limiter ) opts.limiter.acquire()
        super.beforeExecute(thread,task)
        log.trace "Before execute task"
    }

    @Override
    protected void afterExecute(Runnable task, Throwable ex) {
        super.afterExecute(task,ex)
        semaphore?.release()
        log.trace "After execute task"
    }

    private void acquireSemaphore() {
        boolean acquired = false
        while (!acquired) {
            try {
                semaphore?.acquire()
                acquired = true
            }
            catch (InterruptedException e) {
                log.warn("InterruptedException whilst acquiring semaphore", e);
            }
        }
    }


    private void handleSuccess0(Object result) {
        log.trace "Succeed task"

        // invoke the success handler
        opts.successAction?.call(result)
        // reset the retry counter on the first successful execution
        retryCount.set(0)
        // try to speed-up the submission rate
        final count = successCount.incrementAndGet()
        if( opts.autoThrottle && opts.rampUpInterval>0 && (count % opts.rampUpInterval)==0 ) {
            speedUp0(opts.limiter)
        }
    }

    private boolean handleRetry0(Throwable t) {
        log.trace "Failure exception=$t"

        if( opts.retryCondition && opts.retryCondition.call(t) && retryCount.getAndIncrement()<opts.maxRetries ) {
            // slow down the submission throughput
            if( opts.autoThrottle ) slowDown0()
            // now retry
            retry(t)
            return true
        }
        else if( opts.failureAction ) {
            opts.failureAction.call(t)
        }
        else {
            log.error("Unable to process client request", t)
        }
        return false
    }

    private void retry(Throwable t) {
        if( isTerminating() ) {
            log.debug("WARN: Task execution cannot be retried -- executor is terminating")
            return
        }

        // sleep before retry
        try {
            sleep(opts.retryDelay.millis * retryCount.get())
        }
        catch( InterruptedException e ) {
            log.debug "Oops.. retry sleep interrupted", e
        }

        // invoke the retry handler
        opts.retryAction?.call(t)
    }

    @PackageScope
    synchronized void slowDown0() {
        final long now = System.currentTimeMillis()
        if( now-lastFailureMillis > opts.errorBurstDelayMillis ) {
            log.debug "Throttler slowing down event -- delta=${now-lastFailureMillis}; limiter=$opts.limiter"
            lastFailureMillis = now
            backOff0(opts.limiter)
        }
        else {
            log.debug "Throttler skipped slowing down event -- delta=${now-lastFailureMillis}"
        }
    }

    @PackageScope
    void backOff0(RateLimiter limiter) {
        if( limiter == null )
            return

        final double newRate = Math.max( limiter.rate / opts.backOffFactor, opts.backOffMinRate )
        if( newRate == limiter.rate )
            return

        limiter.setRate(newRate)
        opts.rateLimitChangeAction?.call(new RateUnit(newRate))
    }

    @PackageScope
    synchronized void speedUp0(RateLimiter limiter) {
        if( limiter == null )
            return

        final double newRate = Math.min( limiter.rate * opts.rampUpFactor, opts.rampUpMaxRate )
        if( newRate == limiter.rate)
            return

        limiter.setRate(newRate)
        opts.rateLimitChangeAction?.call(new RateUnit(newRate))
    }

    /**
     * Invocation handler used by {@link ClientProxyThrottler}
     *
     * See {@link ClientProxyThrottler#invokeMethod(java.lang.String, java.lang.Object)}
     *
     *
     * @param target
     * @param name
     * @param args
     * @return
     */
    Object doInvoke0(final Object target, final String name, final Object args) {
        assert ((Object[])args).length==1

        final task = new Recoverable() {
            @Override Object invoke() {
                // invoke the closure passing the client object
                ((Closure)((Object[])args)[0]).call(target)
            }
        }

        this.submit((Callable)task)
    }

    /**
     * Invocation handler used by {@link ClientProxyThrottler}
     *
     * See {@link ClientProxyThrottler#invokeMethod(java.lang.String, java.lang.Object)}
     *
     * @param name The name of the method to execute
     * @param args The arguments of the method
     * @return The method invocation result
     */
    Object doInvoke1(final Object target, final String name, final Object args) {

        final task = new Recoverable() {
            @Override  Object invoke() {
                target.invokeMethod(name,args)
            }
        }

        def list = new ArrayList(1)
        list.add(task)
        def result = this.invokeAll(list)

        return result.get(0).get()
    }


}


