/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import java.util.concurrent.Callable
import java.util.concurrent.FutureTask
import java.util.concurrent.PriorityBlockingQueue
import java.util.concurrent.RejectedExecutionException
import java.util.concurrent.RunnableFuture
import java.util.concurrent.Semaphore
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
/**
 * Implements the throttling and retry logic
 *
 * Inspired by
 *   http://fahdshariff.blogspot.com/2013/11/throttling-task-submission-with.html
 *   https://www.logicbig.com/tutorials/core-java-tutorial/java-multi-threading/thread-pools.html
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ThrottlingExecutor extends ThreadPoolExecutor {


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

        /** the thread pool name */
        String poolName

        /**
         * The rate limiter object submitting request at the specified rate
         */
        RateLimiter limiter

        /**
         * Closure invoked when the task is successfully completed
         */
        Closure successAction

        /**
         * Closure invoked when the task execution failed (after all retry attempts failed)
         */
        Closure failureAction

        /**
         * Closure invoked on each retry attempt
         */
        Closure retryAction

        /**
         * Closure invoked to signal that the rate limit has been automatically adjusted
         */
        Closure rateLimitChangeAction

        /**
         * A closure receiving a {@link Throwable} object as parameter a returning {@code true}
         * if one or more retries need to be tempted for a failing task.
         */
        Closure<Boolean> retryCondition

        /**
         * The thread pool size
         */
        int poolSize

        /**
         * The thread pool max size
         */
        int maxPoolSize

        /**
         * The size of the queue where tasks are kept when the thread pool is full
         */
        int queueSize

        /**
         * Number of attempts a task execution is retried when an exception
         * matching the #retryCondition is raised
         */
        int maxRetries = Integer.MAX_VALUE

        /**
         * Timeout period after when idle thread will be evicted
         */
        Duration keepAlive = Duration.of('2 min')

        /**
         * When {@code true} the executor automatically decreases and increases
         * the execution rate to adjust to failure events
         */
        boolean autoThrottle

        /**
         * Period during which error conditions are ignored
         */
        Duration errorBurstDelay = Duration.of('1 sec')

        /**
         * Number of successful task execution after which the execution rate is increased
         */
        int rampUpInterval = 100

        /**
         * Factor by which the rate limit is increased. It must be greater than 1.
         */
        float rampUpFactor = 1.2

        /**
         * Max allowed execution rate
         */
        double rampUpMaxRate = Double.MAX_VALUE

        /**
         * Factor by which the execution rate is decreased on a failed execution
         */
        float backOffFactor = 2.0f

        /**
         * Min allowed execution rate
         */
        double backOffMinRate = RateUnit.of('1 / min').rate

        /**
         * Period after which the task execution should be retried
         */
        Duration retryDelay = Duration.of('1 sec')

        /**
         * Initialise missing attribute to default values
         *
         * @return The {@link Options} object itself
         */
        Options normalise() {
            if( !poolSize )
                poolSize = Runtime.runtime.availableProcessors()
            if( !maxPoolSize || poolSize>maxPoolSize )
                maxPoolSize = poolSize
            if( keepAlive==null )
                keepAlive = Duration.of(0)
            if( !poolName )
                poolName = ThrottlingExecutor.simpleName

            return this
        }

        void setRampUpFactor(float factor) {
            if( factor > 1 )
                rampUpFactor = factor
            else
                log.warn "Throttling executor ramp-up factor must be greater than 1"
        }

        void setBackOffFactor(float factor) {
            if( factor > 1 )
                backOffFactor = factor
            else
                log.warn "Throttling executor back-off factor must be greater than 1"
        }

        Options withPoolName(String n) {
            this.poolName = n
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

        Options withQueueSize(int size) {
            queueSize = size
            return this
        }

        Options withRateLimit(String rate) {
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

        Options retryOn( @ClosureParams(value = SimpleType.class, options = "java.lang.Throwable") Closure<Boolean> condition) {
            retryCondition = condition
            return this
        }

        Options withRampUp(double maxRate, Float factor=null, Integer interval=null) {
            this.rampUpMaxRate = maxRate
            if( factor!=null ) {
                this.setRampUpFactor(factor)
            }
            if( interval!=null )
                this.rampUpInterval = interval
            return this
        }

        Options withBackOff(double minRate, Float factor=null) {
            backOffMinRate = minRate
            if( factor != null )
                setBackOffFactor(factor)
            return this
        }

        Options withBackOff(String minRate, float factor) {
            withBackOff(new RateUnit(minRate).rate, factor)
        }

        Options withErrorBurstDelay(Duration duration) {
            this.errorBurstDelay = duration
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
                errorBurstDelay = opts.errorBurstDelay as Duration

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

    @CompileStatic
    private static class PriorityTask<V> extends FutureTask<V> implements Comparable<PriorityTask> {

        private short priority

        PriorityTask(Callable<V> callable, short p) {
            super(callable)
            this.priority = p
        }

        PriorityTask(Runnable runnable, V result, short p) {
            super(runnable,result)
            this.priority = p
        }

        @Override
        int compareTo(PriorityTask that) {
            // important: the comparator is inverted because it implements
            // an *inverted* natural ordering i.e. higher values comes first
            // therefore are taken first from the head of the queue
            return that.priority <=> this.priority
        }
    }

    abstract static class Recoverable implements Callable, Runnable {

        private ThrottlingExecutor executor

        private Options opts

        private int retryCount

        abstract protected Object invoke()

        protected Recoverable setOwner(ThrottlingExecutor e) {
            this.executor=e
            this.opts=e.opts
            return this
        }

        protected byte getPriority() { 0 }

        @Override
        final void run() { call() }

        @Override
        final Object call() throws Exception {
            while(true)
                try {
                    final result = invoke()
                    handleSuccess0(result)
                    return result
                }
                catch (Throwable t) {
                    try {
                        if( handleRetry0(t) )
                            continue
                        // invoke error handling callback
                        onFailure(t)
                        return null
                    }
                    catch (Throwable e) {
                        log.error("Unexpected failure", e)
                        return null
                    }
                }

        }

        /**
         * Subclasses can provide task specific error handling policy
         * @param the error thrown
         */
        protected void onFailure(Throwable t) {
            if( opts.failureAction ) {
                opts.failureAction.call(t)
            }
            else {
                log.error("Unable to process client request", t)
            }
        }

        private void handleSuccess0(Object result) {
            // invoke the success handler
            opts.successAction?.call(result)
            // try to speed-up the submission rate
            if( executor.canSpeedUp0() ) {
                executor.speedUp0(opts.limiter)
            }
        }

        private boolean handleRetry0(Throwable t) {
            log.trace "Failure exception=$t"
            // reset the success counter
            synchronized (executor.countLock) { executor.successCount=0 }
            // check if can retry the task execution
            if( opts.retryCondition && opts.retryCondition.call(t) && retryCount++ < opts.maxRetries ) {
                // slow down the submission throughput
                if( opts.autoThrottle )
                    executor.slowDown0()
                // wait before retry this task execution
                sleep(opts.retryDelay.millis * retryCount)
                // notify the retry event
                opts.retryAction?.call(t)
                return true
            }

            return false
        }
    }

    // ---- executor settings ----

    private Options opts

    private Semaphore semaphore

    private final Object countLock = new Object()

    private int successCount

    private long lastFailureMillis

    private volatile double timeToAcquire

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
                new PriorityBlockingQueue<Runnable>(),
                new CustomThreadFactory(opts.poolName, new ExceptionHandler0()))

        // the executor options
        this.opts = opts
        // allow core thread to be removed from cache
        allowCoreThreadTimeOut(true)
        // the semaphore is bounding both the number of tasks currently executing
        // and those queued up
        if( opts.queueSize )
            this.semaphore = new Semaphore(opts.maxPoolSize + opts.queueSize)
    }

    @Override
    void execute(final Runnable task) {
        Runnable wrapper = task instanceof PriorityTask ? (Runnable)task : wrap(task)
        semaphore?.acquire()

        try {
            super.execute(wrapper)
        }
        catch (RejectedExecutionException e) {
            semaphore?.release()
            throw e
        }
    }

    protected <T> Recoverable wrap(Callable<T> task) {
        Recoverable result
        if( task instanceof Recoverable )
            result = task

        else
            result = new Recoverable() {
                @Override Object invoke() { return task.call() }
            }

        result.setOwner(this)
    }

    protected Recoverable wrap(Runnable task) {
        Recoverable result
        if( task instanceof Recoverable )
            result = task

        else
            result = new Recoverable() {
                @Override Object invoke() { return task.run() }
            }

        result.setOwner(this)
    }

    @Override
    protected <T> RunnableFuture<T> newTaskFor(Callable<T> task) {
        final recoverable = wrap(task)
        return new PriorityTask<T>(recoverable, recoverable.priority);
    }

    @Override
    protected <T> RunnableFuture<T> newTaskFor(Runnable task, T value) {
        final recoverable = wrap(task)
        return new PriorityTask<T>(recoverable, value, recoverable.priority);
    }

    protected void beforeExecute(Thread thread, Runnable task) {
        super.beforeExecute(thread,task)
        log.trace "Before acquire limit=$opts.limiter"
        if( opts.limiter ) {
            timeToAcquire = opts.limiter.acquire()
            if( timeToAcquire )
            log.trace "Time spent to enforce rate=$timeToAcquire"
        }
    }

    @Override
    protected void afterExecute(Runnable task, Throwable ex) {
        log.trace "After execute task -- ${this}"
        super.afterExecute(task,ex)
        semaphore?.release()
    }

    @PackageScope
    synchronized void slowDown0() {
        final long now = System.currentTimeMillis()
        if( now-lastFailureMillis > opts.errorBurstDelay.millis ) {
            log.trace "Throttler slowing down event -- delta=${now-lastFailureMillis}; limiter=$opts.limiter"
            lastFailureMillis = now
            backOff0(opts.limiter)
        }
        else {
            log.trace "Throttler skipped slowing down event -- delta=${now-lastFailureMillis}"
        }
    }

    @PackageScope
    void backOff0(RateLimiter limiter) {
        if( limiter == null || !opts.backOffFactor )
            return

        final double norm = limiter.rate / opts.backOffFactor
        final double newRate = Math.max( norm, opts.backOffMinRate )
        if( newRate == limiter.rate )
            return

        limiter.setRate(newRate)
        opts.rateLimitChangeAction?.call(new RateUnit(newRate))
    }

    private boolean canSpeedUp0() {
        if( !opts.autoThrottle )
            return false
        // increase execution rate only when some time is spent to enforce the current limit
        // otherwise is useless because there no execution pressure 
        if( !timeToAcquire )
            return false
        
        synchronized (countLock) {
            def result = opts.rampUpInterval>0 && ++successCount >= opts.rampUpInterval
            if( result )
                successCount=0
            return result
        }
    }

    @PackageScope
    synchronized void speedUp0(RateLimiter limiter) {
        if( limiter == null || !opts.rampUpFactor )
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
    @Deprecated
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
    Object doInvoke1(final Object target, final String name, final Object args, final byte priority) {

        final task = new Recoverable() {
            @Override Object invoke() {
                log.trace "Invoke client method [priority=$priority] >> $name($args)"
                target.invokeMethod(name,args)
            }

            @Override byte getPriority() { priority }
        }

        final result = submit((Callable)task)
        // get the result waiting if necessary
        return result.get()
    }
}


