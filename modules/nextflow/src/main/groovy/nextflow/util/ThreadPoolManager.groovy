/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.util

import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.ISession
import nextflow.Session
/**
 * Holder object for file transfer thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ThreadPoolManager {

    private static final AtomicInteger poolCount = new AtomicInteger()
    
    final static public int DEFAULT_MIN_THREAD = 10
    final static public int DEFAULT_MAX_THREAD = Math.max(DEFAULT_MIN_THREAD, Runtime.runtime.availableProcessors()*3)
    final static public int DEFAULT_QUEUE_SIZE = 10_000
    final static public Duration DEFAULT_KEEP_ALIVE =  Duration.of('60sec')
    final static public Duration DEFAULT_MAX_AWAIT = Duration.of('12 hour')

    private Integer minThreads = DEFAULT_MIN_THREAD
    private Integer maxThreads = DEFAULT_MAX_THREAD
    private Integer maxQueueSize = DEFAULT_QUEUE_SIZE
    private Duration keepAlive = DEFAULT_KEEP_ALIVE
    private Boolean allowThreadTimeout
    private Duration maxAwait = DEFAULT_MAX_AWAIT
    private ExecutorService executorService
    private String name
    private String waitMsg = "Waiting for file transfers to complete (%d files)"
    private String exitMsg = "Exiting before file transfers were completed -- Some files may be lost"

    ThreadPoolManager(String name) {
        this.name = name
    }

    ThreadPoolManager withMaxThreads(int maxThreads) {
        if( maxThreads>0 )
            this.maxThreads = maxThreads
        return this
    }

    ThreadPoolManager withConfig(Map config) {
        this.minThreads = config.navigate("threadPool.${name}.minThreads", minThreads) as Integer
        this.maxThreads = config.navigate("threadPool.${name}.maxThreads", maxThreads) as Integer
        this.maxQueueSize = config.navigate("threadPool.${name}.maxQueueSize", maxQueueSize) as Integer
        this.keepAlive = config.navigate("threadPool.${name}.keepAlive", keepAlive) as Duration
        this.allowThreadTimeout = config.navigate("threadPool.${name}.allowThreadTimeout", false) as Boolean
        this.maxAwait = config.navigate("threadPool.${name}.maxAwait", maxAwait) as Duration
        return this
    }

    ThreadPoolManager withShutdownMessage(String waitMsg, String exitMsg) {
        this.waitMsg = waitMsg
        this.exitMsg = exitMsg
        return this
    }

    ExecutorService create() {
        if( minThreads>maxThreads ) {
            log.debug("Thread pool '$name' minThreads ($minThreads) cannot be greater than maxThreads ($maxThreads) - Setting minThreads to $maxThreads")
            minThreads = maxThreads
        }

        return executorService = Threads.useVirtual()
                ? Executors.newThreadPerTaskExecutor(new CustomThreadFactory(name ?: "nf-thread-pool-${poolCount.getAndIncrement()}".toString()))
                : legacyThreadPool()
    }

    protected ExecutorService legacyThreadPool() {
        new ThreadPoolBuilder()
                .withName(name)
                .withMinSize(minThreads)
                .withMaxSize(maxThreads)
                .withQueueSize(maxQueueSize)
                .withKeepAliveTime(keepAlive)
                .withAllowCoreThreadTimeout(allowThreadTimeout)
                .build()
    }

    ExecutorService createAndRegisterShutdownCallback(Session session) {
        final result = create()
        // register the cleanup callback
        Global.onCleanup( (it) -> shutdown(session))
        return result
    }

    void shutdown(ISession session) {
        final sess = (Session) session
        shutdown( sess != null && sess.aborted )
    }

    void shutdown(boolean hard) {
        if( !executorService )
            return

        if( hard ) {
            executorService.shutdownNow()
            return
        }

        executorService.shutdown()
        // wait for remaining threads to complete
        ThreadPoolHelper.await(executorService, maxAwait, waitMsg, exitMsg)
        log.debug "Thread pool '$name' shutdown completed (hard=$hard)"
    }

    static ExecutorService create(String name, int maxThreads=0) {
        final session = Global.session as Session
        new ThreadPoolManager(name)
            .withMaxThreads(maxThreads) // default max threads
            .withConfig(session.config) // config can override maxThread specified above
            .createAndRegisterShutdownCallback(session)
    }
}
