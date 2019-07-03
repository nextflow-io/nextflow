/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.util.concurrent.ExecutorService
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
/**
 * A thread pool executor that block the submitting thread when
 * the pool queue is full.
 *
 * Inspired by
 *   http://fahdshariff.blogspot.com/2013/11/throttling-task-submission-with.html
 *   https://www.logicbig.com/tutorials/core-java-tutorial/java-multi-threading/thread-pools.html
 *
 * See also
 *   https://community.oracle.com/docs/DOC-983726
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BlockingThreadExecutorFactory {

    private static String DEFAULT_NAME = 'BlockingThreadExecutor'
    private static int DEFAULT_POOL_SIZE = Runtime.runtime.availableProcessors()
    private static Duration DEFAULT_KEEP_ALIVE = Duration.of('60sec')
    private static Duration DEFAULT_MAX_AWAIT = Duration.of('36 hour')

    String name
    Integer maxThreads
    Integer maxQueueSize
    Duration keepAlive
    Duration maxAwait

    static ExecutorService create(String name) {
        new BlockingThreadExecutorFactory().withName(name).create()
    }

    BlockingThreadExecutorFactory withName(String name) {
        this.name = name
        return this
    }

    BlockingThreadExecutorFactory withMaxThreads(int max) {
        this.maxThreads = max
        return this
    }

    BlockingThreadExecutorFactory withMaxQueueSize(int max) {
        this.maxQueueSize = max
        return this
    }

    BlockingThreadExecutorFactory withKeepAlive(Duration duration) {
        this.keepAlive = duration
        return this
    }

    ExecutorService create() {
        final session = Global.session as Session
        if( session && name ) {
            maxQueueSize = session.config.navigate("threadPool.${name}.maxQueueSize") as Integer
            keepAlive = session.config.navigate("threadPool.${name}.keepAlive") as Duration
            maxThreads = session.config.navigate("threadPool.${name}.maxThreads") as Integer
            maxAwait = session.config.navigate("threadPool.${name}.maxAwait") as Duration
            if( !maxThreads )
                maxThreads = (session.config.poolSize as Integer ?: 0) *3
        }

        if( !maxThreads )
            maxThreads = DEFAULT_POOL_SIZE
        if( !maxQueueSize )
            maxQueueSize = maxThreads *3
        if( keepAlive == null )
            keepAlive = DEFAULT_KEEP_ALIVE
        if( maxAwait == null )
            maxAwait = DEFAULT_MAX_AWAIT

        log.debug "Thread pool name=$name; maxThreads=$maxThreads; maxQueueSize=$maxQueueSize; keepAlive=$keepAlive"
        create0()
    }

    private ExecutorService create0() {
        assert maxThreads>0, "Thread pool size must be greater than zero"
        assert maxQueueSize>=maxThreads, "Thread queue size must be greater or equal to the pool size"

        if( keepAlive == null )
            keepAlive = DEFAULT_KEEP_ALIVE

        final prefix = (name ?: DEFAULT_NAME) + '-thread'
        final pool = new ThreadPoolExecutor(
                    maxThreads,
                    Integer.MAX_VALUE,
                    keepAlive.getMillis(),
                    TimeUnit.MILLISECONDS,
                    new BlockingBlockingQueue<Runnable>(maxQueueSize-maxThreads),
                    new CustomThreadFactory(prefix),
                    new ThreadPoolExecutor.CallerRunsPolicy() )
                    // ^^^^^
                    // note: the CallerRunsPolicy is only used in the case the queue offer method
                    // (which blocking in the BlockingBlockingQueue implementation) return `false`
                    // because it gets an interrupted execution

        pool.allowCoreThreadTimeOut(true)
        return pool
    }
}
