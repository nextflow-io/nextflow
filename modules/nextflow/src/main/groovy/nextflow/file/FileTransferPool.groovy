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

package nextflow.file

import java.util.concurrent.ExecutorService
import java.util.concurrent.ThreadPoolExecutor

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.util.Duration
import nextflow.util.ThreadPoolBuilder
import nextflow.util.ThreadPoolHelper

/**
 * Holder object for file transfer thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FileTransferPool {

    final static private DEFAULT_MIN_THREAD = 1
    final static private DEFAULT_MAX_THREAD = Math.min(Runtime.runtime.availableProcessors()*2, 10)
    final static private DEFAULT_QUEUE = 10_000
    final static private DEFAULT_KEEP_ALIVE =  Duration.of('60sec')
    final DEFAULT_MAX_AWAIT = Duration.of('12 hour')

    private Integer minThreads
    final private Integer maxThreads
    final private Integer maxQueueSize
    final private Duration keepAlive
    final private Boolean allowThreadTimeout
    final private Duration maxAwait
    final private ThreadPoolExecutor executorService

    static private FileTransferPool instance0

    FileTransferPool(Map config) {
        this.minThreads = config.navigate("threadPool.FileTransfer.minThreads", DEFAULT_MIN_THREAD) as Integer
        this.maxThreads = config.navigate("threadPool.FileTransfer.maxThreads", DEFAULT_MAX_THREAD) as Integer
        this.maxQueueSize = config.navigate("threadPool.FileTransfer.maxQueueSize", DEFAULT_QUEUE) as Integer
        this.keepAlive = config.navigate("threadPool.FileTransfer.keepAlive", DEFAULT_KEEP_ALIVE) as Duration
        this.allowThreadTimeout = config.navigate("threadPool.FileTransfer.allowThreadTimeout", false) as Boolean
        this.maxAwait = config.navigate("threadPool.FileTransfer.maxAwait", DEFAULT_MAX_AWAIT) as Duration

        if( minThreads>maxThreads ) {
            log.debug("FileTransfer minThreads ($minThreads) cannot be greater than maxThreads ($maxThreads) - Setting minThreads to $maxThreads")
            minThreads = maxThreads
        }

        executorService = new ThreadPoolBuilder()
                .withName('FileTransfer')
                .withMinSize(minThreads)
                .withMaxSize(maxThreads)
                .withQueueSize(maxQueueSize)
                .withKeepAliveTime(keepAlive)
                .withAllowCoreThreadTimeout(allowThreadTimeout)
                .build()
    }

    private ExecutorService getExecutorService0() {
        return executorService
    }

    private void shutdown0(boolean hard) {
        if( hard ) {
            executorService.shutdownNow()
            return
        }

        executorService.shutdown()
        // wait for ongoing file transfer to complete
        final waitMsg = "Waiting files transfer to complete (%d files)"
        final exitMsg = "Exiting before FileTransfer thread pool complete -- Some files maybe lost"
        ThreadPoolHelper.await(executorService, maxAwait, waitMsg, exitMsg)
    }

    @Memoized
    static synchronized ExecutorService getExecutorService() {
        final session = Global.session
        if( session == null )
            throw new IllegalStateException("Nextflow session object has not been created yet")
        instance0 = new FileTransferPool(session.getConfig())
        return instance0.getExecutorService0()
    }

    static shutdown(boolean hard) {
        if( instance0 )
            instance0.shutdown0(hard)
    }

}
