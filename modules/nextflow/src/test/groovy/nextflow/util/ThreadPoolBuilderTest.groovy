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

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit

import spock.lang.Specification

class ThreadPoolBuilderTest extends Specification {

    def 'should create thread pool with default' () {

        when:
        def builder = new ThreadPoolBuilder()
                .withMinSize(1)
                .withMaxSize(10)
        then:
        builder.getMinSize() == 1
        builder.getMaxSize() == 10

        when:
        def pool = builder.build()
        then:
        pool.getCorePoolSize() == 1
        pool.getMaximumPoolSize() == 10
        pool.getKeepAliveTime(TimeUnit.MILLISECONDS) == 60_000
        pool.getThreadFactory() instanceof CustomThreadFactory
        pool.getRejectedExecutionHandler() instanceof ThreadPoolExecutor.CallerRunsPolicy
        and:
        builder.getName().startsWith('nf-thread-pool-')
        builder.getWorkQueue() instanceof LinkedBlockingQueue
    }

    def 'should create pool with all settings' () {
        when:
        def builder = new ThreadPoolBuilder()
                .withName('foo')
                .withMinSize(1)
                .withMaxSize(10)
                .withKeepAliveTime(100)
                .withQueueSize(1000)
                .withAllowCoreThreadTimeout(true)
                .withRejectionPolicy(new ThreadPoolExecutor.AbortPolicy())
        then:
        builder.name == 'foo'
        builder.getMinSize() == 1
        builder.getMaxSize() == 10
        builder.keepAliveTime == 100
        builder.queueSize == 1000
        builder.allowCoreThreadTimeout
        builder.rejectionPolicy instanceof ThreadPoolExecutor.AbortPolicy

        when:
        def pool = builder.build()
        then:
        pool.getCorePoolSize() == 1
        pool.getMaximumPoolSize() == 10
        pool.getKeepAliveTime(TimeUnit.MILLISECONDS) == 100
        pool.getThreadFactory() instanceof CustomThreadFactory
        pool.getRejectedExecutionHandler() instanceof ThreadPoolExecutor.AbortPolicy
    }
}
