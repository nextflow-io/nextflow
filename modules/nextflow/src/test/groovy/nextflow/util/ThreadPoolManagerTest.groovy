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

package nextflow.util

import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

import spock.lang.Specification
import spock.util.concurrent.PollingConditions

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ThreadPoolManagerTest extends Specification {

    def 'should create pool with defaults' () {
        when:
        def pool = new ThreadPoolManager('foo').create()
        then:
        pool.getCorePoolSize()  == ThreadPoolManager.DEFAULT_MIN_THREAD
        pool.getMaximumPoolSize() == ThreadPoolManager.DEFAULT_MAX_THREAD
    }

    def 'should create pool with with max threads' () {
        when:
        def pool = new ThreadPoolManager('foo').withMaxThreads(100).create()
        then:
        pool.getCorePoolSize()  == ThreadPoolManager.DEFAULT_MIN_THREAD
        pool.getMaximumPoolSize() == 100
    }

    def 'should create pool with with min threads' () {
        when:
        def pool = new ThreadPoolManager('foo').withMaxThreads(1).create()
        then:
        pool.getCorePoolSize()  == 1
        pool.getMaximumPoolSize() == 1
    }

    def 'should create pool with with config' () {
        given:
        def CONFIG = [threadPool:[foo:[minThreads: 2, maxThreads: 8]]]
        when:
        def pool = new ThreadPoolManager('foo').withConfig(CONFIG).create()
        then:
        pool.getCorePoolSize()  == 2
        pool.getMaximumPoolSize() == 8
    }

    def 'should create pool with with max threads and config' () {
        given:
        def CONFIG = [threadPool:[foo:[minThreads: 2, maxThreads: 8]]]
        when:
        def pool = new ThreadPoolManager('foo').withMaxThreads(5).withConfig(CONFIG).create()
        then:
        pool.getCorePoolSize()  == 2
        pool.getMaximumPoolSize() == 8
    }

    def 'should count running tasks in virtual thread pool' () {
        given:
        def manager = new ThreadPoolManager('foo')
        def pool = manager.virtualThreadService()
        manager.@executorService = pool
        def latch = new CountDownLatch(1)

        when:
        2.times { pool.submit { latch.await() } }
        then:
        new PollingConditions(timeout: 5).eventually { assert manager.pending() == 2 }

        when:
        latch.countDown()
        pool.shutdown()
        pool.awaitTermination(5, TimeUnit.SECONDS)
        then:
        new PollingConditions(timeout: 5).eventually { assert manager.pending() == 0 }
    }
}
