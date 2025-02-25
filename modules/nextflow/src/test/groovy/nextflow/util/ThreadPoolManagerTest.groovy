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

import spock.lang.Specification

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
}
