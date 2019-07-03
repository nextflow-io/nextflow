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

import java.util.concurrent.TimeUnit

import groovy.util.logging.Slf4j
import nextflow.Session
import spock.lang.Specification
/**
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
class BlockingTreadExecutorFactoryTest extends Specification {

    def 'should create thread pool with default config' () {
        given:
        def CPUS = Runtime.runtime.availableProcessors()
        def pool = new BlockingThreadExecutorFactory()

        when:
        pool.create()
        then:
        pool.maxThreads == CPUS
        pool.maxQueueSize == CPUS * 3
        pool.keepAlive == Duration.of('60 sec')
        pool.maxAwait == Duration.of('36 hours')

    }

    def 'should set properties' () {
        when:
        def pool = new BlockingThreadExecutorFactory()
                .withName('foo')
                .withMaxThreads(10)
                .withMaxQueueSize(20)
                .withKeepAlive(Duration.of(30))

        then:
        pool.name == 'foo'
        pool.maxThreads == 10
        pool.maxQueueSize == 20
        pool.keepAlive.millis == 30
    }

    def 'should create thread pool with session config' () {
        given:
        new Session([threadPool: [foo:[
                maxThreads: 5,
                maxQueueSize: 10,
                keepAlive: '15 sec',
                maxAwait: '5 min'
        ]]])
        def pool = new BlockingThreadExecutorFactory().withName('foo')

        when:
        pool.create()
        then:
        pool.maxThreads == 5
        pool.maxQueueSize == 10
        pool.keepAlive == Duration.of('15 sec')
        pool.maxAwait == Duration.of('5 min')
    }

    def 'should execute tasks' () {

        given:
        def factory = new BlockingThreadExecutorFactory()
                    .withMaxThreads(1)
                    .withMaxQueueSize(2)

        when:
        def executor = factory.create()
        def start = System.currentTimeMillis()
        executor.execute( { log.info 'hello 1'; sleep 1_000 } )
        executor.execute( { log.info 'hello 2'; sleep 1_000 } )
        executor.execute( { log.info 'hello 3'; sleep 1_000 } )
        def delta1 = System.currentTimeMillis() - start
        
        executor.shutdown()
        executor.awaitTermination(30, TimeUnit.SECONDS)
        def delta2 = System.currentTimeMillis() - start

        then:
        // since queueSize == 2, up to two tasks can be queued in the thread pool
        // therefore delta1 must be greater than 1 second
        delta1>= 1_000 && delta1 < 2_000
        // since poolSize
        delta2>= 3_000 && delta2 < 4_000

    }

    def 'should create pool' () {
        given:
        def factory = new BlockingThreadExecutorFactory(maxThreads: 10)
        def executor = factory.create()

        when:
        executor.submit( { println 'hello' })
        executor.shutdown()
        executor.awaitTermination(1, TimeUnit.SECONDS)
        then:
        noExceptionThrown()
    }
}
