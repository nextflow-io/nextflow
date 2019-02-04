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

import java.util.concurrent.atomic.AtomicInteger

import groovy.util.logging.Slf4j
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ClientProxyThrottlerTest extends Specification {

    static class MyClient {
        int foo() {
            return 42 // you know way
        }

        String bar( String str ) {
            return str.reverse()
        }

        protected sayHello0() {
            return 'Hello world'
        }

        String sayHello() {
            // note internal method invocation are not forwarded to the execution pool
            return sayHello0()
        }

        def runThis( Closure c ) {
            c.call()
        }
    }

    static class MyClientProxy extends ClientProxyThrottler<MyClient> {

        @Delegate
        MyClient target

        MyClientProxy(MyClient client, ThrottlingExecutor.Options opts) {
            super(client, opts)
            target = client
        }

        MyClientProxy(MyClient client, ThrottlingExecutor executor) {
            super(client, executor, [beta: 10 as byte])
            target = client
        }
    }

    def 'should invoke target client method' () {

        given:
        def success = new AtomicInteger()
        def opts = new ThrottlingExecutor.Options()
                        .onSuccess { success.incrementAndGet() }
        def proxy = new MyClientProxy( new MyClient(), opts )

        expect:
        proxy.foo() == 42
        proxy.bar('hello') == 'olleh'
        proxy.sayHello() == 'Hello world'

        success.get()==3
    }


    def 'should throttle requests' () {

        given:
        def opts = new ThrottlingExecutor.Options().withRateLimit('1/sec')
        def client = new MyClient()
        def proxy = new MyClientProxy(client, opts)
        long begin = System.currentTimeMillis()
        sleep 50

        when:
        int count = 0
        for( int i=0; i<5; i++ ) {
            proxy.runThis( { log.info "tick=${count++}" } )
        }
 
        long delta = System.currentTimeMillis()-begin
        then:
        delta > 3_000
        delta < 5_000
        count == 5
    }

    def 'should invoke onSuccess callback' () {
        given:
        def success = new AtomicInteger()
        def client = new MyClient()
        def opts = new ThrottlingExecutor.Options().onSuccess { success.incrementAndGet() }
        def proxy = new MyClientProxy(client, opts)

        when:
        proxy.runThis { println "task 1" }
        proxy.runThis { println "task 2" }
        then:
        success.get() ==2
    }


    def 'should invoke onFailure callback' () {
        given:
        def success = new AtomicInteger()
        def failure = new AtomicInteger()
        def illegal = new AtomicInteger()
        def retry = new AtomicInteger()

        def opts = new ThrottlingExecutor.Options()
                    .onSuccess { success.incrementAndGet() }
                    .onRetry { retry.incrementAndGet() }
                    .onFailure { t -> failure.incrementAndGet(); if( t instanceof IllegalArgumentException ) illegal.incrementAndGet()}

        def proxy = new MyClientProxy(new MyClient(), opts)

        when:
        proxy.runThis { println "task 1" }
        proxy.runThis { throw new IllegalArgumentException('error 1') }
        proxy.runThis { throw new UnsupportedOperationException('error 2') }
        proxy.runThis( { println "task 2" } )

        then:
        illegal.get() == 1
        failure.get() == 2
        success.get() == 2
        retry.get() == 0
    }

    def 'should retry execution' () {

        given:
        def success = new AtomicInteger()
        def failure = new AtomicInteger()
        def limitChange = new AtomicInteger()
        def execCount = new AtomicInteger()
        def retry = new AtomicInteger()

        def opts = new ThrottlingExecutor.Options()
                        .withRateLimit('1000/sec')
                        .withAutoThrottle()
                        .withErrorBurstDelay(Duration.of('5 sec'))
                        .onSuccess { success.incrementAndGet() }
                        .onFailure { failure.incrementAndGet() }
                        .onRateLimitChange { limitChange.incrementAndGet() }
                        .onRetry { retry.incrementAndGet() }
                        .retryOn(IllegalArgumentException)

        def proxy = new MyClientProxy(new MyClient(), opts)

        when:
        proxy.runThis { execCount.incrementAndGet(); if(execCount.get()<=2) throw new IllegalArgumentException() }

        then:
        success.get() == 1
        failure.get() == 0
        retry.get() == 2
        execCount.get() == 3
        limitChange.get() == 1
    }

    def 'should retry execution with max retries' () {

        given:
        def success = new AtomicInteger()
        def failure = new AtomicInteger()
        def limitChange = new AtomicInteger()
        def execCount = new AtomicInteger()

        def opts = new ThrottlingExecutor.Options()
                        .withRateLimit('1000/sec')
                        .withAutoThrottle()
                        .withMaxRetries(1)
                        .onSuccess { success.incrementAndGet() }
                        .onFailure { failure.incrementAndGet() }
                        .onRateLimitChange { limitChange.incrementAndGet() }
                        .retryOn(IllegalArgumentException)

        def proxy = new MyClientProxy(new MyClient(), opts)

        when:
        proxy.runThis { execCount.incrementAndGet(); if(execCount.get()<=2) throw new IllegalArgumentException() }

        then:
        success.get() == 0
        failure.get() == 1
        execCount.get() == 2
        limitChange.get() == 1
    }

    def 'should block execution' () {

        given:
        def opts = new ThrottlingExecutor.Options()
                                .withPoolSize(1)
                                .withQueueSize(1)

        def proxy = new MyClientProxy(new MyClient(), opts)

        when:
        def count=new AtomicInteger()
        def begin = System.currentTimeMillis()

        for( int i=0; i<10; i++ ) {
            proxy.runThis { count.incrementAndGet(); sleep 500 }
        }

        def delta = System.currentTimeMillis()-begin
        println "delta=$delta"
        then:
        count.get() == 10
        delta>5_000
    }


    def 'should call invoke1' () {
        given:
        def client = new MyClient()
        def exec = Mock(ThrottlingExecutor)
        def proxy = new MyClientProxy(client, exec)
        byte ZERO = 0 as byte
        byte _10 = 10 as byte
        when:
        def result = proxy.alpha('hello world')
        then:
        1 * exec.doInvoke1(client, 'alpha', ['hello world'] as Object[], ZERO)  >> 'OK'
        result == 'OK'

        // it should use a different priority for method `beta`, the priority is defined in the
        // MyClientProxy class definition
        when:
        result = proxy.beta('hola mundo')
        then:
        1 * exec.doInvoke1(client, 'beta', ['hola mundo'] as Object[], _10)  >> 'OK'
        result == 'OK'
    }

}
