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

import java.util.concurrent.atomic.AtomicInteger

import com.google.common.util.concurrent.RateLimiter
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

    static class MyClientProxy extends ClientProxyThrottler {

        @Delegate
        MyClient target

        MyClientProxy(MyClient client, ThrottlingExecutor.Options opts) {
            super(client, opts)
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

    def 'should invoke async' () {

        given:
        def success = new AtomicInteger()
        def opts = new ThrottlingExecutor.Options() .onSuccess { success.incrementAndGet() }
        def proxy = new MyClientProxy( new MyClient(), opts )
        
        when:
        def f = proxy.async { MyClient c -> c.foo() }
        while( !f.isDone() ) sleep 50
        then:
        f.get() == 42
        success.get()==1

    }


    def 'should set properties' () {

        given:
        def success = { 1 }
        def failure = { 2 }
        def change = { 3 }
        def retry = { 4 }
        def retryOn = { false }

        when:
        def builder = new ThrottlingExecutor.Options()
                .withRateLimit('20/sec')
                .withPoolSize(10)
                .withMaxPoolSize(20)
                .withQueueSize(50)
                .withKeepAlive(Duration.of('3 min'))
                .withMaxRetries(11)
                .withRampUp(50.0, 1.5, 500)
                .withBackOff('10/sec', 1.9)
                .withAutoThrottle(true)
                .withErrorBurstDelay(Duration.of('5 sec'))
                .onFailure(failure)
                .onSuccess(success)
                .onRetry(retry)
                .onRateLimitChange(change)
                .retryOn(retryOn)

        then:
        builder.limiter.rate == RateLimiter.create(20).rate
        builder.poolSize == 10
        builder.maxPoolSize == 20
        builder.queueSize == 50
        builder.maxRetries == 11
        builder.keepAlive == Duration.of('3 min')
        builder.failureAction.is(failure)
        builder.successAction.is(success)
        builder.retryAction.is(retry)
        builder.rateLimitChangeAction.is(change)
        builder.retryCondition.is(retryOn)
        builder.errorBurstDelayMillis == 5_000
        builder.rampUpInterval == 500
        builder.rampUpFactor == 1.5
        builder.rampUpMaxRate == 50.0
        builder.backOffMinRate == 10d
        builder.backOffFactor == 1.9f

    }

    def 'should configure with options map' () {

        given:
        def OPTS = [
                poolSize: 1,
                maxPoolSize: 2,
                queueSize: 3,
                maxRetries: 4,
                keepAlive: '5 sec',
                errorBurstDelay: '6 sec',
                rampUpInterval: 7,
                rampUpFactor: 8,
                rampUpMaxRate: 11.1,
                backOffMinRate: '9/sec',
                backOffFactor: 10,
                rateLimit: '20/s',
                autoThrottle: true,

        ]

        when:
        def builder = new ThrottlingExecutor.Options().withOptions(OPTS)

        then:
        builder.poolSize == 1
        builder.maxPoolSize == 2
        builder.queueSize == 3
        builder.maxRetries == 4
        builder.keepAlive == Duration.of('5 sec')
        builder.autoThrottle == true
        builder.errorBurstDelayMillis == 6_000
        builder.rampUpInterval == 7
        builder.rampUpFactor == 8
        builder.rampUpMaxRate == 11.1
        builder.backOffMinRate == 9
        builder.backOffFactor == 10
        builder.limiter.rate == RateLimiter.create(20).rate
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
            proxy.runThis( { log.info "thick=${count++}" } )
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

}
