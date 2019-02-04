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

import java.util.concurrent.Callable
import java.util.concurrent.Future
import java.util.concurrent.PriorityBlockingQueue

import com.google.common.util.concurrent.RateLimiter
import groovy.util.logging.Slf4j
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ThrottlingExecutorTest extends Specification {


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
                .withRetryDelay(Duration.of('7 sec'))
                .onFailure(failure)
                .onSuccess(success)
                .onRetry(retry)
                .onRateLimitChange(change)
                .retryOn(retryOn)
                .withPoolName('foo')

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
        builder.errorBurstDelay == Duration.of('5 sec')
        builder.retryDelay == Duration.of('7 sec')
        builder.rampUpInterval == 500
        builder.rampUpFactor == 1.5
        builder.rampUpMaxRate == 50.0
        builder.backOffMinRate == 10d
        builder.backOffFactor == 1.9f
        builder.poolName == 'foo'

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
        builder.errorBurstDelay == Duration.of('6 sec')
        builder.rampUpInterval == 7
        builder.rampUpFactor == 8
        builder.rampUpMaxRate == 11.1
        builder.backOffMinRate == 9
        builder.backOffFactor == 10
        builder.limiter.rate == RateLimiter.create(20).rate
    }



    def 'should increase rate limit' () {

        given:
        def change=0
        def opts = new ThrottlingExecutor.Options() .onRateLimitChange { change++ } .normalise()
        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(1.0)

        when:
        exec.speedUp0(null)
        then:
        noExceptionThrown()
        limiter.rate == 1.0d

        when:
        exec.speedUp0(limiter)
        then:
        limiter.rate as float == 1.2d as float
        change == 1

    }

    def 'should NOT increase rate limit' () {

        given:
        def change=0
        def opts = new ThrottlingExecutor.Options()
                    .onRateLimitChange { change++ }
                    .withRampUp(5, 2)
                    .normalise()

        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(4)

        when:
        exec.speedUp0(limiter)
        then:
        limiter.rate == 5 // <- use the max rate define in the options
        change == 1
    }


    def 'should decrease submit rate' () {
        given:
        def change=0
        def opts = new ThrottlingExecutor.Options()
                        .onRateLimitChange { change++ }
                        .normalise()
        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(1.0)

        when:
        exec.backOff0(null)
        then:
        limiter.rate == 1.0d
        change == 0

        when:
        exec.backOff0(limiter)
        then:
        limiter.rate == 0.5
        change == 1 
    }

    def 'should not decrease submit rate' () {
        given:
        def change=0
        def opts = new ThrottlingExecutor.Options()
                        .onRateLimitChange { change++ }
                        .withBackOff(0.75)
                        .normalise()
        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(1.0)

        when:
        exec.backOff0(limiter)
        then:
        limiter.rate == 0.75    // <-- use the min rate define in the setting
        change == 1

    }


    def 'should throttle requests' () {

        given:
        def opts = new ThrottlingExecutor.Options().withRateLimit('1/sec')
        def exec = ThrottlingExecutor.create(opts)
        long begin = System.currentTimeMillis()

        when:
        List<Future> futures = []
        int count = 0
        for( int i=0; i<5; i++ ) {
            futures << exec.submit( { log.info "tick=${count++}" } )
        }

        futures.each {  while(!it.isDone()) sleep(500) }
        long delta = System.currentTimeMillis()-begin
        println "delta=$delta"
        then:
        delta > 3_000
        delta < 5_000
        count == 5
    }

    def 'should apply request priority' () {

        given:
        def opts = new ThrottlingExecutor.Options().withPoolSize(1).withQueueSize(10)
        def exec = ThrottlingExecutor.create(opts)

        when:
        int seq = 0
        def f1 = exec.submit( (Callable)new ThrottlingExecutor.Recoverable() {
            @Override protected Object invoke() { sleep 500; return ++seq }
        })
        def f2 = exec.submit( (Callable)new ThrottlingExecutor.Recoverable() {
            @Override protected Object invoke() { sleep 500; return ++seq }
        })
        def f3 = exec.submit( (Callable)new ThrottlingExecutor.Recoverable() {
            @Override protected Object invoke() { sleep 500; return ++seq }
            @Override byte getPriority() { 10 as byte }
        })
        def f4 = exec.submit( (Callable)new ThrottlingExecutor.Recoverable() {
            @Override protected Object invoke() { sleep 500; return ++seq }
            @Override byte getPriority() { 20 as byte }
        })
        def f5 = exec.submit( (Callable)new ThrottlingExecutor.Recoverable() {
            @Override protected Object invoke() { sleep 500; return ++seq }
            @Override byte getPriority() { 30 as byte }
        })

        then:
        f1.get() == 1
        f5.get() == 2   // f3 is executed as second because it has higher priority (30)
        f4.get() == 3
        f3.get() == 4
        f2.get() == 5   // f2 is executed as last because it has lower priority (0)

    }

    def 'check priority blocking queue ordering' () {
        given:
        def q = new PriorityBlockingQueue<Byte>()
        q.offer( 30 as byte )
        q.offer( 10 as byte )
        q.offer( -10 as byte )
        q.offer( 20 as byte )
        q.offer( 0 as byte )

        expect: 
        q.take() == -10
        q.take() == 0
        q.take() == 10
        q.take() == 20
        q.take() == 30
    }

}
