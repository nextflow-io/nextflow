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
 */

package nextflow.processor


import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

import nextflow.Session
import nextflow.util.Duration
import nextflow.util.ThrottlingExecutor
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParallelPollingMonitorTest extends Specification {

    def 'should retry task' () {

        given:
        def success = new AtomicInteger()
        def failure = new AtomicInteger()
        def count = new AtomicInteger()
        def change = new AtomicInteger()
        def retry = new AtomicInteger()

        def session = Mock(Session)
        def handler = Mock(TaskHandler)

        def opts = new ThrottlingExecutor.Options()
                        .retryOn(IllegalArgumentException)
                        .withRateLimit('10/sec')
                        .withErrorBurstDelay(Duration.of('5sec'))
                        .withAutoThrottle()
                        .onSuccess { success.incrementAndGet() }
                        .onRetry { retry.incrementAndGet() }
                        .onFailure { failure.incrementAndGet() }
                        .onRateLimitChange { change.incrementAndGet() }

        def exec = ThrottlingExecutor.create(opts)

        def mon = new ParallelPollingMonitor(exec, [session:session, name:'foo', pollInterval:'1sec'])
        when:
        mon.submitter = exec
        mon.submit(handler)
        sleep 3_000
        exec.shutdown()
        exec.awaitTermination(1, TimeUnit.MINUTES)

        then:
        handler.submit() >> { println "c=$count"; if(count.getAndIncrement()<2) throw new IllegalArgumentException("Ooops!") }

        success.get() == 1
        retry.get() == 2
        failure.get() == 0
        change.get() == 1

    }

    @Unroll
    def 'should validate can submit method' () {
        given:
        def success = new AtomicInteger()
        def failure = new AtomicInteger()
        def count = new AtomicInteger()
        def change = new AtomicInteger()
        def retry = new AtomicInteger()

        def session = Mock(Session)
        def handler = Mock(TaskHandler)

        def opts = new ThrottlingExecutor.Options()
            .retryOn(IllegalArgumentException)
            .withRateLimit('10/sec')
            .withErrorBurstDelay(Duration.of('5sec'))
            .withAutoThrottle()
            .onSuccess { success.incrementAndGet() }
            .onRetry { retry.incrementAndGet() }
            .onFailure { failure.incrementAndGet() }
            .onRateLimitChange { change.incrementAndGet() }

        def exec = ThrottlingExecutor.create(opts)
        def mon = Spy(new ParallelPollingMonitor(exec, [capacity:CAPACITY, session:session, name:'foo', pollInterval:'1sec']))
        and:
        SUBMIT.times { mon.runningQueue.add(Mock(TaskHandler))  }

        when:
        def result = mon.canSubmit(handler)
        then:
        handler.canForkProcess() >> FORK
        and:
        result == EXPECTED

        where:
        CAPACITY    | SUBMIT | FORK  | EXPECTED
        0           |   5    | true  | true
        10          |   5    | true  | true
        10          |   10   | true  | false
        and:
        0           |   1    | false | false
        10          |   1    | false | false
    }

}
