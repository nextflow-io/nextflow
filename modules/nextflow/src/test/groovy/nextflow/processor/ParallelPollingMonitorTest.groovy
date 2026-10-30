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

package nextflow.processor


import java.util.concurrent.ExecutorService
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

import nextflow.Session
import nextflow.exception.ProcessFailedException
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
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
        // the session is healthy in this scenario -- the retries come from the
        // submit operation failing, not from the session being aborted
        _ * session.canSubmitTasks() >> true

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
        handler.isReady() >> true
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

    def 'should not submit a queued task once the session can no longer submit' () {
        given:
        def session = Mock(Session)
        def handler = Mock(TaskHandler)
        and:
        def opts = new ThrottlingExecutor.Options().withRateLimit('100/sec')
        def exec = ThrottlingExecutor.create(opts)
        def mon = new ParallelPollingMonitor(exec, [session:session, name:'foo', pollInterval:'1sec'])

        when:
        mon.submit(handler)
        exec.shutdown()
        exec.awaitTermination(1, TimeUnit.MINUTES)

        then:
        // the session was aborted while the wrapper was sitting in the submitter queue
        _ * session.canSubmitTasks() >> false
        and:
        0 * handler.prepareLauncher()
        0 * handler.submit()
        0 * session.notifyTaskSubmit(_)
    }

    def 'should submit a queued task while the session can still submit' () {
        given:
        def session = Mock(Session)
        def handler = Mock(TaskHandler)
        and:
        def opts = new ThrottlingExecutor.Options().withRateLimit('100/sec')
        def exec = ThrottlingExecutor.create(opts)
        def mon = new ParallelPollingMonitor(exec, [session:session, name:'foo', pollInterval:'1sec'])

        when:
        mon.submit(handler)
        exec.shutdown()
        exec.awaitTermination(1, TimeUnit.MINUTES)

        then:
        _ * session.canSubmitTasks() >> true
        and:
        1 * handler.prepareLauncher()
        1 * handler.submit()
        1 * session.notifyTaskSubmit(handler)
    }

    def 'should inherit array size validation from parent TaskPollingMonitor' () {
        given:
        def session = Mock(Session)
        def throttlingExecutor = Mock(ThrottlingExecutor)
        def monitor = new ParallelPollingMonitor(throttlingExecutor, [name: 'foo', session: session, capacity: 5, pollInterval: Duration.of('1min')])
        and:
        def processor = Mock(TaskProcessor)
        def arrayHandler = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) {
                getName() >> 'oversized_array'
                getArraySize() >> 10
                getProcessor() >> processor
            }
        }

        when:
        monitor.canSubmit(arrayHandler)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Process 'oversized_array' declares array size (10) which exceeds the executor queue size (5)")
    }

    @Unroll
    def 'should release fork slot and submit permit when submission fails [success: #SUCCESS]' () {
        given:
        // success is false when a task was ignored with failOnIgnore, but tasks can still be submitted
        def session = Mock(Session) {
            canSubmitTasks() >> true
            isSuccess() >> SUCCESS
            getExecService() >> Mock(ExecutorService)
        }
        def processor = new TaskProcessorTest.DummyProcessor('foo', session, Mock(BaseScript), new ProcessConfig([maxForks: 1]))
        processor.createStateObj()
        processor.state.update { StateObj it -> it.incSubmitted() }
        and:
        def task = new TaskRun(processor: processor, context: new TaskContext(processor), config: new TaskConfig(errorStrategy: 'ignore', maxErrors: -1))
        def handler = Mock(TaskHandler, constructorArgs: [task]) {
            getTask() >> task
            canForkProcess() >> true
            isReady() >> true
            submit() >> { throw new ProcessFailedException("Error submitting process 'foo' for execution") }
        }
        and:
        def exec = ThrottlingExecutor.create(new ThrottlingExecutor.Options().withRateLimit('100/sec'))
        def mon = new ParallelPollingMonitor(exec, [session: session, name: 'foo', pollInterval: '1sec', capacity: 1])
        mon.getPendingQueue() << handler

        when:
        mon.submitPendingTasks()
        exec.shutdown()
        exec.awaitTermination(1, TimeUnit.MINUTES)
        processor.state.await()
        then:
        task.errorAction == ErrorStrategy.IGNORE
        processor.state.val.completed == 1
        processor.forksCount.sum() == 0
        mon.@semaphore.availablePermits() == 1

        where:
        SUCCESS << [true, false]
    }

}
