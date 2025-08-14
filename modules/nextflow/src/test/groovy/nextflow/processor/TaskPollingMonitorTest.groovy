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


import nextflow.Session
import nextflow.executor.ExecutorConfig
import nextflow.util.Duration
import nextflow.util.RateUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskPollingMonitorTest extends Specification {

    def testCreate() {

        setup:
        def name = 'hello'
        def session = Mock(Session)
        def config = new ExecutorConfig(pollInterval: '1h', queueSize: 11, dumpInterval: '3h')

        def defSize = 99
        def defPollDuration = Duration.of('44s')
        when:
        def monitor = TaskPollingMonitor.create(session, config, name, defSize, defPollDuration)
        then:
        monitor.name == 'hello'
        monitor.pollIntervalMillis == Duration.of('1h').toMillis()
        monitor.capacity == 11
        monitor.dumpInterval ==  Duration.of('3h')

    }

    def 'should create a rate limiter for the given rate format'() {

        given:
        def session = Mock(Session)
        def config = new ExecutorConfig(submitRateLimit: RATE)
        def monitor = new TaskPollingMonitor(name:'local', session: session, config: config, pollInterval: '1s', capacity: 100)
        
        when:
        def limit = monitor.createSubmitRateLimit()
        then:
        limit ? Math.round(limit.getRate()) : null == EXPECTED

        where:
        RATE            | EXPECTED
        '1'             | 1                         // 1 per second
        '5'             | 5                         // 5 per second
        '100 min'       | 100 / 60                  // 100 per minute
        '100 / 1 s'     | 100                       // 100 per second
        '100 / 2 s'     | 50                        // 100 per 2 seconds
        '200 / sec'     | 200                       // 200 per second
        '600 / 5'       | 600i / 5l as double       // 600 per 5 seconds
        '600 / 5min'    | 600 / (5 * 60)            // 600 per 5 minutes

    }


    def 'check equals and hash code' () {
        expect:
        new RateUnit(2.1) == new RateUnit(2.1)
        new RateUnit(2.1).hashCode() == new RateUnit(2.1).hashCode()
        new RateUnit(2.1) != new RateUnit(3.3)
        new RateUnit(2.1).hashCode() != new RateUnit(3.3).hashCode()
    }

    def 'should stringify' () {
        expect:
        new RateUnit(0.1).toString() == '0.10/sec'
        new RateUnit(2.1).toString() == '2.10/sec'
        new RateUnit(123.4).toString() == '123.40/sec'
    }

    def 'should cancel running jobs' () {
        given:
        def session = Mock(Session)
        def monitor = new TaskPollingMonitor(name:'foo', session: session, pollInterval: Duration.of('1min'))
        def spy = Spy(monitor)
        and:
        def handler = Mock(TaskHandler) { getTask() >> Mock(TaskRun) }
        and:
        spy.submit(handler)

        when:
        spy.cleanup()
        then:
        1 * session.disableJobsCancellation >> false
        and:
        1 * handler.kill() >> null
        1 * session.notifyTaskComplete(handler) >> null
    }

    def 'should not cancel running jobs' () {
        given:
        def session = Mock(Session)
        def monitor = new TaskPollingMonitor(name:'foo', session: session, pollInterval: Duration.of('1min'))
        def spy = Spy(monitor)
        and:
        def handler = Mock(TaskHandler) { getTask() >> Mock(TaskRun) }
        and:
        spy.submit(handler)

        when:
        spy.cleanup()
        then:
        1 * session.disableJobsCancellation >> true
        and:
        0 * handler.killTask() >> null
        0 * session.notifyTaskComplete(handler) >> null
    }

    def 'should submit a job array' () {
        given:
        def session = Mock(Session)
        def monitor = Spy(new TaskPollingMonitor(name: 'foo', session: session, pollInterval: Duration.of('1min')))
        and:
        def handler = Mock(TaskHandler) {
            getTask() >> Mock(TaskRun)
        }
        def arrayHandler = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) {
                children >> (1..3).collect( i -> handler )
            }
        }

        when:
        monitor.submit(arrayHandler)
        then:
        1 * arrayHandler.prepareLauncher()
        1 * arrayHandler.submit()
        0 * handler.prepareLauncher()
        0 * handler.submit()
        3 * session.notifyTaskSubmit(handler)
    }

    @Unroll
    def 'should throw error if job array size exceeds queue size [capacity: #CAPACITY, array: #ARRAY_SIZE]' () {
        given:
        def session = Mock(Session)
        def monitor = Spy(new TaskPollingMonitor(name: 'foo', session: session, capacity: CAPACITY, pollInterval: Duration.of('1min')))
        and:
        def processor = Mock(TaskProcessor)
        def arrayHandler = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) {
                getName() >> TASK_NAME
                getArraySize() >> ARRAY_SIZE
                getProcessor() >> processor
            }
        }

        when:
        monitor.canSubmit(arrayHandler)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Process '$TASK_NAME' declares array size ($ARRAY_SIZE) which exceeds the executor queue size ($CAPACITY)")

        where:
        CAPACITY | ARRAY_SIZE | TASK_NAME
        10       | 15         | 'test_array'
        5        | 10         | 'large_array'  
        1        | 2          | 'small_array'
    }

    @Unroll
    def 'should validate array size accounting in queue capacity [capacity: #CAPACITY, running: #RUNNING_COUNT, array: #ARRAY_SIZE]' () {
        given:
        def session = Mock(Session)
        def monitor = Spy(new TaskPollingMonitor(name: 'foo', session: session, capacity: CAPACITY, pollInterval: Duration.of('1min')))
        and:
        def processor = Mock(TaskProcessor)
        def regularHandler = Mock(TaskHandler) {
            getTask() >> Mock(TaskRun) {
                getProcessor() >> processor
            }
            canForkProcess() >> CAN_FORK
            isReady() >> IS_READY
        }
        def arrayHandler = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) {
                getArraySize() >> ARRAY_SIZE
                getProcessor() >> processor
            }
            canForkProcess() >> CAN_FORK
            isReady() >> IS_READY
        }

        and:
        RUNNING_COUNT.times { monitor.runningQueue.add(regularHandler) }

        expect:
        monitor.runningQueue.size() == RUNNING_COUNT
        monitor.canSubmit(regularHandler) == REGULAR_EXPECTED
        monitor.canSubmit(arrayHandler) == ARRAY_EXPECTED

        where:
        CAPACITY | RUNNING_COUNT | ARRAY_SIZE | CAN_FORK | IS_READY | REGULAR_EXPECTED | ARRAY_EXPECTED
        10       | 6             | 5          | true     | true     | true             | false     // Array too big (6+5>10)
        10       | 5             | 5          | true     | true     | true             | true      // Array fits exactly (5+5=10)
        10       | 9             | 1          | true     | true     | true             | true      // Both fit (9+1=10)
        10       | 10            | 1          | true     | true     | false            | false     // Queue full (10+1>10)
        5        | 4             | 1          | true     | true     | true             | true      // Both fit (4+1=5)
        5        | 4             | 2          | true     | true     | true             | false     // Array too big (4+2>5)
        0        | 5             | 10         | true     | true     | true             | true      // Unlimited capacity
        10       | 5             | 3          | false    | true     | false            | false     // Cannot fork
        10       | 5             | 3          | true     | false    | false            | false     // Not ready
    }

}
