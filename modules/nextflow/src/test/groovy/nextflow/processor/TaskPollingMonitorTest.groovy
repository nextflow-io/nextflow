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


import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.ExecutorConfig
import nextflow.util.Duration
import nextflow.util.MemoryUnit
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


    def 'should create a resource account from executor config'() {
        when:
        def config = new ExecutorConfig([:])
        then:
        TaskPollingMonitor.resourceAccount(config, 'slurm').isUnlimited()

        when:
        config = new ExecutorConfig(cpus: 8, memory: '16GB')
        def account = TaskPollingMonitor.resourceAccount(config, 'slurm')
        then:
        account.maxCpus == 8
        account.maxMemory == MemoryUnit.of('16GB').toBytes()

        when: 'the setting is scoped to a specific executor'
        config = new ExecutorConfig(cpus: 8, '$slurm': [cpus: 4])
        account = TaskPollingMonitor.resourceAccount(config, 'slurm')
        then:
        account.maxCpus == 4

        when: 'a negative value is specified'
        TaskPollingMonitor.resourceAccount(new ExecutorConfig(cpus: -1), 'slurm')
        then:
        thrown(AssertionError)

        when: 'a negative memory value is specified'
        TaskPollingMonitor.resourceAccount(new ExecutorConfig(memory: -1), 'slurm')
        then:
        thrown(AssertionError)
    }

    def 'should limit the submission based on executor cpus and memory'() {
        given:
        def session = Mock(Session)
        def account = new ResourceAccount(8, MemoryUnit.of('16GB').toBytes())
        def monitor = Spy(new TaskPollingMonitor(name: 'foo', session: session, pollInterval: Duration.of('1min'), resourceAccount: account))
        and:
        def handler = Mock(TaskHandler) {
            getTask() >> new TaskRun(config: new TaskConfig(cpus: 4, memory: MemoryUnit.of('8GB')))
            canForkProcess() >> true
            isReady() >> true
        }

        expect:
        monitor.canSubmit(handler)

        when: 'a first task is submitted'
        monitor.submit(handler)
        then:
        1 * handler.prepareLauncher()
        1 * handler.submit()
        account.availableCpus() == 4
        account.availableMemory() == MemoryUnit.of('8GB').toBytes()

        when: 'a second task is requested'
        then:
        monitor.canSubmit(handler)

        when: 'the second task is submitted'
        monitor.submit(handler)
        then:
        account.availableCpus() == 0
        account.availableMemory() == 0

        when: 'a third task is requested'
        then:
        !monitor.canSubmit(handler)

        when: 'a task completes'
        monitor.remove(handler)
        then:
        account.availableCpus() == 4
        account.availableMemory() == MemoryUnit.of('8GB').toBytes()
        monitor.canSubmit(handler)
    }

    @Unroll
    def 'should fail when a task requirement exceeds the executor resources'() {
        given:
        def session = Mock(Session)
        def monitor = Spy(new TaskPollingMonitor(name: 'foo', session: session, pollInterval: Duration.of('1min'), resourceAccount: new ResourceAccount(8, MemoryUnit.of('16GB').toBytes())))
        and:
        def handler = Mock(TaskHandler) {
            getTask() >> new TaskRun(config: new TaskConfig(cpus: CPUS, memory: MemoryUnit.of(MEMORY)))
            canForkProcess() >> true
            isReady() >> true
        }

        when:
        monitor.canSubmit(handler)
        then:
        def e = thrown(ProcessUnrecoverableException)
        e.message.contains(EXPECTED)

        where:
        CPUS | MEMORY  | EXPECTED
        10   | '8GB'   | 'Process requirement exceeds available CPUs -- req: 10; avail: 8'
        4    | '20GB'  | 'Process requirement exceeds available memory -- req: 20 GB; avail: 16 GB'
    }

    def 'should reserve the total resources of a job array'() {
        given:
        def session = Mock(Session)
        def account = new ResourceAccount(8, MemoryUnit.of('16GB').toBytes())
        def monitor = Spy(new TaskPollingMonitor(name: 'foo', session: session, pollInterval: Duration.of('1min'), resourceAccount: account))
        and:
        def children = (1..3).collect {
            Mock(TaskHandler) {
                getTask() >> new TaskRun(config: new TaskConfig(cpus: 2, memory: MemoryUnit.of('2GB')))
            }
        }
        def arrayHandler = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) { getChildren() >> children }
            canForkProcess() >> true
            isReady() >> true
        }

        when:
        monitor.submit(arrayHandler)
        then:
        1 * arrayHandler.prepareLauncher()
        1 * arrayHandler.submit()
        account.availableCpus() == 2
        account.availableMemory() == MemoryUnit.of('10GB').toBytes()

        when: 'one child completes'
        monitor.remove(children[0])
        then: 'the resources are still reserved'
        account.availableCpus() == 2
        account.availableMemory() == MemoryUnit.of('10GB').toBytes()

        when: 'all children complete'
        monitor.remove(children[1])
        monitor.remove(children[2])
        then:
        account.availableCpus() == 8
        account.availableMemory() == MemoryUnit.of('16GB').toBytes()
    }

}
