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
 *
 */

package io.seqera.executor

import java.util.concurrent.atomic.AtomicInteger

import io.seqera.sched.api.schema.v1a1.CreateTasksResponse
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.client.SchedClient
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Timeout

/**
 * Tests for SeqeraBatchSubmitter
 *
 * @author Lorenzo Fontana <fontanalorenz@gmail.com>
 */
@Timeout(30)
class SeqeraBatchSubmitterTest extends Specification {

    def 'should batch multiple tasks submitted within the interval'() {
        given:
        def taskIds = ['task-1', 'task-2', 'task-3']
        def capturedTasks = []
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                capturedTasks.addAll(args[0])
                Stub(CreateTasksResponse) {
                    getTaskIds() >> taskIds.take(args[0].size())
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('500ms'))
        def handlers = (1..3).collect { createMockHandler() }
        def tasks = (1..3).collect { new Task().sessionId("session-$it") }

        when: 'start submitter and enqueue tasks quickly'
        submitter.start()
        handlers.eachWithIndex { handler, i ->
            submitter.enqueue(handler, tasks[i])
        }
        // Wait for the batch interval to elapse and tasks to be submitted
        sleep(800)
        submitter.shutdown()

        then: 'all tasks should be submitted in a single batch'
        capturedTasks.size() == 3
        capturedTasks*.sessionId == ['session-1', 'session-2', 'session-3']

        and: 'task IDs should be assigned to handlers'
        handlers.each { handler ->
            1 * handler.setBatchTaskId(_)
        }
    }

    def 'should submit tasks in separate batches when interval elapses between them'() {
        given:
        def batchCount = new AtomicInteger()
        def batches = [].asSynchronized()  // thread-safe list to capture batches
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                def batchNum = batchCount.incrementAndGet()
                def tasksInBatch = args[0].collect { it.sessionId }
                batches << [batch: batchNum, tasks: tasksInBatch]
                Stub(CreateTasksResponse) {
                    getTaskIds() >> args[0].collect { "task-${it.sessionId}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('200ms'))

        when: 'enqueue first batch, wait for flush, then enqueue second batch'
        submitter.start()
        // First batch - tasks s1 and s2
        submitter.enqueue(createMockHandler(), new Task().sessionId('s1'))
        submitter.enqueue(createMockHandler(), new Task().sessionId('s2'))
        // Wait for first batch to flush (interval + buffer)
        sleep(400)
        // Second batch - task s3
        submitter.enqueue(createMockHandler(), new Task().sessionId('s3'))
        // Wait for second batch to flush
        sleep(400)
        submitter.shutdown()

        then: 'should have two separate batches with correct tasks'
        batches.size() == 2
        and: 'first batch contains s1 and s2'
        batches[0].batch == 1
        batches[0].tasks == ['s1', 's2']
        and: 'second batch contains only s3'
        batches[1].batch == 2
        batches[1].tasks == ['s3']
    }

    def 'should flush batch immediately when reaching max size'() {
        given:
        def batchSizes = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                batchSizes << args[0].size()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> args[0].collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('10s')) // Long interval

        when: 'enqueue exactly TASKS_PER_REQUEST tasks'
        submitter.start()
        (1..SeqeraBatchSubmitter.TASKS_PER_REQUEST).each {
            submitter.enqueue(createMockHandler(), new Task().sessionId("s$it"))
        }
        // Give a small delay for the batch to be processed
        sleep(200)
        submitter.shutdown()

        then: 'should flush immediately due to batch size limit'
        batchSizes.size() >= 1
        batchSizes[0] == SeqeraBatchSubmitter.TASKS_PER_REQUEST
    }

    def 'should split into multiple batches when exceeding max size'() {
        given:
        def totalTasks = SeqeraBatchSubmitter.TASKS_PER_REQUEST + 20  // 120 tasks
        def batchSizes = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                batchSizes << args[0].size()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> args[0].collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('10s')) // Long interval

        when: 'enqueue more than TASKS_PER_REQUEST tasks'
        submitter.start()
        (1..totalTasks).each {
            submitter.enqueue(createMockHandler(), new Task().sessionId("s$it"))
        }
        // Wait for batches to be processed
        sleep(500)
        submitter.shutdown()

        then: 'should create two batches: 100 + 20'
        batchSizes.size() == 2
        batchSizes[0] == SeqeraBatchSubmitter.TASKS_PER_REQUEST  // 100
        batchSizes[1] == 20
    }

    def 'should flush remaining tasks on shutdown'() {
        given:
        def capturedTasks = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                capturedTasks.addAll(args[0])
                Stub(CreateTasksResponse) {
                    getTaskIds() >> args[0].collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('10s')) // Long interval

        when: 'enqueue tasks and immediately shutdown'
        submitter.start()
        submitter.enqueue(createMockHandler(), new Task().sessionId('s1'))
        submitter.enqueue(createMockHandler(), new Task().sessionId('s2'))
        // Shutdown without waiting for interval
        submitter.shutdown()

        then: 'tasks should still be submitted'
        capturedTasks.size() == 2
    }

    def 'should throw exception when enqueueing after shutdown'() {
        given:
        def client = Stub(SchedClient)
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('1s'))

        when:
        submitter.start()
        submitter.shutdown()
        submitter.enqueue(createMockHandler(), new Task().sessionId('s1'))

        then:
        thrown(IllegalStateException)
    }

    def 'should propagate failure to handlers on API error'() {
        given:
        def apiError = new RuntimeException('API error')
        def client = Stub(SchedClient) {
            createTasks(_) >> { throw apiError }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('100ms'))
        def handler1 = Mock(SeqeraTaskHandler)
        def handler2 = Mock(SeqeraTaskHandler)

        when:
        submitter.start()
        submitter.enqueue(handler1, new Task().sessionId('s1'))
        submitter.enqueue(handler2, new Task().sessionId('s2'))
        sleep(300)
        submitter.shutdown()

        then: 'all handlers should receive the failure'
        1 * handler1.onBatchSubmitFailure(apiError)
        1 * handler2.onBatchSubmitFailure(apiError)
    }

    def 'should propagate failure to handlers when API returns wrong number of task IDs'() {
        given:
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                Stub(CreateTasksResponse) {
                    // Return fewer IDs than tasks submitted
                    getTaskIds() >> ['task-1']
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('100ms'))
        def handler1 = Mock(SeqeraTaskHandler)
        def handler2 = Mock(SeqeraTaskHandler)
        def handler3 = Mock(SeqeraTaskHandler)

        when:
        submitter.start()
        submitter.enqueue(handler1, new Task().sessionId('s1'))
        submitter.enqueue(handler2, new Task().sessionId('s2'))
        submitter.enqueue(handler3, new Task().sessionId('s3'))
        sleep(300)
        submitter.shutdown()

        then: 'handlers should receive failure due to mismatched IDs'
        1 * handler1.onBatchSubmitFailure({ it instanceof IllegalStateException })
        1 * handler2.onBatchSubmitFailure({ it instanceof IllegalStateException })
        1 * handler3.onBatchSubmitFailure({ it instanceof IllegalStateException })
    }

    def 'should start batch timer only when first task arrives not when thread starts'() {
        given:
        def submitTimes = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_) >> { List<List<Task>> args ->
                submitTimes << System.currentTimeMillis()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> args[0].collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, Duration.of('300ms'))

        when: 'start submitter, wait longer than interval, then enqueue tasks'
        submitter.start()
        // Wait longer than the batch interval before enqueueing
        sleep(500)
        def enqueueTime = System.currentTimeMillis()
        submitter.enqueue(createMockHandler(), new Task().sessionId('s1'))
        submitter.enqueue(createMockHandler(), new Task().sessionId('s2'))
        // Wait for batch to flush
        sleep(500)
        submitter.shutdown()

        then: 'batch should be submitted ~300ms after enqueue, not immediately'
        submitTimes.size() == 1
        // The submit should happen ~300ms after enqueue
        def timeSinceEnqueue = submitTimes[0] - enqueueTime
        timeSinceEnqueue >= 250 // Allow some tolerance
        timeSinceEnqueue < 600  // But not too long
    }

    def 'should use default interval when not specified'() {
        given:
        def client = Stub(SchedClient)

        when:
        def submitter = new SeqeraBatchSubmitter(client)

        then:
        submitter.requestInterval == SeqeraBatchSubmitter.REQUEST_INTERVAL
    }

    /**
     * Creates a mock handler that can track setBatchTaskId and onBatchSubmitFailure calls
     */
    private SeqeraTaskHandler createMockHandler() {
        Mock(SeqeraTaskHandler)
    }
}
