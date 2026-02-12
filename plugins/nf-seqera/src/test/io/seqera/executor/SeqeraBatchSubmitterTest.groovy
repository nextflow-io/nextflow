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

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger

import io.seqera.sched.api.schema.v1a1.CreateTasksResponse
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.client.SchedClient
import nextflow.file.FileHolder
import nextflow.processor.TaskRun
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Timeout

/**
 * Tests for SeqeraBatchSubmitter
 *
 * @author Lorenzo Fontana <fontanalorenz@gmail.com>
 */
@Timeout(30)
class SeqeraBatchSubmitterTest extends Specification {

    static final String TEST_RUN = 'run-test:local'

    def 'should batch multiple tasks submitted within the interval'() {
        given:
        def taskIds = ['task-1', 'task-2', 'task-3']
        def capturedTasks = []
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                capturedTasks.addAll(tasks)
                Stub(CreateTasksResponse) {
                    getTaskIds() >> taskIds.take(tasks.size())
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('500ms'))
        def handlers = (1..3).collect { createMockHandler() }
        def tasks = (1..3).collect { new Task().image("image-$it") }

        when: 'start submitter and enqueue tasks quickly'
        submitter.start()
        handlers.eachWithIndex { handler, i ->
            submitter.submit(handler, tasks[i])
        }
        // Wait for the batch interval to elapse and tasks to be submitted
        sleep(800)
        submitter.shutdown()

        then: 'all tasks should be submitted in a single batch'
        capturedTasks.size() == 3
        capturedTasks*.image == ['image-1', 'image-2', 'image-3']

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
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                def batchNum = batchCount.incrementAndGet()
                def tasksInBatch = tasks.collect { it.image }
                batches << [batch: batchNum, tasks: tasksInBatch]
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${it.image}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('200ms'))

        when: 'enqueue first batch, wait for flush, then enqueue second batch'
        submitter.start()
        // First batch - tasks s1 and s2
        submitter.submit(createMockHandler(), new Task().image('s1'))
        submitter.submit(createMockHandler(), new Task().image('s2'))
        // Wait for first batch to flush (interval + buffer)
        sleep(400)
        // Second batch - task s3
        submitter.submit(createMockHandler(), new Task().image('s3'))
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
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                batchSizes << tasks.size()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('10s')) // Long interval

        when: 'enqueue exactly TASKS_PER_REQUEST tasks'
        submitter.start()
        (1..SeqeraBatchSubmitter.TASKS_PER_REQUEST).each {
            submitter.submit(createMockHandler(), new Task().image("img$it"))
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
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                batchSizes << tasks.size()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('10s')) // Long interval

        when: 'enqueue more than TASKS_PER_REQUEST tasks'
        submitter.start()
        (1..totalTasks).each {
            submitter.submit(createMockHandler(), new Task().image("img$it"))
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
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                capturedTasks.addAll(tasks)
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('10s')) // Long interval

        when: 'enqueue tasks and immediately shutdown'
        submitter.start()
        submitter.submit(createMockHandler(), new Task().image('img1'))
        submitter.submit(createMockHandler(), new Task().image('img2'))
        // Shutdown without waiting for interval
        submitter.shutdown()

        then: 'tasks should still be submitted'
        capturedTasks.size() == 2
    }

    def 'should throw exception when enqueueing after shutdown'() {
        given:
        def client = Stub(SchedClient)
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('1s'))

        when:
        submitter.start()
        submitter.shutdown()
        submitter.submit(createMockHandler(), new Task().image('img1'))

        then:
        thrown(IllegalStateException)
    }

    def 'should propagate failure to handlers on API error'() {
        given:
        def apiError = new RuntimeException('API error')
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { throw apiError }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'))
        def handler1 = Mock(SeqeraTaskHandler)
        def handler2 = Mock(SeqeraTaskHandler)

        when:
        submitter.start()
        submitter.submit(handler1, new Task().image('img1'))
        submitter.submit(handler2, new Task().image('img2'))
        sleep(300)
        submitter.shutdown()

        then: 'all handlers should receive the failure'
        1 * handler1.onBatchSubmitFailure(apiError)
        1 * handler2.onBatchSubmitFailure(apiError)
    }

    def 'should propagate failure to handlers when API returns wrong number of task IDs'() {
        given:
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                Stub(CreateTasksResponse) {
                    // Return fewer IDs than tasks submitted
                    getTaskIds() >> ['task-1']
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'))
        def handler1 = Mock(SeqeraTaskHandler)
        def handler2 = Mock(SeqeraTaskHandler)
        def handler3 = Mock(SeqeraTaskHandler)

        when:
        submitter.start()
        submitter.submit(handler1, new Task().image('img1'))
        submitter.submit(handler2, new Task().image('img2'))
        submitter.submit(handler3, new Task().image('img3'))
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
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                submitTimes << System.currentTimeMillis()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('300ms'))

        when: 'start submitter, wait longer than interval, then enqueue tasks'
        submitter.start()
        // Wait longer than the batch interval before enqueueing
        sleep(500)
        def enqueueTime = System.currentTimeMillis()
        submitter.submit(createMockHandler(), new Task().image('img1'))
        submitter.submit(createMockHandler(), new Task().image('img2'))
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
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN)

        then:
        submitter.requestInterval == SeqeraBatchSubmitter.REQUEST_INTERVAL
        submitter.keepAliveInterval == SeqeraBatchSubmitter.KEEP_ALIVE_INTERVAL
    }

    def 'should send keep-alive when no tasks received within keep-alive interval'() {
        given:
        def submissions = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                submissions << [runId: runId, taskCount: tasks.size()]
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        // Short keep-alive interval for testing
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('10s'), Duration.of('200ms'))

        when: 'start submitter and wait for keep-alive interval without enqueueing tasks'
        submitter.start()
        // Wait for keep-alive to trigger (interval + buffer)
        sleep(400)
        submitter.shutdown()

        then: 'should have sent at least one keep-alive (empty submission)'
        submissions.size() >= 1
        submissions.any { it.taskCount == 0 }
    }

    def 'should not send keep-alive when tasks are being submitted regularly'() {
        given:
        def submissions = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                submissions << [runId: runId, taskCount: tasks.size()]
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        // Short intervals for testing
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'), Duration.of('300ms'))

        when: 'start submitter and continuously enqueue tasks'
        submitter.start()
        // Enqueue tasks at intervals shorter than keep-alive
        5.times { i ->
            submitter.submit(createMockHandler(), new Task().image("img$i"))
            sleep(80)
        }
        sleep(200) // Wait for final batch to flush
        submitter.shutdown()

        then: 'all submissions should have tasks, no keep-alive (empty) submissions'
        submissions.size() >= 1
        submissions.every { it.taskCount > 0 }
    }

    def 'should invoke error callback on fatal error in thread'() {
        given:
        // Use an Error to simulate a truly fatal condition that escapes flushBatch
        def fatalError = new OutOfMemoryError('Fatal thread error')
        def errorReceived = null
        def onError = { Throwable t -> errorReceived = t }
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                // Throw a fatal Error that will escape flushBatch
                throw fatalError
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'), Duration.of('10s'), onError)
        def handler = Mock(SeqeraTaskHandler)

        when:
        submitter.start()
        submitter.submit(handler, new Task().image('img1'))
        sleep(300)
        submitter.shutdown()

        then: 'handler should receive failure (wrapped in RuntimeException since it was an Error)'
        1 * handler.onBatchSubmitFailure({ it instanceof RuntimeException && it.cause == fatalError })
        and: 'error callback should be invoked with original error'
        errorReceived == fatalError
    }

    def 'should drain and fail pending tasks on fatal error'() {
        given:
        def fatalError = new RuntimeException('Fatal thread error')
        def failedHandlers = [].asSynchronized()
        def onError = { Throwable t -> /* no-op */ }
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                // Always throw to simulate fatal error
                throw fatalError
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'), Duration.of('10s'), onError)

        // Create handlers that track when they receive failure notification
        def handlers = (1..3).collect {
            def h = Mock(SeqeraTaskHandler) {
                onBatchSubmitFailure(_) >> { args -> failedHandlers << it }
            }
            h
        }

        when:
        submitter.start()
        handlers.each { h ->
            submitter.submit(h, new Task().image('img'))
        }
        sleep(300)
        submitter.shutdown()

        then: 'all handlers should receive failure notification'
        handlers.each { h ->
            1 * h.onBatchSubmitFailure(fatalError)
        }
    }

    def 'should continue after keep-alive failure'() {
        given:
        def keepAliveFailCount = new AtomicInteger()
        def taskSubmissions = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                if (tasks.isEmpty()) {
                    // Keep-alive call - fail it
                    keepAliveFailCount.incrementAndGet()
                    throw new RuntimeException('Keep-alive failed')
                }
                // Normal task submission
                taskSubmissions << tasks.size()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        // Short keep-alive interval to trigger failures quickly
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('500ms'), Duration.of('100ms'))

        when: 'start submitter, let keep-alive fail, then submit tasks'
        submitter.start()
        // Wait for a few keep-alive failures
        sleep(350)
        // Now submit a task
        submitter.submit(createMockHandler(), new Task().image('img1'))
        // Wait for task to be submitted
        sleep(700)
        submitter.shutdown()

        then: 'keep-alive should have failed but thread should continue'
        keepAliveFailCount.get() >= 1
        and: 'task submission should still succeed'
        taskSubmissions.size() >= 1
    }

    def 'should work without error callback'() {
        given:
        def taskSubmissions = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                taskSubmissions << tasks.size()
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        // Constructor without error callback (null)
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'), Duration.of('10s'), null)

        when:
        submitter.start()
        submitter.submit(createMockHandler(), new Task().image('img1'))
        sleep(300)
        submitter.shutdown()

        then: 'should work normally'
        taskSubmissions.size() >= 1
    }

    // -- input files metrics tests --

    @TempDir
    Path tempDir

    def 'should attach input files metrics to task before submission'() {
        given:
        def file1 = tempDir.resolve('a.txt')
        Files.write(file1, new byte[500])
        def file2 = tempDir.resolve('b.txt')
        Files.write(file2, new byte[3000])
        def capturedTasks = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                capturedTasks.addAll(tasks)
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'))
        def taskRun = Mock(TaskRun) {
            getInputFiles() >> [new FileHolder(file1), new FileHolder(file2)]
        }
        def handler = Mock(SeqeraTaskHandler) {
            getTask() >> taskRun
        }

        when:
        submitter.start()
        submitter.submit(handler, new Task().image('img1'))
        sleep(300)
        submitter.shutdown()

        then: 'metrics should be set on the submitted task'
        capturedTasks.size() == 1
        def metrics = capturedTasks[0].getInputFiles()
        metrics != null
        metrics.count == 2
        metrics.totalBytes == 3500
        metrics.maxFileBytes == 3000
        metrics.minFileBytes == 500
    }

    def 'should submit task without metrics when no input files'() {
        given:
        def capturedTasks = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                capturedTasks.addAll(tasks)
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'))
        def taskRun = Mock(TaskRun) {
            getInputFiles() >> []
        }
        def handler = Mock(SeqeraTaskHandler) {
            getTask() >> taskRun
        }

        when:
        submitter.start()
        submitter.submit(handler, new Task().image('img1'))
        sleep(300)
        submitter.shutdown()

        then: 'task should be submitted without metrics'
        capturedTasks.size() == 1
        capturedTasks[0].getInputFiles() == null
    }

    def 'should submit task even when metrics computation fails'() {
        given:
        def capturedTasks = [].asSynchronized()
        def client = Stub(SchedClient) {
            createTasks(_, _) >> { String runId, List<Task> tasks ->
                capturedTasks.addAll(tasks)
                Stub(CreateTasksResponse) {
                    getTaskIds() >> tasks.collect { "task-${System.nanoTime()}" }
                }
            }
        }
        def submitter = new SeqeraBatchSubmitter(client, TEST_RUN, Duration.of('100ms'))
        def taskRun = Mock(TaskRun) {
            getInputFiles() >> { throw new RuntimeException('Cannot access input files') }
        }
        def handler = Mock(SeqeraTaskHandler) {
            getTask() >> taskRun
        }

        when:
        submitter.start()
        submitter.submit(handler, new Task().image('img1'))
        sleep(300)
        submitter.shutdown()

        then: 'task should still be submitted despite metrics failure'
        capturedTasks.size() == 1
    }

    /**
     * Creates a mock handler that can track setBatchTaskId and onBatchSubmitFailure calls
     */
    private SeqeraTaskHandler createMockHandler() {
        Mock(SeqeraTaskHandler)
    }
}
