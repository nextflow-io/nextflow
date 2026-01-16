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

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.client.SchedClient
import nextflow.SysEnv
import nextflow.util.Duration
import nextflow.util.Threads

/**
 * Batches task submissions to the Seqera scheduler API.
 *
 * @author Lorenzo Fontana <fontanalorenz@gmail.com>
 */
@Slf4j
@CompileStatic
class SeqeraBatchSubmitter {

    /** Maximum tasks per API call */
    static final int TASKS_PER_REQUEST = SysEnv.getInteger('NXF_SEQERA_TASK_PER_REQUEST', 100)

    /** Default flush interval */
    static final Duration REQUEST_INTERVAL = SysEnv.get('NXF_SEQERA_REQUEST_INTERVAL', '1 sec') as Duration

    /**
     * Holds a task handler and its prepared Task object pending submission
     */
    @TupleConstructor
    static class PendingTask {
        SeqeraTaskHandler handler
        Task task
    }

    private final SchedClient client
    private final Duration requestInterval
    private final LinkedBlockingQueue<PendingTask> pendingQueue = new LinkedBlockingQueue<>()
    private Thread sender
    private volatile boolean completed = false

    SeqeraBatchSubmitter(SchedClient client) {
        this(client, REQUEST_INTERVAL)
    }

    SeqeraBatchSubmitter(SchedClient client, Duration requestInterval) {
        this.client = client
        this.requestInterval = requestInterval
    }

    /**
     * Start the sender thread that processes the batch queue
     */
    void start() {
        log.debug "[SEQERA] Starting batch submitter - interval=${requestInterval}"
        this.sender = Threads.start('Seqera-batch-submitter', this.&sendTasks0)
    }

    /**
     * Enqueue a task for batch submission.
     */
    void enqueue(SeqeraTaskHandler handler, Task task) {
        if (completed) {
            throw new IllegalStateException("Batch submitter has been shutdown")
        }
        pendingQueue.add(new PendingTask(handler, task))
    }

    /**
     * Signal completion and wait for sender thread to finish
     */
    void shutdown() {
        log.debug "[SEQERA] Shutting down batch submitter"
        completed = true
        if (sender) {
            sender.join()
        }
        log.debug "[SEQERA] Batch submitter shutdown complete"
    }

    /**
     * Sender thread loop
     */
    protected void sendTasks0(dummy) {
        final List<PendingTask> batch = new ArrayList<>(TASKS_PER_REQUEST)
        long previous = System.currentTimeMillis()
        final long period = requestInterval.millis
        final long delay = period / 10 as long

        while (!completed || !pendingQueue.isEmpty()) {
            // Poll with timeout
            final PendingTask pending = pendingQueue.poll(delay, TimeUnit.MILLISECONDS)
            if (pending) {
                // Start the batch timer when first task arrives
                if (batch.isEmpty()) {
                    previous = System.currentTimeMillis()
                }
                batch.add(pending)
            }

            // Check if we should flush
            final now = System.currentTimeMillis()
            final delta = now - previous

            if (!batch.isEmpty()) {
                // Flush if: time elapsed OR batch full OR shutting down
                if (delta > period || batch.size() >= TASKS_PER_REQUEST || completed) {
                    flushBatch(batch)
                    previous = System.currentTimeMillis()
                    batch.clear()
                }
            }
        }

        // Final flush of any remaining tasks
        if (!batch.isEmpty()) {
            flushBatch(batch)
        }
    }

    /**
     * Submit a batch of tasks to the scheduler API
     */
    protected void flushBatch(List<PendingTask> batch) {
        log.debug "[SEQERA] Submitting batch of ${batch.size()} tasks"

        try {
            // Extract Task objects for API call
            final List<Task> tasks = batch.collect { it.task }

            // Submit batch to API
            final response = client.createTasks(tasks)
            final List<String> taskIds = response.getTaskIds()

            // Validate response
            if (taskIds.size() != batch.size()) {
                throw new IllegalStateException("Seqera Scheduler API returned ${taskIds.size()} task IDs but submitted ${batch.size()} tasks")
            }

            // Map task IDs back to handlers
            for (int i = 0; i < batch.size(); i++) {
                final handler = batch[i].handler
                final taskId = taskIds[i]
                handler.setBatchTaskId(taskId)
            }

            log.debug "[SEQERA] Batch submission complete: ${taskIds.size()} tasks submitted"

        } catch (Exception e) {
            log.error "[SEQERA] Batch submission failed for ${batch.size()} tasks", e

            // Propagate failure to all handlers in this batch
            for (PendingTask pending : batch) {
                try {
                    pending.handler.onBatchSubmitFailure(e)
                } catch (Exception ex) {
                    log.warn "[SEQERA] Error handling batch failure for task", ex
                }
            }
        }
    }
}
