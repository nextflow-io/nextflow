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

import java.util.concurrent.CompletableFuture
import java.util.concurrent.ExecutorService
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException

import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.InputFilesMetrics
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.client.SchedClient
import nextflow.SysEnv
import nextflow.util.Duration
import nextflow.util.ThreadPoolBuilder
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

    /** Keep-alive interval - send empty submission to maintain session */
    static final Duration KEEP_ALIVE_INTERVAL = SysEnv.get('NXF_SEQERA_KEEP_ALIVE_INTERVAL', '60 sec') as Duration

    /** Timeout for waiting on metrics computation */
    static final Duration METRICS_TIMEOUT = SysEnv.get('NXF_SEQERA_METRICS_TIMEOUT', '30 sec') as Duration

    /**
     * Holds a task handler, its prepared Task object, and async metrics computation
     */
    @TupleConstructor
    static class PendingTask {
        SeqeraTaskHandler handler
        Task task
        CompletableFuture<InputFilesMetrics> metricsFuture
    }

    private final SchedClient client
    private final String sessionId
    private final Duration requestInterval
    private final Duration keepAliveInterval
    private final Closure onError
    private final LinkedBlockingQueue<PendingTask> pendingQueue = new LinkedBlockingQueue<>()
    private Thread sender
    private volatile boolean completed = false

    /** Executor pool for async input file metrics computation */
    private final ExecutorService metricsExecutor

    SeqeraBatchSubmitter(SchedClient client, String sessionId) {
        this(client, sessionId, REQUEST_INTERVAL, KEEP_ALIVE_INTERVAL)
    }

    SeqeraBatchSubmitter(SchedClient client, String sessionId, Duration requestInterval) {
        this(client, sessionId, requestInterval, KEEP_ALIVE_INTERVAL)
    }

    SeqeraBatchSubmitter(SchedClient client, String sessionId, Duration requestInterval, Duration keepAliveInterval, Closure onError=null) {
        this.client = client
        this.sessionId = sessionId
        this.requestInterval = requestInterval
        this.keepAliveInterval = keepAliveInterval
        this.onError = onError
        // Create a thread pool for metrics computation
        this.metricsExecutor = new ThreadPoolBuilder()
                .withName('seqera-metrics')
                .withMinSize(0)
                .withMaxSize(10)
                .withKeepAliveTime(60_000L)
                .withAllowCoreThreadTimeout(true)
                .build()
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
     * Starts async computation of input files metrics immediately.
     */
    void enqueue(SeqeraTaskHandler handler, Task task) {
        if (completed) {
            throw new IllegalStateException("Batch submitter has been shutdown")
        }

        // Start async metrics computation
        final taskRun = handler.task
        final metricsFuture = CompletableFuture.supplyAsync(
            ()-> InputFilesComputer.compute(taskRun),
            metricsExecutor
        )

        pendingQueue.add(new PendingTask(handler, task, metricsFuture))
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
        // Shutdown metrics executor
        metricsExecutor.shutdown()
        try {
            if (!metricsExecutor.awaitTermination(30, TimeUnit.SECONDS)) {
                metricsExecutor.shutdownNow()
            }
        }
        catch (InterruptedException e) {
            metricsExecutor.shutdownNow()
            Thread.currentThread().interrupt()
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

        try {
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
                else if (delta > keepAliveInterval.millis) {
                    // Keep-alive: send empty submission to maintain session
                    try {
                        log.debug "[SEQERA] Sending keep-alive for session ${sessionId}"
                        client.createTasks(sessionId, Collections.emptyList())
                    }
                    catch (Exception e) {
                        log.warn "[SEQERA] Keep-alive failed: ${e.message}"
                        // Don't crash the thread for keep-alive failures
                    }
                    // Always update timestamp to avoid rapid retry on failure
                    previous = System.currentTimeMillis()
                }
            }

            // Final flush of any remaining tasks
            if (!batch.isEmpty()) {
                flushBatch(batch)
            }
        }
        catch (Throwable e) {
            log.error "[SEQERA] Fatal error in batch submitter thread", e
            // Convert Throwable to Exception for handler API
            final Exception exception = e instanceof Exception ? (Exception) e : new RuntimeException(e)
            // Fail any tasks in the current batch
            for (PendingTask pending : batch) {
                try {
                    pending.handler.onBatchSubmitFailure(exception)
                }
                catch (Exception ex) {
                    log.warn "[SEQERA] Error failing batch task", ex
                }
            }
            // Drain and fail any remaining pending tasks
            drainAndFailPendingTasks(exception)
            // Invoke error callback to abort session
            if (onError) {
                try {
                    onError.call(e)
                }
                catch (Throwable t) {
                    log.warn "[SEQERA] Error in failure callback", t
                }
            }
        }
    }

    /**
     * Drain the pending queue and fail all tasks with the given error
     */
    private void drainAndFailPendingTasks(Exception cause) {
        PendingTask pending
        while ((pending = pendingQueue.poll()) != null) {
            try {
                pending.handler.onBatchSubmitFailure(cause)
            }
            catch (Exception e) {
                log.warn "[SEQERA] Error failing pending task", e
            }
        }
    }

    /**
     * Submit a batch of tasks to the scheduler API
     */
    protected void flushBatch(List<PendingTask> batch) {
        log.debug "[SEQERA] Submitting batch of ${batch.size()} tasks"

        try {
            // Resolve async metrics for all tasks in batch
            resolveMetrics(batch)

            // Extract Task objects for API call
            final List<Task> tasks = batch.collect { it.task }

            // Submit batch to API
            final response = client.createTasks(sessionId, tasks)
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

    /**
     * Wait for and attach metrics to all tasks in the batch.
     * Uses timeout to avoid blocking indefinitely on slow computations.
     */
    private void resolveMetrics(List<PendingTask> batch) {
        final timeout = METRICS_TIMEOUT.millis

        for (PendingTask pending : batch) {
            try {
                final metrics = pending.metricsFuture.get(timeout, TimeUnit.MILLISECONDS)
                if (metrics) {
                    pending.task.inputFilesMetrics(metrics)
                }
            }
            catch (TimeoutException e) {
                log.warn "[SEQERA] Timeout computing input files metrics for task: ${pending.handler.task.name}"
                pending.metricsFuture.cancel(true)
            }
            catch (Exception e) {
                log.warn "[SEQERA] Failed to compute input files metrics for task: ${pending.handler.task.name} - ${e.message}"
            }
        }
    }
}
