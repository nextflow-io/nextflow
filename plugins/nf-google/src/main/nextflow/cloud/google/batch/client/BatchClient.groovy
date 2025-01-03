/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.client

import java.time.temporal.ChronoUnit
import java.util.concurrent.TimeoutException
import java.util.function.Predicate

import com.google.api.gax.core.CredentialsProvider
import com.google.api.gax.retrying.RetrySettings
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.api.gax.rpc.NotFoundException
import com.google.api.gax.rpc.UnavailableException
import com.google.auth.Credentials
import com.google.cloud.batch.v1.BatchServiceClient
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.JobName
import com.google.cloud.batch.v1.JobStatus
import com.google.cloud.batch.v1.LocationName
import com.google.cloud.batch.v1.Task
import com.google.cloud.batch.v1.TaskGroupName
import com.google.cloud.batch.v1.TaskName
import com.google.cloud.batch.v1.TaskStatus
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.threeten.bp.Duration
/**
 * Implements Google Batch HTTP client
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchClient {

    protected String projectId
    protected String location
    protected BatchServiceClient batchServiceClient
    protected BatchConfig config

    BatchClient(BatchConfig config) {
        this.config = config
        this.projectId = config.googleOpts.projectId
        this.location = config.googleOpts.location
        this.batchServiceClient = createBatchService(config)
    }

    /** Only for testing - do not use */
    protected BatchClient() {}

    protected CredentialsProvider createCredentialsProvider(BatchConfig config) {
        if( !config.getCredentials() )
            return null
        return new CredentialsProvider() {
            @Override
            Credentials getCredentials() throws IOException {
                return config.getCredentials()
            }
        }
    }

    protected BatchServiceSettings.Builder applyRetrySettings(BatchConfig batchConfig, BatchServiceSettings.Builder settings) {
        log.info "[GOOGLE BATCH] Applying custom retry configuration"
        log.info "[GOOGLE BATCH] google.grpc.initialRetryDelay = " + batchConfig.grpcInitialRetryDelay
        log.info "[GOOGLE BATCH] google.grpc.retryDelayMultiplier = " + batchConfig.grpcRetryDelayMultiplier
        log.info "[GOOGLE BATCH] google.grpc.maxRetryDelay = " + batchConfig.grpcMaxRetryDelay
        log.info "[GOOGLE BATCH] google.grpc.initialRpcTimeout = " + batchConfig.grpcInitialRpcTimeout
        log.info "[GOOGLE BATCH] google.grpc.timeoutMultiplier = " + batchConfig.grpcTimeoutMultiplier
        log.info "[GOOGLE BATCH] google.grpc.maxRpcTimeout = " + batchConfig.grpcMaxRpcTimeout
        log.info "[GOOGLE BATCH] google.grpc.totalTimeout = " + batchConfig.grpcTotalTimeout
        log.info "[GOOGLE BATCH] google.grpc.retryableCodes = " + batchConfig.grpcRetryableCodes
        final retrySettings = RetrySettings.newBuilder()
                .setInitialRetryDelay(Duration.ofMillis(batchConfig.grpcInitialRetryDelay.toMillis()))
                .setRetryDelayMultiplier(batchConfig.grpcRetryDelayMultiplier)
                .setMaxRetryDelay(Duration.ofMillis(batchConfig.grpcMaxRetryDelay.toMillis()))
                .setInitialRpcTimeout(Duration.ofMillis(batchConfig.grpcInitialRpcTimeout.toMillis()))
                .setRpcTimeoutMultiplier(batchConfig.grpcTimeoutMultiplier)
                .setMaxRpcTimeout(Duration.ofMillis(batchConfig.grpcMaxRpcTimeout.toMillis()))
                .setTotalTimeout(Duration.ofMillis(batchConfig.grpcTotalTimeout.toMillis()))
                .build()
        settings.getJobSettings()
                .setRetryableCodes(batchConfig.getGrpcRetryableCodes())
                .setRetrySettings(retrySettings)
        settings.listJobsSettings()
                .setRetryableCodes(batchConfig.getGrpcRetryableCodes())
                .setRetrySettings(retrySettings)
        settings.getTaskSettings()
                .setRetryableCodes(batchConfig.getGrpcRetryableCodes())
                .setRetrySettings(retrySettings)
        settings.listTasksSettings()
                .setRetryableCodes(batchConfig.getGrpcRetryableCodes())
                .setRetrySettings(retrySettings)
        return settings
    }

    protected BatchServiceClient createBatchService(BatchConfig config) {
        final provider = createCredentialsProvider(config)
        final settings = applyRetrySettings(config, BatchServiceSettings.newBuilder())
        if( provider ) {
            log.debug "[GOOGLE BATCH] Creating service client with config credentials"
            final userAgent = FixedHeaderProvider.create('user-agent', 'Nextflow')
            settings.setHeaderProvider(userAgent)
            settings.setCredentialsProvider(provider)
        }

        return BatchServiceClient.create(settings.build())
    }

    Job submitJob(String jobId, Job job) {
        final parent = LocationName.of(projectId, location)
        return apply(()-> batchServiceClient.createJob(parent, job, jobId))
    }

    Job describeJob(String jobId) {
        final name = JobName.of(projectId, location, jobId)
        return apply(()-> batchServiceClient.getJob(name))
    }

    Iterable<Task> listTasks(String jobId) {
        final parent = TaskGroupName.of(projectId, location, jobId, 'group0')
        return apply(()-> batchServiceClient.listTasks(parent).iterateAll())
    }

    Task describeTask(String jobId, String taskId) {
        final name = TaskName.of(projectId, location, jobId, 'group0', taskId)
        return apply(()-> batchServiceClient.getTask(name))
    }

    void deleteJob(String jobId) {
        final name = JobName.of(projectId, location, jobId).toString()
        apply(()-> batchServiceClient.deleteJobAsync(name))
    }

    TaskStatus getTaskStatus(String jobId, String taskId) {
        return describeTask(jobId, taskId).getStatus()
    }

    JobStatus getJobStatus(String jobId) {
        return describeJob(jobId).getStatus()
    }

    String getTaskState(String jobId, String taskId) {
        final status = getTaskStatus(jobId, taskId)
        return status ? status.getState().toString() : null
    }

    void shutdown() {
        batchServiceClient.close()
    }

    String getLocation() {
        return location
    }

    /**
     * Creates a retry policy using the configuration specified by {@link BatchRetryConfig}
     *
     * @param cond A predicate that determines when a retry should be triggered
     * @return The {@link dev.failsafe.RetryPolicy} instance
     */
    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond) {
        final cfg = config.getRetryConfig()
        final listener = new EventListener<ExecutionAttemptedEvent<T>>() {
            @Override
            void accept(ExecutionAttemptedEvent<T> event) throws Throwable {
                log.debug("[GOOGLE BATCH] response error - attempt: ${event.attemptCount}; reason: ${event.lastFailure.message}")
            }
        }
        return RetryPolicy.<T>builder()
            .handleIf(cond)
            .withBackoff(cfg.delay.toMillis(), cfg.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(cfg.maxAttempts)
            .withJitter(cfg.jitter)
            .onRetry(listener)
            .build()
    }

    /**
     * Carry out the invocation of the specified action using a retry policy
     * when an API UnavailableException is thrown
     *
     * see also https://github.com/nextflow-io/nextflow/issues/4537
     *
     * @param action A {@link dev.failsafe.function.CheckedSupplier} instance modeling the action to be performed in a safe manner
     * @return The result of the supplied action
     */
    protected <T> T apply(CheckedSupplier<T> action) {
        // define the retry condition
        final cond = new Predicate<? extends Throwable>() {
            @Override
            boolean test(Throwable t) {
                if( t instanceof UnavailableException )
                    return true
                if( t instanceof IOException || t.cause instanceof IOException )
                    return true
                if( t instanceof TimeoutException || t.cause instanceof TimeoutException )
                    return true
                if( t instanceof NotFoundException || t.cause instanceof NotFoundException )
                    return true
                return false
            }
        }
        // create the retry policy object
        final policy = retryPolicy(cond)
        // apply the action with
        return Failsafe.with(policy).get(action)
    }
}
