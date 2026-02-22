/*
 * Copyright 2013-2025, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.config.SeqeraConfig
import io.seqera.config.ExecutorOpts
import io.seqera.util.MapperUtil
import io.seqera.sched.api.schema.v1a1.CreateRunRequest
import io.seqera.sched.api.schema.v1a1.PredictionModel
import io.seqera.sched.client.SchedClient
import io.seqera.sched.api.schema.v1a1.TerminateRunRequest
import io.seqera.sched.client.SchedClientConfig
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.fusion.FusionHelper
import nextflow.platform.PlatformHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.SysEnv
import nextflow.util.Duration
import nextflow.util.ServiceName
import org.pf4j.ExtensionPoint

/**
 * Nextflow executor that delegates task execution to the Seqera scheduler API.
 *
 * <p>This executor creates a run on the Seqera scheduler, submits tasks in batches
 * via {@link SeqeraBatchSubmitter}, and monitors their lifecycle through the scheduler API.
 * It requires Fusion file system to be enabled and all processes to specify a container image.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName(SEQERA)
@CompileStatic
class SeqeraExecutor extends Executor implements ExtensionPoint {

    public static final String SEQERA = 'seqera'

    private ExecutorOpts seqeraConfig

    private SchedClient client

    private volatile String runId

    private SeqeraBatchSubmitter batchSubmitter

    @Override
    protected void register() {
        createClient()
    }

    @Override
    void shutdown() {
        // Flush any pending batch jobs before terminating run
        session.error
        batchSubmitter?.shutdown()
        terminateRun()
    }

    protected void createClient() {
        final seqera = new SeqeraConfig(session.config.seqera as Map ?: Collections.<String,Object>emptyMap())
        this.seqeraConfig = seqera.executor
        if (!seqeraConfig)
            throw new IllegalArgumentException("Missing Seqera executor configuration - make sure to specify 'seqera.executor' settings")
        // Get access token and refresh token from tower config (shares authentication with Platform)
        def towerConfig = session.config.tower as Map ?: Collections.emptyMap()
        def accessToken = PlatformHelper.getAccessToken(towerConfig, SysEnv.get())
        def refreshToken = PlatformHelper.getRefreshToken(towerConfig, SysEnv.get())
        def platformUrl = PlatformHelper.getEndpoint(towerConfig, SysEnv.get())
        def clientConfig = SchedClientConfig.builder()
                .endpoint(seqeraConfig.endpoint)
                .platformUrl(platformUrl)
                .accessToken(accessToken)
                .refreshToken(refreshToken)
                .retryConfig(seqeraConfig.retryOpts())
                .build()
        this.client = new SchedClient(clientConfig)
    }

    protected void createRun() {
        final towerConfig = session.config.tower as Map ?: Collections.emptyMap()
        final workflowId = session.workflowMetadata?.platform?.workflowId
        final workflowUrl = session.workflowMetadata?.platform?.workflowUrl
        final labels = new Labels()
        if( seqeraConfig.autoLabels )
            labels.withWorkflowMetadata(session.workflowMetadata)
        labels.withUserLabels(seqeraConfig.labels)
        final predictionModel = seqeraConfig.predictionModel ? PredictionModel.fromValue(seqeraConfig.predictionModel) : null
        final request = new CreateRunRequest()
                .region(seqeraConfig.region)
                .name(session.runName)
                .machineRequirement(MapperUtil.toMachineRequirement(seqeraConfig.machineRequirement))
                .labels(labels.entries)
                .workspaceId(PlatformHelper.getWorkspaceId(towerConfig, SysEnv.get()) as Long)
                .workflowId(workflowId)
                .workflowUrl(workflowUrl)
                .predictionModel(predictionModel)
        log.debug "[SEQERA] Creating run: ${request}"
        final response = client.createRun(request)
        this.runId = response.getRunId()
        log.debug "[SEQERA] Run created id: ${runId}; workflowId: '${workflowId}'; workflowUrl: '${workflowUrl}'"
        // Initialize and start batch submitter with error callback to abort on fatal errors
        this.batchSubmitter = new SeqeraBatchSubmitter(
            client,
            runId,
            seqeraConfig.batchFlushInterval,
            SeqeraBatchSubmitter.KEEP_ALIVE_INTERVAL,
            { Throwable t -> session.abort(t) }
        )
        this.batchSubmitter.start()
    }

    protected void terminateRun() {
        if (!runId) {
            return
        }
        final stopReason = truncate(session.fault?.report, 10_000)
        log.debug "[SEQERA] Terminating run: ${runId}; stopReason: ${stopReason}"
        client.terminateRun(runId, new TerminateRunRequest().stopReason(stopReason))
        log.debug "[SEQERA] Run terminated"
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, config, name, 1000, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new SeqeraTaskHandler(task, this)
    }

    /**
     * @return {@code true} whenever the containerization is managed by the executor itself
     */
    boolean isContainerNative() {
        return true
    }

    @Override
    boolean isFusionEnabled() {
        final enabled = FusionHelper.isFusionEnabled(session)
        if (!enabled)
            throw new AbortOperationException("Seqera executor requires the use of Fusion file system")
        return true
    }

    /**
     * Lazily creates the run on first access, ensuring workflowId and labels
     * are available (they are set by TowerClient.onFlowCreate before tasks are submitted).
     */
    void ensureRunCreated() {
        if (runId) return
        synchronized (this) {
            if (runId) return
            createRun()
        }
    }

    SchedClient getClient() {
        return client
    }

    String getRunId() {
        return runId
    }

    SeqeraBatchSubmitter getBatchSubmitter() {
        return batchSubmitter
    }

    ExecutorOpts getSeqeraConfig() {
        return seqeraConfig
    }

    protected static String truncate(String value, int maxLen) {
        if (!value || value.length() <= maxLen)
            return value
        return value.take(maxLen) + '\n.. [TRUNCATED]'
    }
}
