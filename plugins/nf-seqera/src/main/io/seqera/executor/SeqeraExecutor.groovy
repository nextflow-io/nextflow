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

import groovy.util.logging.Slf4j
import io.seqera.client.SeqeraClient
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName
import org.pf4j.ExtensionPoint

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('seqera')
class SeqeraExecutor extends Executor implements ExtensionPoint {

    private SeqeraClient client

    private String clusterId

    @Override
    protected void register() {
        createClient()
        createCluster()
    }

    @Override
    public void shutdown() {
        deleteCluster()
    }

    protected void createCluster() {
        log.debug "[SEQERA] Creating cluster for workflow"
        final cluster = client.createCluster()
        this.clusterId = cluster.clusterId
        log.debug "[SEQERA] Cluster created id: " + cluster.clusterId
    }

    protected void deleteCluster() {
        log.debug "[SEQERA] Deleting cluster: " + clusterId
        client.deleteCluster(this.clusterId)
        log.debug "[SEQERA] Cluster id deleted"
    }

    protected void createClient() {
        this.client = new SeqeraClient(session)
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 1000, Duration.of('10 sec'))
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

    SeqeraClient getClient() {
        return client
    }

    String getClusterId() {
        return clusterId
    }
}
