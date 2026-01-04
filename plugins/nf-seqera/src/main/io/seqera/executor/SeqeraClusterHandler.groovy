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
import io.seqera.sched.client.SchedClient

/**
 * Handles cluster management operations for Seqera executor
 */
@Slf4j
class SeqeraClusterHandler {

    private SchedClient client
    private String clusterId

    SeqeraClusterHandler(SchedClient client) {
        this.client = client
    }

    void createCluster() {
        log.debug "[SEQERA] Creating cluster for workflow"
        final cluster = client.createCluster()
        this.clusterId = cluster.getClusterId()
        log.debug "[SEQERA] Cluster created id: ${cluster.getClusterId()}"
    }

    void deleteCluster() {
        if (!clusterId) {
            return
        }
        log.debug "[SEQERA] Terminating cluster: ${clusterId}"
        client.terminateCluster(this.clusterId)
        log.debug "[SEQERA] Cluster terminated"
    }

    String getClusterId() {
        return clusterId
    }
}