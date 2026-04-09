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

package io.seqera.tower.plugin

import groovy.json.JsonSlurper
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.SysEnv
import nextflow.config.ConfigBuilder
import nextflow.util.Duration

@Slf4j
class BaseCommandImpl {

    protected static final int API_TIMEOUT_MS = 10_000

    /**
     * Creates a TowerClient instance with optional authentication token.
     *
     * @param apiUrl Seqera Platform API url
     * @param accessToken Optional personal access token for authentication (PAT)
     * @return Configured TowerClient instance with timeout settings
     */
    @Memoized
    protected TowerClient createTowerClient(String apiUrl, String accessToken) {
        final env = SysEnv.get()
        return new TowerClient( new TowerConfig( [accessToken: accessToken, endpoint: apiUrl, httpConnectTimeout: Duration.of(API_TIMEOUT_MS)], env), env)
    }

    /**
     * Convert API endpoint to web URL
     * e.g., https://api.cloud.seqera.io -> https://cloud.seqera.io
     *      https://cloud.seqera.io/api -> https://cloud.seqera.io
     */
    protected String getWebUrlFromApiEndpoint(String apiEndpoint) {
        return apiEndpoint.replace('://api.', '://').replace('/api', '')
    }

    protected Map readConfig() {
        final builder = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
        return builder.buildConfigObject().flatten()
    }

    protected List<Map> listUserWorkspaces(TowerClient client, String userId) {
        return client.listUserWorkspacesAndOrgs(userId).findAll { ((Map) it).workspaceId != null }
    }

    protected List listComputeEnvironments(TowerClient client, String workspaceId) {
        try {
            final json = client.apiGet("/compute-envs", workspaceId ? [workspaceId: workspaceId] : [:])
            return json.computeEnvs as List ?: []
        } catch ( Exception e ) {
             throw new RuntimeException("Failed to get compute environments: ${e.message}", e)
        }
    }

    protected Map getComputeEnvironment(TowerClient client, String computeEnvId, String workspaceId) {
        try {
            final json = client.apiGet(workspaceId ? "/compute-envs/${computeEnvId}" : "/compute-envs",  workspaceId ? [workspaceId: workspaceId] : [:])
            return unifyComputeEnvDescription(json.computeEnv as Map ?: [:])
        } catch ( Exception e ) {
             throw new RuntimeException("Failed to get compute environments: ${e.message}", e)
        }
    }

    private Map unifyComputeEnvDescription(Map computeEnv) {
        if (computeEnv && !computeEnv.workDir) {
            final config = computeEnv?.config as Map
            log.debug("Config $config")
            computeEnv.workDir = config?.workDir as String
        }
        return computeEnv
    }
}
