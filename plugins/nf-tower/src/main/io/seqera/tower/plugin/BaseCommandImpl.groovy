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

package io.seqera.tower.plugin

import groovy.json.JsonSlurper
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import nextflow.Const
import nextflow.config.ConfigBuilder

import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration

@Slf4j
class BaseCommandImpl extends TowerCommonApi {

    private static final int API_TIMEOUT_MS = 10_000

    /**
     * Creates an HxClient instance with optional authentication token.
     *
     * @param accessToken Optional personal access token for authentication (PAT)
     * @return Configured HxClient instance with timeout settings
     */
    protected HxClient createHttpClient(String accessToken = null) {
        return HxClient.newBuilder()
            .connectTimeout(Duration.ofMillis(API_TIMEOUT_MS))
            .bearerToken(accessToken)
            .build()
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

    protected List listUserWorkspaces(HxClient client, String endpoint, String userId) {
        final url = "${endpoint}/user/${userId}/workspaces"
        log.debug "Platform list workspaces - GET ${url}"
        final request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get workspaces: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        final orgsAndWorkspaces = json.orgsAndWorkspaces as List

        return orgsAndWorkspaces.findAll { ((Map) it).workspaceId != null }
    }

    protected List listComputeEnvironments(HxClient client, String endpoint, String workspaceId) {
        final uri = workspaceId
            ? "${endpoint}/compute-envs?workspaceId=${workspaceId}"
            : "${endpoint}/compute-envs"
        log.debug "Platform list compute env - GET ${uri}"

        final request = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get compute environments: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        return json.computeEnvs as List ?: []
    }

    protected Map getComputeEnvironment(HxClient client, String endpoint, String computeEnvId, String workspaceId) {
        final uri = workspaceId ?
            "${endpoint}/compute-envs/${computeEnvId}?workspaceId=${workspaceId}" :
            "${endpoint}/compute-envs"
        log.debug "Platform get compute env - GET ${uri}"

        final request = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get compute environment: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        return unifyComputeEnvDescription(json.computeEnv as Map ?: [:])
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
