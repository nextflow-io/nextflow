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

import java.net.http.HttpRequest
import java.net.http.HttpResponse

/**
 * Class with common API calls used in different classes
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
class TowerCommonApi {

    /**
     * Calls the Seqera Platform user-info API to retrieve user information.
     *
     * @param client HTTP client to perform the API calls
     * @param endpoint Seqera Platform API endpoint
     * @return Map containing user information (id, userName, email, etc.)
     * @throws RuntimeException if the API call fails
     */
    Map getUserInfo(HxClient client, String endpoint) {
        final json = apiGet(client, endpoint, "/user-info")
        return json.user as Map
    }

    /**
     * Calls the Seqera Platform to retrieve a the user's workspaces information
     * and select the one matching with the workspace Id.
     *
     * @param client HTTP client to perform the API calls
     * @param userId Id of the workspace user
     * @param endpoint Seqera Platform API endpoint
     * @param workspaceId Id of the workspace
     * @return Map containing workspace informatiion
     * @throws RuntimeException if the API call fails
     */
    Map getUserWorkspaceDetails(HxClient client, String userId, String endpoint, String workspaceId) {
        if( !userId || !workspaceId ) {
            return null
        }
        try {
            final json = apiGet(client, endpoint, "/user/${userId}/workspaces")

            final orgsAndWorkspaces = json.orgsAndWorkspaces as List

            final workspace = orgsAndWorkspaces.find { ((Map) it).workspaceId?.toString() == workspaceId }
            if( workspace ) {
                final ws = workspace as Map
                return [
                    orgName          : ws.orgName,
                    workspaceId      : ws.workspaceId,
                    workspaceName    : ws.workspaceName,
                    workspaceFullName: ws.workspaceFullName,
                    roles            : ws.roles
                ]
            }

            return null
        } catch( Exception e ) {
            log.debug("Failed to get workspace details for workspace ${workspaceId}: ${e.message}", e)
            return null
        }
    }

    /**
     * Calls the Seqera Platform to retrieve a the workflow information.
     *
     * @param client HTTP client to perform the API calls
     * @param endpoint Seqera Platform API endpoint
     * @param workflowId Id of the workflow
     * @return Map containing workflow information
     * @throws RuntimeException if the API call fails
     */
    Map getWorkflowDetails(HxClient client, String endpoint, String workflowId, Map queryParams = [:]) {
        final json = apiGet(client, endpoint, "/workflow/${workflowId}", queryParams)
        return json.workflow as Map
    }

    Map apiGet(HxClient client, String apiEndpoint, String path, Map queryParams = [:]) {
        final url = buildUrl(apiEndpoint, path, queryParams)
        log.debug "Platform API - GET ${url}"
        final request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("API GET request ${url} failed: ${error}")
        }

        return new JsonSlurper().parseText(response.body()) as Map
    }

    String buildUrl(String endpoint, String path, Map queryParams) {
        def url = new StringBuilder(endpoint)
        if( !path.startsWith('/') ) {
            url.append('/')
        }
        url.append(path)

        if( queryParams && !queryParams.isEmpty() ) {
            url.append('?')
            url.append(queryParams.collect { k, v -> "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}" }.join('&'))
        }

        return url.toString()
    }
}
