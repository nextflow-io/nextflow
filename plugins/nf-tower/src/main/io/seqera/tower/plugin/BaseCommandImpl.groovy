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
class BaseCommandImpl {

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

    /**
     * Calls the Seqera Platform user-info API to retrieve user information.
     *
     * @param accessToken Authentication token
     * @param apiUrl Seqera Platform API endpoint
     * @return Map containing user information (id, userName, email, etc.)
     * @throws RuntimeException if the API call fails
     */
    protected Map getUserInfo(String accessToken, String apiUrl) {
        final client = createHttpClient(accessToken)
        final url = "${apiUrl}/user-info"
        log.debug "Platform get user info - GET ${url}"
        final request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if (response.statusCode() != 200) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get user info: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        return json.user as Map
    }

    protected Map getWorkspaceDetails(String accessToken, String endpoint, String workspaceId) {
        try {
            final userInfo = getUserInfo(accessToken, endpoint)
            final userId = userInfo.id as String

            final client = createHttpClient(accessToken)
            final url = "${endpoint}/user/${userId}/workspaces"
            log.debug "Platform get workdspace - GET ${url}"
            final request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build()

            final response = client.send(request, HttpResponse.BodyHandlers.ofString())

            if (response.statusCode() != 200) {
                return null
            }

            final json = new JsonSlurper().parseText(response.body()) as Map
            final orgsAndWorkspaces = json.orgsAndWorkspaces as List

            final workspace = orgsAndWorkspaces.find { ((Map)it).workspaceId?.toString() == workspaceId }
            if (workspace) {
                final ws = workspace as Map
                return [
                    orgName: ws.orgName,
                    workspaceName: ws.workspaceName,
                    workspaceFullName: ws.workspaceFullName,
                    roles: ws.roles
                ]
            }

            return null
        } catch (Exception e) {
            log.debug("Failed to get workspace details for workspace ${workspaceId}: ${e.message}", e)
            return null
        }
    }

    protected List listUserWorkspaces(String accessToken, String endpoint, String userId) {
        final client = createHttpClient(accessToken)
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

    protected List listComputeEnvironments(String accessToken, String endpoint, String workspaceId) {
        final client = createHttpClient(accessToken)
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

    protected Map getComputeEnvironment(String accessToken, String endpoint, String computeEnvId, String workspaceId) {
        final client = createHttpClient(accessToken)
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
