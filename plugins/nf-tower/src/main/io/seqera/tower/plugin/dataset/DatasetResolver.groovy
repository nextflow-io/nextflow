/*
 * Copyright 2013-2024, Seqera Labs
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

package io.seqera.tower.plugin.dataset

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.platform.PlatformHelper

/**
 * Resolves a Seqera Platform dataset reference to its backing cloud storage path.
 * <p>
 * Resolution chain:
 * 1. Dataset name → GET /datasets?workspaceId=X → DatasetDto.id
 * 2. Dataset id + version → GET /datasets/{id}/versions → DatasetVersionDto.url
 * 3. Cloud URL string → FileHelper.asPath() → concrete cloud Path (S3/GCS/Azure)
 *
 * @author Edmund Miller
 */
@Slf4j
@CompileStatic
class DatasetResolver {

    /**
     * Resolve a dataset name (and optional version) to the backing cloud storage Path.
     *
     * @param datasetName The dataset name as shown in Seqera Platform
     * @param version     The version number (null = latest)
     * @return A concrete cloud storage Path (e.g. S3Path, GcsPath)
     */
    static Path resolve(String datasetName, String version) {
        if (!datasetName)
            throw new IllegalArgumentException("Dataset name cannot be null or empty")

        final String endpoint = getEndpoint()
        final String accessToken = getAccessToken()
        final String workspaceId = getWorkspaceId()

        if (!accessToken)
            throw new AbortOperationException("Missing Seqera Platform access token -- set TOWER_ACCESS_TOKEN or tower.accessToken in config")

        final HttpClient httpClient = HttpClient.newHttpClient()

        // Step 1: Resolve dataset name → dataset ID
        final String datasetId = resolveDatasetId(httpClient, endpoint, accessToken, workspaceId, datasetName)

        // Step 2: Resolve dataset ID + version → cloud storage URL
        final String cloudUrl = resolveCloudUrl(httpClient, endpoint, accessToken, workspaceId, datasetId, version)

        log.debug "Dataset '{}' resolved to cloud URL: {}", datasetName, cloudUrl

        // Step 3: Convert cloud URL → Path via Nextflow's FileHelper
        return FileHelper.asPath(cloudUrl)
    }

    /**
     * Look up a dataset by name, return its ID.
     */
    static private String resolveDatasetId(HttpClient httpClient, String endpoint, String accessToken, String workspaceId, String datasetName) {
        String url = "${endpoint}/datasets"
        if (workspaceId) {
            url += "?workspaceId=${workspaceId}"
        }

        log.debug "Listing datasets from: {}", url

        final Map json = httpGet(httpClient, url, accessToken)
        final List<Map> datasets = json.datasets as List<Map>

        if (!datasets) {
            throw new AbortOperationException("No datasets found in workspace")
        }

        // Find dataset by name (case-sensitive match)
        final Map dataset = datasets.find { Map it -> it.name == datasetName }
        if (!dataset) {
            final String available = datasets.collect { Map it -> it.name }.join(', ')
            throw new AbortOperationException("Dataset '${datasetName}' not found. Available: ${available}")
        }

        return dataset.id as String
    }

    /**
     * Look up the dataset version's backing cloud URL.
     * If version is null, uses the latest version.
     */
    static private String resolveCloudUrl(HttpClient httpClient, String endpoint, String accessToken, String workspaceId, String datasetId, String version) {
        String url = "${endpoint}/datasets/${datasetId}/versions"
        if (workspaceId) {
            url += "?workspaceId=${workspaceId}"
        }

        log.debug "Listing dataset versions from: {}", url

        final Map json = httpGet(httpClient, url, accessToken)
        final List<Map> versions = json.versions as List<Map>

        if (!versions) {
            throw new AbortOperationException("No versions found for dataset ID: ${datasetId}")
        }

        Map targetVersion
        if (version) {
            // Find specific version
            targetVersion = versions.find { Map it -> it.version?.toString() == version }
            if (!targetVersion) {
                throw new AbortOperationException("Version '${version}' not found for dataset ID: ${datasetId}")
            }
        }
        else {
            // Use the latest version (highest version number)
            targetVersion = versions.max { Map it -> (it.version as Integer) ?: 0 } as Map
        }

        final String cloudUrl = targetVersion.url as String
        if (!cloudUrl) {
            throw new AbortOperationException("Dataset version has no backing storage URL -- dataset ID: ${datasetId}, version: ${targetVersion.version}")
        }

        return cloudUrl
    }

    /**
     * Execute a GET request against the Platform API and parse the JSON response.
     */
    static private Map httpGet(HttpClient httpClient, String url, String accessToken) {
        final request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header("Authorization", "Bearer ${accessToken}")
                .GET()
                .build()

        final HttpResponse<String> response = httpClient.send(request, HttpResponse.BodyHandlers.ofString())

        if (response.statusCode() < 200 || response.statusCode() >= 300) {
            if (response.statusCode() == 401 || response.statusCode() == 403) {
                throw new AbortOperationException("Access denied to Seqera Platform API -- check your access token")
            }
            throw new AbortOperationException("Seqera Platform API error -- HTTP ${response.statusCode()}: ${response.body()}")
        }

        return new JsonSlurper().parseText(response.body()) as Map
    }

    // -- config helpers --

    static private String getEndpoint() {
        final Map opts = getTowerOpts()
        return PlatformHelper.getEndpoint(opts, System.getenv())
    }

    static private String getAccessToken() {
        final Map opts = getTowerOpts()
        return PlatformHelper.getAccessToken(opts, System.getenv())
    }

    static private String getWorkspaceId() {
        final Map opts = getTowerOpts()
        return PlatformHelper.getWorkspaceId(opts, System.getenv())
    }

    static private Map getTowerOpts() {
        final session = Global.session as Session
        if (!session)
            throw new AbortOperationException("Nextflow session not initialized -- dataset:// URIs require an active session")

        return session.config?.navigate('tower') as Map ?: Collections.emptyMap()
    }
}
