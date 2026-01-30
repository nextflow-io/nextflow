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

package nextflow.module

import com.google.gson.Gson
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import io.seqera.npr.api.schema.v1.Module
import io.seqera.npr.api.schema.v1.ModuleRelease
import io.seqera.npr.api.schema.v1.PublishModuleResponse
import io.seqera.npr.api.schema.v1.SearchModulesResponse
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.serde.gson.GsonEncoder
import nextflow.util.RetryConfig

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path

/**
 * REST API client for Nextflow module registry using npr-api models
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class ModuleRegistryClient {

    private final RegistryConfig config
    private final HxClient httpClient

    ModuleRegistryClient(RegistryConfig config) {
        this.config = config ?: new RegistryConfig()
        this.httpClient = HxClient.newBuilder()
                .retryConfig(RetryConfig.config())
                .followRedirects(HttpClient.Redirect.NORMAL)
                .build()
    }

    private String encodeName(String name) {
        return URLEncoder.encode(
            name.startsWith('@') ? name.substring(1) : name,
            'UTF-8'
        )
    }

    /**
     * Fetch module metadata from the registry
     *
     * @param name The module name (e.g., "nf-core/fastqc")
     * @return Module object with metadata
     */
    Module fetchModule(String name) {
        def registryUrls = config.allUrls

        Exception lastError = null
        for (String registryUrl : registryUrls) {
            log.debug "Trying to fetch from $registryUrl"
            try {
                return fetchModuleFromRegistry(registryUrl, name)
            } catch (Exception e) {
                log.debug "Failed to fetch module from ${registryUrl}: ${e.message}"
                lastError = e
            }
        }

        throw new AbortOperationException(
            "Unable to fetch module ${name} from any configured registry",
            lastError
        )
    }

    /**
     * Fetch module from a specific registry URL
     */
    private Module fetchModuleFromRegistry(String registryUrl, String name) {
        def endpoint = "${registryUrl}/api/modules/${encodeName(name)}"
        def uri = URI.create(endpoint)

        def requestBuilder = HttpRequest.newBuilder()
            .uri(uri)
            .GET()

        // Add authentication if available
        log.debug "Getting auth from: ${registryUrl}"
        def token = config.getAuthTokenResolved(registryUrl)
        if (token) {
            requestBuilder.header("Authorization", "Bearer ${token}")
        }
        log.debug "Building request: ${registryUrl}"
        def request = requestBuilder.build()

        try {
            log.debug "Fetching module from: ${uri}"
            def response = httpClient.send(request, HttpResponse.BodyHandlers.ofString())
            def body = response.body()

            log.debug "Registry request: ${response.uri()}\n- code: ${response.statusCode()}\n- body: ${body}"

            if (response.statusCode() == 404) {
                throw new AbortOperationException("Module not found: ${name}")
            }

            if (response.statusCode() != 200) {
                throw new AbortOperationException(
                    "Invalid response from registry: ${uri}\n" +
                    "- http status: ${response.statusCode()}\n" +
                    "- response: ${body}"
                )
            }

            // Parse response using npr-api Module model
            def encoder = new GsonEncoder<Module>() {}
            return encoder.decode(body)
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            e.printStackTrace()
            throw new AbortOperationException("Failed to fetch module from: ${uri}", e)
        }
    }

    /**
     * Fetch specific module version/release
     *
     * @param name The module name
     * @param version The version string
     * @return ModuleRelease object from npr-api
     */
    ModuleRelease fetchRelease(String name, String version) {
        def registryUrls = config.allUrls

        Exception lastError = null
        for (String registryUrl : registryUrls) {
            try {
                return fetchReleaseFromRegistry(registryUrl, name, version)
            } catch (Exception e) {
                log.debug "Failed to fetch release from ${registryUrl}: ${e.message}"
                lastError = e
            }
        }

        throw new AbortOperationException(
            "Unable to fetch module ${name}@${version} from any configured registry",
            lastError
        )
    }

    /**
     * Fetch release from a specific registry URL
     */
    private ModuleRelease fetchReleaseFromRegistry(String registryUrl, String name, String version) {
        def endpoint = "${registryUrl}/api/modules/${encodeName(name)}/${version}"
        def uri = URI.create(endpoint)

        def requestBuilder = HttpRequest.newBuilder()
            .uri(uri)
            .GET()

        def token = config.getAuthTokenResolved(registryUrl)
        if (token) {
            requestBuilder.header("Authorization", "Bearer ${token}")
        }

        def request = requestBuilder.build()

        try {
            log.debug "Fetching module release from: ${uri}"
            def response = httpClient.send(request, HttpResponse.BodyHandlers.ofString())
            def body = response.body()

            if (response.statusCode() == 404) {
                throw new AbortOperationException("Module version not found: ${name}@${version}")
            }

            if (response.statusCode() != 200) {
                throw new AbortOperationException(
                    "Invalid response from registry: ${uri}\n" +
                    "- http status: ${response.statusCode()}\n" +
                    "- response: ${body}"
                )
            }

            // Parse response using npr-api ModuleRelease model
             return new GsonEncoder<ModuleRelease>() {}.decode(body)
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to fetch module release from: ${uri}", e)
        }
    }

    /**
     * Download a module bundle from the registry
     *
     * @param name The module name
     * @param version The module version
     * @param targetPath The target path to download to
     * @return Path with the downloaded file path
     */
    Path downloadModule(String name, String version, Path targetPath) {
        def registryUrls = config.allUrls
        if (targetPath.exists()){
            targetPath.delete()
        }
        Exception lastError = null
        for (String registryUrl : registryUrls) {
            try {
                return downloadModuleFromRegistry(registryUrl, name, version, targetPath)
            } catch (Exception e) {
                log.debug "Failed to download from ${registryUrl}: ${e.message}"
                lastError = e
            }
        }

        throw new AbortOperationException(
            "Unable to download module ${name}@${version} from any configured registry",
            lastError
        )
    }

    /**
     * Download module from a specific registry URL
     */
    private Path downloadModuleFromRegistry(String registryUrl, String name, String version, Path targetPath) {
        def endpoint = "${registryUrl}/api/modules/${encodeName(name)}/${version}/download"
        def uri = URI.create(endpoint)

        def requestBuilder = HttpRequest.newBuilder()
            .uri(uri)
            .GET()

        def token = config.getAuthTokenResolved(registryUrl)
        if (token) {
            requestBuilder.header("Authorization", "Bearer ${token}")
        }

        def request = requestBuilder.build()

        try {
            log.debug "Downloading module from: ${uri}"
            def response = httpClient.send(request, HttpResponse.BodyHandlers.ofInputStream())

            if (response.statusCode() == 404) {
                throw new AbortOperationException("Module bundle not found: ${name}@${version}")
            }

            if (response.statusCode() != 200) {
                throw new AbortOperationException(
                    "Invalid response from registry: ${uri}\n" +
                    "- http status: ${response.statusCode()}"
                )
            }

            // Create parent directories if needed
            if (targetPath.parent) {
                Files.createDirectories(targetPath.parent)
            }

            // Write response body to file
            Files.copy(response.body(), targetPath)
            log.debug "Downloaded module to: ${targetPath}"

            validateDownloadIntegrity(response, uri, targetPath, name, version)

            return targetPath
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to download module from: ${uri}", e)
        }
    }

    private void validateDownloadIntegrity(HttpResponse<InputStream> response, uri, Path targetPath, String name, String version) {
        // Get checksum from headers (X-Checksum or Docker-Content-Digest)
        def checksumType = ModuleChecksum.CHECKSUM_ALGORITHM
        def checksum = response.headers().firstValue("X-Checksum").orElse(null)

        if( !checksum ) {
            checksum = response.headers().firstValue("Docker-Content-Digest").orElse(null)
        }

        if( !checksum ) {
            log.warn "No X-Checksum or Docker-Content-Digest header found in response from ${uri}"
            return
        }

        // Check if checksum has a digest format including algorithm: "sha256:abc123..."
        def parts = checksum.split(':', 2)
        if( parts.length == 2 ) {
            checksumType = parts[0].toLowerCase()
            checksum = parts[1]
        }
        log.debug "Using checksum: ${checksumType}:${checksum}"

        def actualChecksum = ModuleChecksum.computeFile(targetPath, checksumType)
        if( actualChecksum != checksum ) {
            // Clean up downloaded file
            Files.delete(targetPath)
            throw new AbortOperationException(
                "Downloaded module checksum mismatch for ${name}@${version}:\n" +
                    "- expected (${checksumType}): ${checksum}\n" +
                    "- actual:                    ${actualChecksum}\n" +
                    "The download may be corrupted or tampered with."
            )
        }
        log.debug "Checksum validated successfully: ${checksumType}:${checksum}"
    }

    /**
     * Search for modules in the registry
     *
     * @param query The search query
     * @param limit Maximum number of results (default: 20)
     * @return SearchModulesResponse with results
     */
    SearchModulesResponse search(String query, int limit = 20) {
        def registryUrls = config.allUrls

        Exception lastError = null
        for (String registryUrl : registryUrls) {
            try {
                return searchInRegistry(registryUrl, query, limit)
            } catch (Exception e) {
                log.debug "Failed to search in ${registryUrl}: ${e.message}"
                lastError = e
            }
        }

        throw new AbortOperationException(
            "Unable to search modules in any configured registry",
            lastError
        )
    }

    /**
     * Search in a specific registry
     */
    private SearchModulesResponse searchInRegistry(String registryUrl, String query, int limit) {
        def endpoint = "${registryUrl}/api/modules?query=${URLEncoder.encode(query, 'UTF-8')}&limit=${limit}"
        def uri = URI.create(endpoint)

        def requestBuilder = HttpRequest.newBuilder()
            .uri(uri)
            .GET()

        def token = config.getAuthTokenResolved(registryUrl)
        if (token) {
            requestBuilder.header("Authorization", "Bearer ${token}")
        }

        def request = requestBuilder.build()

        try {
            log.debug "Searching modules: ${uri}"
            def response = httpClient.send(request, HttpResponse.BodyHandlers.ofString())
            def body = response.body()

            if (response.statusCode() != 200) {
                throw new AbortOperationException(
                    "Invalid response from registry: ${uri}\n" +
                    "- http status: ${response.statusCode()}\n" +
                    "- response: ${body}"
                )
            }

            // Parse response using npr-api SearchModulesResponse model
            def encoder = new GsonEncoder<SearchModulesResponse>() {}
            return encoder.decode(body)
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to search modules in: ${uri}", e)
        }
    }

    /**
     * Publish a module to the registry (authenticated)
     *
     * @param name The module name
     * @param request The publish request from npr-api
     * @param authToken The authentication token
     * @return PublishModuleResponse from npr-api
     */
    PublishModuleResponse publishModule(String name, def request, String registry = null) {
        final registryUrl = registry ?: config.url
        final authToken = config.getAuthTokenResolved(registryUrl)

        if (!authToken) {
            throw new AbortOperationException(
                "Authentication required to publish modules.\n" +
                    "Please set NXF_REGISTRY_TOKEN environment variable or configure registry.auth in nextflow.config:\n\n" +
                    "  registry {\n" +
                    "    auth {\n" +
                    "      '${registryUrl}' = '\${NXF_REGISTRY_TOKEN}'\n" +
                    "    }\n" +
                    "  }\n"
            )
        }
        try {
            return publishModuleToRegistry(registryUrl, name, request, authToken)
        } catch( Exception e ) {
            throw new AbortOperationException("Failed to publish to ${registryUrl}", e)
        }
    }

    /**
     * Publish module to a specific registry
     */
    private PublishModuleResponse publishModuleToRegistry(
            String registryUrl,
            String name,
            def request,
            String authToken) {

        String endpoint = "${registryUrl}/api/modules/${encodeName(name)}".toString()
        URI uri = URI.create(endpoint)

        // Serialize request to JSON
        def gson = new Gson()
        String requestBody = gson.toJson(request)

        HttpRequest httpRequest = HttpRequest.newBuilder()
            .uri(uri)
            .header("Content-Type", "application/json")
            .header("Authorization", "Bearer ${authToken}".toString())
            .POST(HttpRequest.BodyPublishers.ofString(requestBody))
            .build()

        try {
            log.debug "Publishing module to: ${uri}"
            log.trace "Request: \n\t${httpRequest}\n\theaders: ${httpRequest.headers()}\n\tbody: ${requestBody}"

            HttpResponse<String> response = httpClient.send(httpRequest, HttpResponse.BodyHandlers.ofString())
            String body = response.body()

            if (response.statusCode() != 201) {
                throw new AbortOperationException(
                    "Failed to publish module: ${uri}\n" +
                    "- http status: ${response.statusCode()}\n" +
                    "- response: ${body}"
                )
            }

            // Parse response using npr-api PublishModuleResponse model
            return new GsonEncoder<PublishModuleResponse>() {}.decode(body)
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to publish module to: ${uri}", e)
        }
    }
}
