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

package nextflow.plugin

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.time.Duration
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import org.pf4j.update.SimpleFileDownloader
/**
 * FileDownloader extension that enables the download of OCI compliant artifact that require a token authorization.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class OciAwareFileDownloader extends SimpleFileDownloader {

    private static final Pattern WWW_AUTH_PATTERN = ~/Bearer realm="([^"]+)",\s*service="([^"]+)",\s*scope="([^"]+)"/
    private static final Set<Integer> REDIRECT_STATUS_CODES = [301, 302, 303, 307, 308] as Set
    private static final Duration REQUEST_TIMEOUT = Duration.ofSeconds(90)

    /**
     * OCI aware download with token authorization. Tries to download the artifact and if it fails checks the headers to get the
     * @param fileUrl source file
     * @return
     */
    private final HttpClient httpClient = HttpClient.newBuilder()
        .followRedirects(HttpClient.Redirect.NEVER)
        .connectTimeout(Duration.ofSeconds(30))
        .build()

    @Override
    protected Path downloadFileHttp(URL fileUrl) {
        final destination = Files.createTempFile("nf-plugin", ".zip")

        // sendRequest now handles redirects and authentication internally
        HttpResponse<InputStream> response = sendRequest(fileUrl)

        // Check for HTTP error status codes
        if (response.statusCode() >= 400) {
            closeResponse(response)
            throw new IOException("HTTP error ${response.statusCode()} downloading from $fileUrl")
        }

        // Save the byte stream response directly to file
        return downloadFileFromResponse(response, destination)
    }

    private HttpRequest.Builder createRequestBuilder(URL url) {
        return HttpRequest.newBuilder()
            .uri(url.toURI())
            .header("User-Agent", "Nextflow/$BuildInfo.version")
            .header("X-Nextflow-Version", BuildInfo.version)
            .timeout(REQUEST_TIMEOUT)
            .GET()
    }
    
    private HttpResponse<InputStream> sendRequest0(URL url, String token = null) {
        final requestBuilder = createRequestBuilder(url)
        if (token) {
            requestBuilder.header("Authorization", "Bearer $token")
        }
        final HttpRequest request = requestBuilder.build()
        final HttpResponse<InputStream> response = httpClient.send(request, HttpResponse.BodyHandlers.ofInputStream())
        log.debug "HTTP response [${response.statusCode()}] from ${url}"
        return response
    }
    
    private HttpResponse<InputStream> sendRequest(URL url) {
        Set<String> attemptedUrls = new HashSet<>()
        URL currentUrl = url
        String token = null

        while (true) {
            // submit the request
            HttpResponse<InputStream> response = sendRequest0(currentUrl, token)
            
            // Handle redirects
            if (response.statusCode() in REDIRECT_STATUS_CODES) {
                // Prevent infinite redirect loops
                if (attemptedUrls.contains(currentUrl.toString())) {
                    closeResponse(response)
                    throw new IOException("Redirect loop detected for URL: $currentUrl")
                }
                attemptedUrls.add(currentUrl.toString())
                
                final newUrl = response.headers().firstValue("Location").orElse(null)
                if (newUrl) {
                    log.trace "Following redirect from $currentUrl to $newUrl"
                    closeResponse(response)
                    currentUrl = URI.create(newUrl).toURL()
                    token = null
                    continue
                }
            }
            // Handle authentication - retry once with token
            if (response.statusCode() == 401 && token == null) {
                token = handleAuthentication(response)
                log.trace "Retrying request with authentication token"
                closeResponse(response)
                continue
            }
            
            return response
        }
    }
    
    private String handleAuthentication(HttpResponse<?> response) {
        log.debug("Received 401 - attempting OCI token auth")
        def wwwAuth = response.headers().firstValue("WWW-Authenticate").orElse(null)
        if (!wwwAuth?.contains("Bearer")) {
            throw new IOException("Unsupported authentication method")
        }
        
        def matcher = WWW_AUTH_PATTERN.matcher(wwwAuth)
        if (!matcher.find()) {
            throw new IOException("Invalid WWW-Authenticate header: $wwwAuth")
        }
        
        def (realm, service, scope) = [matcher.group(1), matcher.group(2), matcher.group(3)]
        def tokenUrl = "${realm}?service=${URLEncoder.encode(service, 'UTF-8')}&scope=${URLEncoder.encode(scope, 'UTF-8')}"
        return fetchToken(tokenUrl)
    }
    
    private Path downloadFileFromResponse(HttpResponse<InputStream> response, Path destination) {
        log.debug "Saving downloaded file to: $destination"
        
        try (InputStream inputStream = response.body()) {
            Files.copy(inputStream, destination, StandardCopyOption.REPLACE_EXISTING)
        }
        finally {
            closeResponse(response)
        }
        return destination
    }
    

    private String fetchToken(String tokenUrl) {
        HttpRequest request = HttpRequest.newBuilder()
            .uri(URI.create(tokenUrl))
            .header("Accept", "application/json")
            .header("User-Agent", "Nextflow/$BuildInfo.version")
            .header("X-Nextflow-Version", BuildInfo.version)
            .timeout(REQUEST_TIMEOUT)
            .GET()
            .build()
        HttpResponse<String> response = httpClient.send(request, HttpResponse.BodyHandlers.ofString())
        log.debug "HTTP token response from $tokenUrl: status=${response.statusCode()}"
        
        if (response.statusCode() >= 400) {
            throw new IOException("HTTP error ${response.statusCode()} fetching token from $tokenUrl")
        }
        
        def json = response.body()
        def matcher = json =~ /"token"\s*:\s*"([^"]+)"/
        if (matcher.find()) {
            return matcher.group(1)
        }
        throw new IOException("Token not found in response: $json")
    }

    static void closeResponse(HttpResponse<?> response) {
        log.trace "Closing HttpClient response: $response"
        try {
            // close the httpclient response to prevent leaks
            // https://bugs.openjdk.org/browse/JDK-8308364
            final b0 = response.body()
            if( b0 instanceof Closeable )
                b0.close()
        }
        catch (Throwable e) {
            log.debug "Unexpected error while closing http response - cause: ${e.message}", e
        }
    }
}


