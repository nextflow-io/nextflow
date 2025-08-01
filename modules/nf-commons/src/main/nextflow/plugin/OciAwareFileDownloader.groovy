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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.util.retry.Retryable
import nextflow.util.HttpRetryableClient
import org.pf4j.update.SimpleFileDownloader

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path
import java.time.Duration
import java.util.regex.Pattern

/**
 * FileDownloader extension that enables the download of OCI compliant artifact that require a token authorization.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class OciAwareFileDownloader extends SimpleFileDownloader {

    private static final Pattern WWW_AUTH_PATTERN = ~/Bearer realm="([^"]+)",\s*service="([^"]+)",\s*scope="([^"]+)"/

    /**
     * OCI aware download with token authorization. Tries to download the artifact and if it fails checks the headers to get the
     * @param fileUrl source file
     * @return
     */
    @Override
    protected Path downloadFileHttp(URL fileUrl) {
        def retriableHttpClient = HttpRetryableClient.create(
            HttpClient.newBuilder()
                .followRedirects(HttpClient.Redirect.NEVER)
                .connectTimeout(Duration.ofSeconds(30))
                .build(),
            Retryable.ofDefaults().config()
        )


        Path destination = Files.createTempDirectory("pf4j-update-downloader")
        destination.toFile().deleteOnExit()

        String path = fileUrl.getPath()
        String fileName = path.substring(path.lastIndexOf('/') + 1)
        HttpRequest request = HttpRequest.newBuilder()
            .uri(fileUrl.toURI())
            .timeout(Duration.ofMinutes(5))
            .GET()
            .build()

        HttpResponse<String> response = retriableHttpClient.send(request, HttpResponse.BodyHandlers.ofString())

        // Handle redirects manually because of filename is got from path
        if (response.statusCode() in [301, 302, 303]) {
            String newUrl = response.headers().firstValue("Location").orElse(null)
            if (newUrl) {
                log.debug("Managing redirection to $newUrl")
                fileUrl = URI.create(newUrl).toURL()
                path = fileUrl.getPath()
                fileName = path.substring(path.lastIndexOf('/') + 1)
                request = HttpRequest.newBuilder()
                    .uri(fileUrl.toURI())
                    .timeout(Duration.ofMinutes(5))
                    .GET()
                    .build()
                response = retriableHttpClient.send(request, HttpResponse.BodyHandlers.ofString())
            }
        }

        if (response.statusCode() == 401) {
            def wwwAuth = response.headers().firstValue("WWW-Authenticate").orElse(null)
            if (wwwAuth?.contains("Bearer")) {
                log.debug("Received 401 — attempting OCI token auth")

                def matcher = WWW_AUTH_PATTERN.matcher(wwwAuth)
                if (!matcher.find()) {
                    throw new IOException("Invalid WWW-Authenticate header: $wwwAuth")
                }

                def (realm, service, scope) = [matcher.group(1), matcher.group(2), matcher.group(3)]
                def tokenUrl = "${realm}?service=${URLEncoder.encode(service, 'UTF-8')}&scope=${URLEncoder.encode(scope, 'UTF-8')}"
                def token = fetchToken(tokenUrl)

                // Retry download with Bearer token
                HttpRequest authRequest = HttpRequest.newBuilder()
                    .uri(fileUrl.toURI())
                    .header("Authorization", "Bearer $token")
                    .timeout(Duration.ofMinutes(5))
                    .GET()
                    .build()
                HttpResponse<byte[]> authResponse = retriableHttpClient.send(authRequest, HttpResponse.BodyHandlers.ofByteArray())
                Path file = destination.resolve(fileName)
                Files.write(file, authResponse.body())
                return file
            }
        }
        // Fallback to default behavior - download with initial response
        HttpResponse<byte[]> downloadResponse = retriableHttpClient.send(request, HttpResponse.BodyHandlers.ofByteArray())
        Path file = destination.resolve(fileName)
        Files.write(file, downloadResponse.body())
        return file
    }

    private String fetchToken(String tokenUrl) {
        def retriableHttpClient = HttpRetryableClient.create(
            HttpClient.newBuilder()
                .connectTimeout(Duration.ofSeconds(30))
                .build(),
            Retryable.ofDefaults().config()
        )

        HttpRequest request = HttpRequest.newBuilder()
            .uri(URI.create(tokenUrl))
            .header("Accept", "application/json")
            .timeout(Duration.ofMinutes(1))
            .GET()
            .build()
        HttpResponse<String> response = retriableHttpClient.send(request, HttpResponse.BodyHandlers.ofString())
        def json = response.body()
        def matcher = json =~ /"token"\s*:\s*"([^"]+)"/
        if (matcher.find()) {
            return matcher.group(1)
        }
        throw new IOException("Token not found in response: $json")
    }
}


