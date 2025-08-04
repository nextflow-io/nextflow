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
import org.pf4j.update.SimpleFileDownloader

import java.nio.file.Files
import java.nio.file.Path
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

        Path destination = Files.createTempDirectory("pf4j-update-downloader");
        destination.toFile().deleteOnExit();

        String path = fileUrl.getPath();
        String fileName = path.substring(path.lastIndexOf('/') + 1);
        Path file = destination.resolve(fileName);
        HttpURLConnection conn = (HttpURLConnection) fileUrl.openConnection()
        conn.instanceFollowRedirects = true

        if (conn.responseCode == HttpURLConnection.HTTP_UNAUTHORIZED) {
            def wwwAuth = conn.getHeaderField("WWW-Authenticate")
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
                def authConn = (HttpURLConnection) fileUrl.openConnection()
                authConn.setRequestProperty("Authorization", "Bearer $token")
                authConn.instanceFollowRedirects = true

                authConn.inputStream.withStream { input ->
                    file.withOutputStream { out -> out << input }
                }

                return file
            }
        }

        // Fallback to default behavior
        conn.inputStream.withStream { input ->
            file.withOutputStream { out -> out << input }
        }
        return file
    }

    private String fetchToken(String tokenUrl) {
        def conn = (HttpURLConnection) URI.create(tokenUrl).toURL().openConnection()
        conn.setRequestProperty("Accept", "application/json")

        def json = conn.inputStream.getText("UTF-8")
        def matcher = json =~ /"token"\s*:\s*"([^"]+)"/
        if (matcher.find()) {
            return matcher.group(1)
        }
        throw new IOException("Token not found in response: $json")
    }
}


