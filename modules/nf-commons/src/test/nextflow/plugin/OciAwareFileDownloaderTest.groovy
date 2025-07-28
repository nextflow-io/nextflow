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
 */

package nextflow.plugin

import com.github.tomakehurst.wiremock.junit.WireMockRule
import org.junit.Rule
import spock.lang.Specification

import java.nio.file.Files

import static com.github.tomakehurst.wiremock.client.WireMock.*

class OciAwareFileDownloaderTest extends Specification {

    @Rule
    WireMockRule wiremock = new WireMockRule(0)

    OciAwareFileDownloader downloader

    def setup() {
        downloader = new OciAwareFileDownloader()
    }

    def 'should download file successfully without authentication'() {
        given:
        def fileContent = "test plugin archive content"
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent
        downloadedFile.fileName.toString() == "plugin.zip"

        cleanup:
        if (downloadedFile) Files.deleteIfExists(downloadedFile)
    }

    def 'should handle OCI token authentication when receiving 401 unauthorized'() {
        given:
        def fileContent = "authenticated plugin archive content"
        def tokenResponse = '{"token": "test-bearer-token-12345"}'
        def authServerUrl = "${wiremock.baseUrl()}/token"
        def wwwAuthHeader = "Bearer realm=\"${authServerUrl}\",service=\"registry.example.com\",scope=\"repository:plugins/nf-test:pull\""

        // Token endpoint returns authentication token

        wiremock.stubFor( get(urlPathEqualTo("/token"))
                .withQueryParam('service', equalTo('registry.example.com'))
                .withQueryParam('scope', equalTo('repository:plugins/nf-test:pull'))
                .willReturn( okJson(tokenResponse))
        )
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(unauthorized().withHeader("WWW-Authenticate", wwwAuthHeader))
        )

        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .withHeader("Authorization", equalTo("Bearer test-bearer-token-12345") )
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent

        cleanup:
        Files.deleteIfExists(downloadedFile)
    }

    def 'should throw IOException when WWW-Authenticate header is malformed'() {
        given:
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(unauthorized().withHeader("WWW-Authenticate", "Bearer invalid-header-format"))
        )

        when:
        downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        def ex = thrown(IOException)
        ex.message.contains("Invalid WWW-Authenticate header")
    }

    def 'should throw IOException when token is not found in response'() {
        given:
        def invalidTokenResponse = '{"access_token": "wrong-field-name"}'
        def authServerUrl = "${wiremock.baseUrl()}/token"
        def wwwAuthHeader = "Bearer realm=\"${authServerUrl}\",service=\"registry.example.com\",scope=\"repository:plugins/nf-test:pull\""

        wiremock.stubFor( get(urlPathEqualTo("/token"))
                .withQueryParam('service', equalTo('registry.example.com'))
                .withQueryParam('scope', equalTo('repository:plugins/nf-test:pull'))
                .willReturn( okJson(invalidTokenResponse))
        )
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(unauthorized().withHeader("WWW-Authenticate", wwwAuthHeader))
        )

        when:
        downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        def ex = thrown(IOException)
        ex.message.contains("Token not found in response")
    }
}