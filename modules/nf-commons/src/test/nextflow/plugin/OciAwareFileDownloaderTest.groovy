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
import nextflow.BuildInfo
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
        downloadedFile.fileName.toString().endsWith(".zip")

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

    def 'should include User-Agent header in all HTTP requests'() {
        given:
        def expectedUserAgent = "Nextflow/$BuildInfo.version"
        def fileContent = "test plugin content"
        
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .withHeader("User-Agent", equalTo(expectedUserAgent))
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent
        
        // Verify the request was made with correct User-Agent
        wiremock.verify(getRequestedFor(urlPathEqualTo("/plugin.zip"))
            .withHeader("User-Agent", equalTo(expectedUserAgent))
        )

        cleanup:
        Files.deleteIfExists(downloadedFile)
    }

    def 'should include X-Nextflow-Version header in all HTTP requests'() {
        given:
        def expectedVersion = BuildInfo.version
        def fileContent = "test plugin content"
        
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .withHeader("X-Nextflow-Version", equalTo(expectedVersion))
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent
        
        // Verify the request was made with correct X-Nextflow-Version
        wiremock.verify(getRequestedFor(urlPathEqualTo("/plugin.zip"))
            .withHeader("X-Nextflow-Version", equalTo(expectedVersion))
        )

        cleanup:
        Files.deleteIfExists(downloadedFile)
    }

    def 'should include both User-Agent and X-Nextflow-Version headers in token requests'() {
        given:
        def expectedUserAgent = "Nextflow/$BuildInfo.version"
        def expectedVersion = BuildInfo.version
        def fileContent = "authenticated plugin archive content"
        def tokenResponse = '{"token": "test-bearer-token-12345"}'
        def authServerUrl = "${wiremock.baseUrl()}/token"
        def wwwAuthHeader = "Bearer realm=\"${authServerUrl}\",service=\"registry.example.com\",scope=\"repository:plugins/nf-test:pull\""

        // Token endpoint expects both headers
        wiremock.stubFor(get(urlPathEqualTo("/token"))
            .withQueryParam('service', equalTo('registry.example.com'))
            .withQueryParam('scope', equalTo('repository:plugins/nf-test:pull'))
            .withHeader("User-Agent", equalTo(expectedUserAgent))
            .withHeader("X-Nextflow-Version", equalTo(expectedVersion))
            .willReturn(okJson(tokenResponse))
        )
        
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(unauthorized().withHeader("WWW-Authenticate", wwwAuthHeader))
        )

        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .withHeader("Authorization", equalTo("Bearer test-bearer-token-12345"))
            .withHeader("User-Agent", equalTo(expectedUserAgent))
            .withHeader("X-Nextflow-Version", equalTo(expectedVersion))
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent
        
        // Verify token request had both headers
        wiremock.verify(getRequestedFor(urlPathEqualTo("/token"))
            .withHeader("User-Agent", equalTo(expectedUserAgent))
            .withHeader("X-Nextflow-Version", equalTo(expectedVersion))
        )
        
        // Verify authenticated download request had both headers
        wiremock.verify(getRequestedFor(urlPathEqualTo("/plugin.zip"))
            .withHeader("Authorization", equalTo("Bearer test-bearer-token-12345"))
            .withHeader("User-Agent", equalTo(expectedUserAgent))
            .withHeader("X-Nextflow-Version", equalTo(expectedVersion))
        )

        cleanup:
        Files.deleteIfExists(downloadedFile)
    }

    def 'should handle redirects automatically'() {
        given:
        def fileContent = "redirected plugin content"
        def redirectUrl = "${wiremock.baseUrl()}/redirected/plugin.zip"
        
        // Initial request returns redirect
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(temporaryRedirect(redirectUrl))
        )
        
        // Redirected URL returns the file
        wiremock.stubFor(get(urlPathEqualTo("/redirected/plugin.zip"))
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent
        
        // Verify both requests were made
        wiremock.verify(getRequestedFor(urlPathEqualTo("/plugin.zip")))
        wiremock.verify(getRequestedFor(urlPathEqualTo("/redirected/plugin.zip")))

        cleanup:
        Files.deleteIfExists(downloadedFile)
    }

    def 'should detect and prevent redirect loops'() {
        given:
        def redirectUrl1 = "${wiremock.baseUrl()}/redirect1"
        def redirectUrl2 = "${wiremock.baseUrl()}/redirect2"
        
        // Create a redirect loop
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(temporaryRedirect(redirectUrl1))
        )
        wiremock.stubFor(get(urlPathEqualTo("/redirect1"))
            .willReturn(temporaryRedirect(redirectUrl2))
        )
        wiremock.stubFor(get(urlPathEqualTo("/redirect2"))
            .willReturn(temporaryRedirect("${wiremock.baseUrl()}/plugin.zip"))
        )

        when:
        downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        def ex = thrown(IOException)
        ex.message.contains("Redirect loop detected")
    }

    def 'should handle authentication after redirects'() {
        given:
        def fileContent = "authenticated after redirect content"
        def tokenResponse = '{"token": "redirect-auth-token"}'
        def authServerUrl = "${wiremock.baseUrl()}/token"
        def wwwAuthHeader = "Bearer realm=\"${authServerUrl}\",service=\"registry.example.com\",scope=\"repository:plugins/nf-test:pull\""
        def redirectUrl = "${wiremock.baseUrl()}/auth/plugin.zip"
        
        // Initial request redirects
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(temporaryRedirect(redirectUrl))
        )
        
        // Redirected URL requires authentication
        wiremock.stubFor(get(urlPathEqualTo("/auth/plugin.zip"))
            .willReturn(unauthorized().withHeader("WWW-Authenticate", wwwAuthHeader))
        )
        
        // Token endpoint
        wiremock.stubFor(get(urlPathEqualTo("/token"))
            .withQueryParam('service', equalTo('registry.example.com'))
            .withQueryParam('scope', equalTo('repository:plugins/nf-test:pull'))
            .willReturn(okJson(tokenResponse))
        )
        
        // Authenticated request succeeds
        wiremock.stubFor(get(urlPathEqualTo("/auth/plugin.zip"))
            .withHeader("Authorization", equalTo("Bearer redirect-auth-token"))
            .willReturn(ok().withBody(fileContent))
        )

        when:
        def downloadedFile = downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        downloadedFile != null
        Files.exists(downloadedFile)
        downloadedFile.text == fileContent
        
        // Verify the sequence: initial request, redirect, auth challenge, token request, authenticated download
        wiremock.verify(getRequestedFor(urlPathEqualTo("/plugin.zip")))
        wiremock.verify(2, getRequestedFor(urlPathEqualTo("/auth/plugin.zip"))) // Auth challenge request + authenticated download
        wiremock.verify(getRequestedFor(urlPathEqualTo("/token")))

        cleanup:
        Files.deleteIfExists(downloadedFile)
    }

    def 'should throw IOException for HTTP error status codes'() {
        given:
        wiremock.stubFor(get(urlPathEqualTo("/plugin.zip"))
            .willReturn(serverError())
        )

        when:
        downloader.downloadFileHttp(new URL("${wiremock.baseUrl()}/plugin.zip"))

        then:
        def ex = thrown(IOException)
        ex.message.contains("HTTP error 500")
    }
}