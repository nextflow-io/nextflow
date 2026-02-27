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

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import groovy.json.JsonOutput
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPOutputStream

import static com.github.tomakehurst.wiremock.client.WireMock.*
import static com.github.tomakehurst.wiremock.core.WireMockConfiguration.wireMockConfig

/**
 * Integration tests for ModuleRegistryClient using WireMock
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleRegistryClientTest extends Specification {

    @TempDir
    Path tempDir

    WireMockServer wireMock
    String url
    static final String MODULES_API_PATH = "/api/v1/modules"
    def setup() {
        wireMock = new WireMockServer(wireMockConfig().dynamicPort())
        wireMock.start()
        WireMock.configureFor("localhost", wireMock.port())
        url = "http://localhost:${wireMock.port()}/api"
    }

    def cleanup() {
        wireMock?.stop()
    }

    def 'should fetch module metadata from registry'() {
        given:
        def moduleResponse = [
            name: 'nf-core/fastqc',
            description: 'FastQC quality control',
            latest: [
                version: '1.1.0',
                createdAt: '2024-02-01T00:00:00Z'
            ]
        ]

        // Note: nf-core/fastqc is URL-encoded as nf-core%2Ffastqc
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody(JsonOutput.toJson(moduleResponse))))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        def result = client.fetchModule('nf-core/fastqc')

        then:
        result != null
        result.name == 'nf-core/fastqc'
        result.latest != null
        result.latest.version == '1.1.0'

        and: 'verify request was made'
        verify(getRequestedFor(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc')))
    }

    def 'should search modules in registry'() {
        given:
        def searchResponse = [
            query: 'fastqc',
            totalResults: 2,
            results: [
                [
                    name: 'nf-core/fastqc',
                    description: 'FastQC quality control',
                    relevanceScore: 0.95,
                    keywords: ['quality-control'],
                    tools: ['fastqc'],
                    revoked: false
                ],
                [
                    name: 'other/fastqc',
                    description: 'Another FastQC module',
                    relevanceScore: 0.75,
                    keywords: ['qc'],
                    tools: ['fastqc'],
                    revoked: false
                ]
            ]
        ]

        stubFor(get(urlPathEqualTo(MODULES_API_PATH))
            .withQueryParam('query', equalTo('fastqc'))
            .withQueryParam('limit', equalTo('10'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody(JsonOutput.toJson(searchResponse))))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        def result = client.search('fastqc', 10)

        then:
        result != null
        result.query == 'fastqc'
        result.totalResults == 2
        result.results.size() == 2
        result.results[0].name == 'nf-core/fastqc'
        // Use closeTo for Float/BigDecimal comparison
        Math.abs(result.results[0].relevanceScore - 0.95) < 0.001

        and: 'verify query parameters'
        verify(getRequestedFor(urlPathEqualTo(MODULES_API_PATH))
            .withQueryParam('query', equalTo('fastqc'))
            .withQueryParam('limit', equalTo('10')))
    }

    def 'should download module package from registry'() {
        given:
        def modulePackage = createTestModulePackage()
        def expectedChecksum = "${computeSha256(modulePackage)}"

        // Note: URL-encoded path
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc/1.0.0/download'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/gzip')
                .withHeader('X-Checksum', "sha256:${expectedChecksum}")
                .withBody(modulePackage)))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)
        def destFile = tempDir.resolve('module.tgz')

        when:
        def result = client.downloadModule('nf-core/fastqc', '1.0.0', destFile)

        then:
        result == destFile
        Files.exists(destFile)
        Files.size(destFile) == modulePackage.length

        and:
        verify(getRequestedFor(urlEqualTo(MODULES_API_PATH +'/nf-core%2Ffastqc/1.0.0/download')))
    }

    def 'should successfully fetch module without authentication'() {
        given:
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody(JsonOutput.toJson([name: 'nf-core/fastqc', latest: [version: '1.0.0']]))))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        def result = client.fetchModule('nf-core/fastqc')

        then:
        result != null
        result.name == 'nf-core/fastqc'

        and: 'verify request was made'
        verify(getRequestedFor(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc')))
    }

    def 'should handle 404 not found error'() {
        given:
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Fnonexistent'))
            .willReturn(aResponse()
                .withStatus(404)
                .withBody('Module not found')))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        client.fetchModule('nf-core/nonexistent')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Unable to fetch module') || ex.message.contains('Module not found')

        and:
        verify(getRequestedFor(urlEqualTo(MODULES_API_PATH + '/nf-core%2Fnonexistent')))
    }

    def 'should handle 500 server error'() {
        given:
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc'))
            .willReturn(aResponse()
                .withStatus(500)
                .withBody('Internal Server Error')))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        client.fetchModule('nf-core/fastqc')

        then:
        thrown(AbortOperationException)
    }

    def 'should send user agent header'() {
        given:
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody(JsonOutput.toJson([name: 'nf-core/fastqc', latest: [version: '1.0.0'], releases: []]))))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        client.fetchModule('nf-core/fastqc')

        then:
        verify(getRequestedFor(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc'))
            .withHeader('User-Agent', matching('.*')))
    }

    def 'should handle empty search results'() {
        given:
        stubFor(get(urlPathEqualTo(MODULES_API_PATH))
            .withQueryParam('query', equalTo('nonexistent'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody(JsonOutput.toJson([
                    query: 'nonexistent',
                    totalResults: 0,
                    results: []
                ]))))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        def result = client.search('nonexistent', 10)

        then:
        result != null
        result.totalResults == 0
        result.results.isEmpty()
    }

    def 'should respect custom search limit'() {
        given:
        stubFor(get(urlPathEqualTo(MODULES_API_PATH))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody(JsonOutput.toJson([
                    query: 'test',
                    totalResults: 0,
                    results: []
                ]))))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        client.search('test', 25)

        then:
        verify(getRequestedFor(urlPathEqualTo(MODULES_API_PATH ))
            .withQueryParam('limit', equalTo('25')))
    }

    def 'should handle malformed JSON response'() {
        given:
        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/json')
                .withBody('not valid json {]')))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)

        when:
        client.fetchModule('nf-core/fastqc')

        then:
        thrown(AbortOperationException)
    }

    def 'should verify download includes checksum header'() {
        given:
        def modulePackage = createTestModulePackage()
        def checksum = "sha256:${computeSha256(modulePackage)}"

        stubFor(get(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc/1.0.0/download'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Type', 'application/gzip')
                .withHeader('X-Checksum', checksum)
                .withHeader('Docker-Content-Digest', checksum)
                .withBody(modulePackage)))

        and:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)
        def destFile = tempDir.resolve('module.tgz')

        when:
        def result = client.downloadModule('nf-core/fastqc', '1.0.0', destFile)

        then:
        result == destFile
        Files.exists(destFile)

        and: 'verify checksum header was present'
        verify(getRequestedFor(urlEqualTo(MODULES_API_PATH + '/nf-core%2Ffastqc/1.0.0/download')))
    }

    def 'should handle network errors gracefully'() {
        given:
        // Stub will not be set up, causing connection refused
        def config = new RegistryConfig([url: "http://localhost:9999"]) // Invalid port
        def client = new ModuleRegistryClient(config)

        when:
        client.fetchModule('nf-core/fastqc')

        then:
        thrown(AbortOperationException)
    }

    def 'should require authentication for publish'() {
        given:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)
        def publishRequest = [name: 'nf-core/mymodule', version: '1.0.0']

        when:
        client.publishModule('nf-core/mymodule', publishRequest)

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Authentication required')
    }

    def 'should handle publish failure with no auth token'() {
        given:
        def config = new RegistryConfig([url: url])
        def client = new ModuleRegistryClient(config)
        def publishRequest = [name: 'nf-core/mymodule', version: '1.0.0']

        when:
        client.publishModule('nf-core/mymodule', publishRequest)

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Authentication required')
    }

    // Helper methods

    private byte[] createTestModulePackage() {
        def baos = new ByteArrayOutputStream()

        new GZIPOutputStream(baos).withCloseable { gzos ->
            new TarArchiveOutputStream(gzos).withCloseable { tos ->
                // Add main.nf
                def mainContent = 'process TEST { script: "echo test" }'
                addTarEntry(tos, 'main.nf', mainContent.bytes)

                // Add meta.yml
                def metaContent = 'name: test\nversion: 1.0.0'
                addTarEntry(tos, 'meta.yml', metaContent.bytes)
            }
        }

        return baos.toByteArray()
    }

    private void addTarEntry(TarArchiveOutputStream tos, String name, byte[] content) {
        def entry = new TarArchiveEntry(name)
        entry.setSize(content.length)
        tos.putArchiveEntry(entry)
        tos.write(content)
        tos.closeArchiveEntry()
    }

    private String computeSha256(byte[] data) {
        def digest = java.security.MessageDigest.getInstance('SHA-256')
        def hash = digest.digest(data)
        return hash.encodeHex().toString()
    }
}
