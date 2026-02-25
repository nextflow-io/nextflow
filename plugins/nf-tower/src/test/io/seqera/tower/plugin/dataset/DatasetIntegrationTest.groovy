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

import java.nio.file.Files

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import nextflow.Global
import nextflow.Session
import spock.lang.AutoCleanup
import spock.lang.Specification
import spock.lang.TempDir

import static com.github.tomakehurst.wiremock.client.WireMock.*

/**
 * Integration tests that exercise the full dataset:// flow:
 * URI string → DatasetPathFactory → DatasetPath → DatasetResolver → resolved Path → I/O
 *
 * Uses WireMock for Platform API and local file:// paths as the resolved "cloud" storage.
 *
 * @author Edmund Miller
 */
class DatasetIntegrationTest extends Specification {

    @AutoCleanup('stop')
    WireMockServer wireMock

    @TempDir
    File tempDir

    def setup() {
        wireMock = new WireMockServer(0)
        wireMock.start()
        WireMock.configureFor(wireMock.port())
    }

    def cleanup() {
        Global.session = null
    }

    private void mockSession() {
        def endpoint = "http://localhost:${wireMock.port()}"
        def config = [tower: [endpoint: endpoint, accessToken: 'test-token', workspaceId: '100']]
        Global.session = Mock(Session) {
            getConfig() >> config
        }
    }

    def 'should resolve dataset path and read file contents'() {
        given: 'a local file simulating cloud storage'
        def csvFile = new File(tempDir, 'samples.csv')
        csvFile.text = 'sample,fastq_1,fastq_2\nSAMPLE1,s1_R1.fq.gz,s1_R2.fq.gz\n'
        def fileUri = csvFile.toURI().toString()

        and: 'Platform API mocks'
        mockSession()
        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .withQueryParam('workspaceId', equalTo('100'))
            .willReturn(okJson("""{"datasets": [
                {"id": "42", "name": "my-samplesheet", "mediaType": "text/csv"}
            ]}""")))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/42/versions'))
            .withQueryParam('workspaceId', equalTo('100'))
            .willReturn(okJson("""{"versions": [
                {"version": 1, "url": "${fileUri}", "fileName": "samples.csv"}
            ]}""")))

        when: 'resolving the dataset path'
        def path = DatasetResolver.resolve('my-samplesheet', null)

        then: 'resolved to the local file and readable'
        path != null
        Files.exists(path)
        Files.readString(path).contains('SAMPLE1')

        and: 'correct API calls were made'
        wireMock.verify(getRequestedFor(urlPathEqualTo('/datasets'))
            .withHeader('Authorization', equalTo('Bearer test-token'))
            .withQueryParam('workspaceId', equalTo('100')))
        wireMock.verify(getRequestedFor(urlPathEqualTo('/datasets/42/versions'))
            .withQueryParam('workspaceId', equalTo('100')))
    }

    def 'should resolve specific version'() {
        given:
        def v1File = new File(tempDir, 'v1.csv')
        v1File.text = 'version1'
        def v2File = new File(tempDir, 'v2.csv')
        v2File.text = 'version2'

        and:
        mockSession()
        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "42", "name": "my-data"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/42/versions'))
            .willReturn(okJson("""{"versions": [
                {"version": 1, "url": "${v1File.toURI()}"},
                {"version": 2, "url": "${v2File.toURI()}"}
            ]}""")))

        when: 'requesting version 1'
        def path = DatasetResolver.resolve('my-data', '1')

        then:
        Files.readString(path) == 'version1'
    }

    def 'should resolve latest version when multiple exist'() {
        given:
        def v1File = new File(tempDir, 'v1.csv')
        v1File.text = 'old'
        def v3File = new File(tempDir, 'v3.csv')
        v3File.text = 'latest'

        and:
        mockSession()
        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "42", "name": "my-data"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/42/versions'))
            .willReturn(okJson("""{"versions": [
                {"version": 1, "url": "${v1File.toURI()}"},
                {"version": 3, "url": "${v3File.toURI()}"}
            ]}""")))

        when: 'requesting latest (no version specified)'
        def path = DatasetResolver.resolve('my-data', null)

        then: 'gets version 3 (highest)'
        Files.readString(path) == 'latest'
    }

    def 'should read dataset file via provider newInputStream'() {
        given:
        def csvFile = new File(tempDir, 'data.csv')
        csvFile.text = 'col1,col2\na,b\n'

        and:
        mockSession()
        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "7", "name": "test-ds"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/7/versions'))
            .willReturn(okJson("""{"versions": [
                {"version": 1, "url": "${csvFile.toURI()}"}
            ]}""")))

        and: 'a DatasetPath created via the provider'
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test-ds'))

        when:
        def content = provider.newInputStream(path).text

        then:
        content == 'col1,col2\na,b\n'
    }

    def 'should read attributes via provider'() {
        given:
        def csvFile = new File(tempDir, 'data.csv')
        csvFile.text = 'hello'

        and:
        mockSession()
        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "7", "name": "test-ds"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/7/versions'))
            .willReturn(okJson("""{"versions": [
                {"version": 1, "url": "${csvFile.toURI()}"}
            ]}""")))

        and:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test-ds'))

        when:
        def attrs = provider.readAttributes(path, java.nio.file.attribute.BasicFileAttributes, new java.nio.file.LinkOption[0])

        then:
        attrs.size() == 5
        !attrs.isDirectory()
        attrs.isRegularFile()
    }

    def 'should cache resolved path across multiple reads'() {
        given:
        def csvFile = new File(tempDir, 'data.csv')
        csvFile.text = 'cached'

        and:
        mockSession()
        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "7", "name": "test-ds"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/7/versions'))
            .willReturn(okJson("""{"versions": [
                {"version": 1, "url": "${csvFile.toURI()}"}
            ]}""")))

        and:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test-ds'))

        when: 'reading twice'
        def content1 = provider.newInputStream(path).text
        def content2 = provider.newInputStream(path).text

        then: 'both reads succeed'
        content1 == 'cached'
        content2 == 'cached'

        and: 'API was only called once (path cached on DatasetPath)'
        wireMock.verify(1, getRequestedFor(urlPathEqualTo('/datasets')))
        wireMock.verify(1, getRequestedFor(urlPathEqualTo('/datasets/7/versions')))
    }
}
