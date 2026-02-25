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
import java.nio.file.LinkOption
import java.nio.file.attribute.BasicFileAttributes

import spock.lang.TempDir
import spock.lang.Unroll

import static com.github.tomakehurst.wiremock.client.WireMock.*

/**
 * Integration tests that exercise the full dataset:// flow:
 * URI string → DatasetPathFactory → DatasetPath → DatasetResolver → resolved Path → I/O
 *
 * Uses WireMock for Platform API and local file:// paths as the resolved "cloud" storage.
 *
 * @author Edmund Miller
 */
class DatasetIntegrationTest extends DatasetWireMockSpec {

    @TempDir
    File tempDir

    def 'should resolve dataset path and read file contents'() {
        given: 'a local file simulating cloud storage'
        def csvFile = makeFile('samples.csv', 'sample,fastq_1,fastq_2\nSAMPLE1,s1_R1.fq.gz,s1_R2.fq.gz\n')

        and: 'Platform API mocks'
        mockSession(workspaceId: '100')
        stubDatasets([[id: '42', name: 'my-samplesheet', mediaType: 'text/csv']], '100')
        stubDatasetVersions('42', [[version: 1, url: csvFile.toURI().toString(), fileName: 'samples.csv']], '100')

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

    @Unroll
    def 'should resolve #scenario'() {
        given:
        mockSession(workspaceId: '100')
        stubDatasets([[id: '42', name: 'my-data']])

        and:
        def versions = versionRows.collect { row ->
            def file = makeFile(row.fileName as String, row.content as String)
            [version: row.version, url: file.toURI().toString()]
        }
        stubDatasetVersions('42', versions)

        when:
        def path = DatasetResolver.resolve('my-data', requestedVersion)

        then:
        Files.readString(path) == expectedContent

        where:
        scenario                         | versionRows                                                                                                                           | requestedVersion | expectedContent
        'specific dataset version'       | [[version: 1, fileName: 'v1.csv', content: 'version1'], [version: 2, fileName: 'v2.csv', content: 'version2']]                    | '1'              | 'version1'
        'latest dataset version'         | [[version: 1, fileName: 'v1.csv', content: 'old'], [version: 3, fileName: 'v3.csv', content: 'latest']]                            | null             | 'latest'
    }

    @Unroll
    def 'should delegate #operation through provider'() {
        given:
        def dataFile = makeFile('data.csv', fileContent)

        and:
        mockSession(workspaceId: '100')
        stubDatasets([[id: '7', name: 'test-ds']])
        stubDatasetVersions('7', [[version: 1, url: dataFile.toURI().toString()]])

        and:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test-ds'))

        when:
        def result = operationFn.call(provider, path)

        then:
        assert assertFn.call(result)

        where:
        operation         | fileContent         | operationFn                                                                                                                            | assertFn
        'newInputStream'  | 'col1,col2\na,b\n' | { DatasetFileSystemProvider providerRef, dsPath -> providerRef.newInputStream(dsPath).text }                                          | { value -> value == 'col1,col2\na,b\n' }
        'readAttributes'  | 'hello'             | { DatasetFileSystemProvider providerRef, dsPath -> providerRef.readAttributes(dsPath, BasicFileAttributes, new LinkOption[0]) }      | { attrs -> attrs.size() == 5 && !attrs.isDirectory() && attrs.isRegularFile() }
    }

    def 'should cache resolved path across multiple reads'() {
        given:
        def dataFile = makeFile('data.csv', 'cached')

        and:
        mockSession(workspaceId: '100')
        stubDatasets([[id: '7', name: 'test-ds']])
        stubDatasetVersions('7', [[version: 1, url: dataFile.toURI().toString()]])

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

    private File makeFile(String name, String content) {
        def file = new File(tempDir, name)
        file.text = content
        return file
    }
}
