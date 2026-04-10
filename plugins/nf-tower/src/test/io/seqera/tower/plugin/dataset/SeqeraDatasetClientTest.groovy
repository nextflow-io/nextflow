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

package io.seqera.tower.plugin.dataset

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import groovy.json.JsonOutput
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.exception.ForbiddenException
import io.seqera.tower.plugin.exception.NotFoundException
import io.seqera.tower.plugin.exception.UnauthorizedException
import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 * Tests for {@link SeqeraDatasetClient} using a mock {@link TowerClient}.
 */
class SeqeraDatasetClientTest extends Specification {

    private TowerClient mockTower(String endpoint = 'https://api.example.com') {
        def tc = Mock(TowerClient)
        tc.endpoint >> endpoint
        return tc
    }
    private TowerClient spyTower(String endpoint = 'https://api.example.com') {
        def tc = Spy(TowerClient)
        tc.@endpoint = endpoint
        return tc
    }

    private static TowerClient.Response ok(String body) {
        new TowerClient.Response(200, body)
    }

    private static TowerClient.Response error(int code) {
        new TowerClient.Response(code, "error $code")
    }

    // ---- getUserInfo ----

    def "getUserInfo returns parsed user map"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest('https://api.example.com/user-info') >> ok(JsonOutput.toJson([user: [id: 42, userName: 'testuser']]))
        def client = new SeqeraDatasetClient(tc)

        when:
        def info = client.getUserInfo()

        then:
        info.id == 42
        info.userName == 'testuser'
    }

    def "getUserInfo throws AbortOperationException on 401"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest(_) >> error(401)
        def client = new SeqeraDatasetClient(tc)

        when:
        client.getUserInfo()

        then:
        thrown(AbortOperationException)
    }

    // ---- listUserWorkspacesAndOrgs ----

    def "listUserWorkspacesAndOrgs returns parsed DTOs"() {
        given:
        def body = JsonOutput.toJson([orgsAndWorkspaces: [
            [orgId: 1, orgName: 'acme', workspaceId: 10, workspaceName: 'research', workspaceFullName: 'acme/research']
        ]])
        def tc = spyTower()
        tc.sendApiRequest('https://api.example.com/user/42/workspaces') >> ok(body)
        def client = new SeqeraDatasetClient(tc)

        when:
        def list = client.listUserWorkspacesAndOrgs(42L)

        then:
        list.size() == 1
        list[0].orgName == 'acme'
        list[0].workspaceId == 10L
        list[0].workspaceName == 'research'
    }

    // ---- listDatasets ----

    def "listDatasets returns parsed DatasetDto list"() {
        given:
        def body = JsonOutput.toJson([datasets: [
            [id: 'ds-1', name: 'samples', version: 2, mediaType: 'text/csv',
             dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z']
        ], totalSize: 1])
        def tc = mockTower()
        tc.sendApiRequest('https://api.example.com/datasets?workspaceId=99') >> ok(body)
        def client = new SeqeraDatasetClient(tc)

        when:
        def list = client.listDatasets(99L)

        then:
        list.size() == 1
        list[0].id == 'ds-1'
        list[0].name == 'samples'
        list[0].version == 2L
    }

    def "listDatasets returns empty list when no datasets"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest('https://api.example.com/datasets?workspaceId=99') >>
            ok(JsonOutput.toJson([datasets: [], totalSize: 0]))
        def client = new SeqeraDatasetClient(tc)

        when:
        def list = client.listDatasets(99L)

        then:
        list.isEmpty()
    }

    // ---- listVersions ----

    def "listVersions returns parsed DatasetVersionDto list"() {
        given:
        def body = JsonOutput.toJson([versions: [
            [datasetId: 'ds-1', version: 1, fileName: 'samples.csv',
             mediaType: 'text/csv', hasHeader: true, dateCreated: '2024-01-01T00:00:00Z', disabled: false]
        ]])
        def tc = mockTower()
        tc.sendApiRequest('https://api.example.com/datasets/ds-1/versions?workspaceId=1234') >> ok(body)
        def client = new SeqeraDatasetClient(tc)

        when:
        def list = client.listVersions('ds-1', 1234)

        then:
        list.size() == 1
        list[0].version == 1L
        list[0].fileName == 'samples.csv'
        list[0].hasHeader
        !list[0].disabled
    }

    // ---- downloadDataset ----

    def "downloadDataset returns InputStream with correct content"() {
        given:
        def content = 'col1,col2\n1,2\n'
        def tc = mockTower()
        tc.sendStreamingRequest('https://api.example.com/datasets/ds-1/v/1/n/samples.csv?workspaceId=1234') >> new ByteArrayInputStream(content.getBytes('UTF-8'))
        def client = new SeqeraDatasetClient(tc)

        when:
        def stream = client.downloadDataset('ds-1', '1', 'samples.csv', 1234)

        then:
        stream.text == content
    }

    def "downloadDataset URL-encodes the filename"() {
        given:
        def tc = mockTower()
        def client = new SeqeraDatasetClient(tc)

        when:
        client.downloadDataset('ds-1', '1', 'my file.csv',1234)

        then:
        1 * tc.sendStreamingRequest('https://api.example.com/datasets/ds-1/v/1/n/my%20file.csv?workspaceId=1234') >> new ByteArrayInputStream('data'.getBytes('UTF-8'))
    }

    def "downloadDataset throws NoSuchFileException on 404"() {
        given:
        def tc = mockTower()
        tc.sendStreamingRequest(_) >> { throw new NotFoundException("not found") }
        def client = new SeqeraDatasetClient(tc)

        when:
        client.downloadDataset('ds-missing', '1', 'file.csv', 1234)

        then:
        thrown(NoSuchFileException)
    }

    def "downloadDataset throws AccessDeniedException on 403"() {
        given:
        def tc = mockTower()
        tc.sendStreamingRequest(_) >> { throw new ForbiddenException("forbidden") }
        def client = new SeqeraDatasetClient(tc)

        when:
        client.downloadDataset('ds-1', '1', 'file.csv', 1234)

        then:
        thrown(AccessDeniedException)
    }

    def "downloadDataset throws AbortOperationException on 401"() {
        given:
        def tc = mockTower()
        tc.sendStreamingRequest(_) >> { throw new UnauthorizedException("unauthorized") }
        def client = new SeqeraDatasetClient(tc)

        when:
        client.downloadDataset('ds-1', '1', 'file.csv', 1234)

        then:
        thrown(AbortOperationException)
    }

    // ---- createDataset ----

    def "createDataset posts and returns created dataset"() {
        given:
        def responseBody = JsonOutput.toJson([dataset: [id: 'ds-new', name: 'results']])
        def tc = mockTower()
        tc.sendApiRequest('https://api.example.com/datasets?workspaceId=10', [name: 'results'], 'POST') >> ok(responseBody)
        def client = new SeqeraDatasetClient(tc)

        when:
        def dto = client.createDataset(10L, 'results')

        then:
        dto.id == 'ds-new'
        dto.name == 'results'
    }
}
