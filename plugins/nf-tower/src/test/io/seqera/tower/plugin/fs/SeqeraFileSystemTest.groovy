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

package io.seqera.tower.plugin.fs

import java.nio.file.NoSuchFileException

import groovy.json.JsonOutput
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import spock.lang.Specification

/**
 * Tests for {@link SeqeraFileSystem} caching and workspace resolution using a mock {@link TowerClient}.
 */
class SeqeraFileSystemTest extends Specification {

    private static final String ENDPOINT = 'https://api.example.com'

    private TowerClient mockTower() {
        def tc = Mock(TowerClient)
        tc.endpoint >> ENDPOINT
        return tc
    }

    private static TowerClient.Response ok(String body) {
        new TowerClient.Response(200, body)
    }

    private static String userInfoJson() {
        JsonOutput.toJson([user: [id: 42L, userName: 'testuser']])
    }

    private static String workspacesJson() {
        JsonOutput.toJson([orgsAndWorkspaces: [
            [orgId: 1L, orgName: 'acme', workspaceId: 10L, workspaceName: 'research', workspaceFullName: 'acme/research'],
            [orgId: 1L, orgName: 'acme', workspaceId: 20L, workspaceName: 'dev', workspaceFullName: 'acme/dev'],
            [orgId: 2L, orgName: 'other', workspaceId: 30L, workspaceName: 'ws', workspaceFullName: 'other/ws']
        ]])
    }

    private SeqeraFileSystem buildFs(TowerClient tc) {
        new SeqeraFileSystem(new SeqeraFileSystemProvider(), new SeqeraDatasetClient(tc))
    }

    // ---- cache loading ----

    def "loadOrgWorkspaceCache is called only once across multiple invocations"() {
        given:
        def tc = mockTower()
        final fs = buildFs(tc)

        when:
        fs.loadOrgWorkspaceCache()
        fs.loadOrgWorkspaceCache()

        then:
        1 * tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        1 * tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
    }

    def "listOrgNames returns distinct org names from cache"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        final fs = buildFs(tc)

        when:
        def orgs = fs.listOrgNames()

        then:
        orgs.size() == 2
        orgs.contains('acme')
        orgs.contains('other')
    }

    def "listWorkspaceNames returns workspace names for the given org"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        final fs = buildFs(tc)

        when:
        def names = fs.listWorkspaceNames('acme')

        then:
        names.size() == 2
        names.containsAll(['research', 'dev'])
    }

    // ---- resolveWorkspaceId ----

    def "resolveWorkspaceId returns correct ID for known org and workspace"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        final fs = buildFs(tc)

        when:
        def id = fs.resolveWorkspaceId('acme', 'research')

        then:
        id == 10L
    }

    def "resolveWorkspaceId throws NoSuchFileException for unknown org"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        final fs = buildFs(tc)

        when:
        fs.resolveWorkspaceId('unknown-org', 'research')

        then:
        thrown(NoSuchFileException)
    }

    def "resolveWorkspaceId throws NoSuchFileException for unknown workspace within known org"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        final fs = buildFs(tc)

        when:
        fs.resolveWorkspaceId('acme', 'no-such-ws')

        then:
        thrown(NoSuchFileException)
    }

    // ---- dataset cache ----

    def "resolveDatasets populates cache and returns datasets"() {
        given:
        def tc = mockTower()
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >>
            ok(JsonOutput.toJson([datasets: [
                [id: 'ds-1', name: 'samples', version: 1L, mediaType: 'text/csv',
                 dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z']
            ], totalSize: 1]))
        final fs = buildFs(tc)

        when:
        def datasets = fs.resolveDatasets(10L)

        then:
        datasets.size() == 1
        datasets[0].name == 'samples'
    }

    def "resolveDatasets returns cached result on second call without extra API request"() {
        given:
        def tc = mockTower()
        final datasetsJson = JsonOutput.toJson([datasets: [
            [id: 'ds-1', name: 'samples', version: 1L, mediaType: 'text/csv',
             dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z']
        ], totalSize: 1])
        final fs = buildFs(tc)

        when:
        fs.resolveDatasets(10L)
        fs.resolveDatasets(10L)

        then:
        1 * tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson)
    }

    def "invalidateDatasetCache forces re-fetch on next resolveDatasets call"() {
        given:
        def tc = mockTower()
        final datasetsJson = JsonOutput.toJson([datasets: [
            [id: 'ds-1', name: 'samples', version: 1L, mediaType: 'text/csv',
             dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z']
        ], totalSize: 1])
        final fs = buildFs(tc)

        when:
        fs.resolveDatasets(10L)
        fs.invalidateDatasetCache(10L)
        fs.resolveDatasets(10L)

        then:
        2 * tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson)
    }
}
