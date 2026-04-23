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
 * Tests for {@link SeqeraFileSystem} org/workspace cache and handler registry.
 * Resource-specific caches (datasets, data-links) are tested against their handlers.
 */
class SeqeraFileSystemTest extends Specification {

    private static final String ENDPOINT = 'https://api.example.com'

    private TowerClient spyTower() {
        def tc = Spy(TowerClient)
        tc.@endpoint = ENDPOINT
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
        final fs = new SeqeraFileSystem(new SeqeraFileSystemProvider())
        fs.setOrgWorkspaceClient(new SeqeraDatasetClient(tc))
        return fs
    }

    // ---- cache loading ----

    def "loadOrgWorkspaceCache is called only once across multiple invocations"() {
        given:
        def tc = spyTower()
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
        def tc = spyTower()
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
        def tc = spyTower()
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
        def tc = spyTower()
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
        def tc = spyTower()
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
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        final fs = buildFs(tc)

        when:
        fs.resolveWorkspaceId('acme', 'no-such-ws')

        then:
        thrown(NoSuchFileException)
    }

    // ---- handler registry ----

    def "registerHandler stores and looks up by resource type"() {
        given:
        def fs = new SeqeraFileSystem(new SeqeraFileSystemProvider())
        def handler = Mock(ResourceTypeHandler) {
            getResourceType() >> 'datasets'
        }

        when:
        fs.registerHandler(handler)

        then:
        fs.getHandler('datasets') === handler
        fs.getHandler('unknown') == null
        fs.getResourceTypes() == ['datasets'] as Set
    }
}
