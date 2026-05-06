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

import java.nio.file.AccessDeniedException
import java.nio.file.DirectoryStream
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.InvalidPathException
import java.nio.file.NoSuchFileException
import java.nio.file.attribute.BasicFileAttributes

import groovy.json.JsonOutput
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 * Tests for {@link SeqeraFileSystemProvider} using a mock {@link TowerClient}.
 */
class SeqeraFileSystemProviderTest extends Specification {

    private static final String ENDPOINT = 'https://api.example.com'

    private TowerClient spyTower() {
        def tc = Spy(TowerClient)
        tc.@endpoint = ENDPOINT
        return tc
    }

    private static TowerClient.Response ok(String body) {
        new TowerClient.Response(200, body)
    }

    private static TowerClient.Response error(int code) {
        new TowerClient.Response(code, "error $code")
    }

    private SeqeraFileSystem buildFs(TowerClient tc) {
        final client = new SeqeraDatasetClient(tc)
        final provider = new SeqeraFileSystemProvider()
        return new SeqeraFileSystem(provider, client)
    }

    private static String userInfoJson() {
        JsonOutput.toJson([user: [id: 42L, userName: 'testuser']])
    }

    private static String workspacesJson() {
        JsonOutput.toJson([orgsAndWorkspaces: [
            [orgId: 1L, orgName: 'acme', workspaceId: 10L, workspaceName: 'research', workspaceFullName: 'acme/research']
        ]])
    }

    private static String datasetsJson() {
        JsonOutput.toJson([datasets: [
            [id: 'ds-1', name: 'samples', version: 2L, mediaType: 'text/csv',
             workspaceId: 10L,
             dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z']
        ], totalSize: 1])
    }

    private static String versionsJson() {
        JsonOutput.toJson([versions: [
            [datasetId: 'ds-1', version: 1L, fileName: 'samples.csv',
             mediaType: 'text/csv', hasHeader: true, dateCreated: '2024-01-01T00:00:00Z', disabled: false],
            [datasetId: 'ds-1', version: 2L, fileName: 'samples_v2.csv',
             mediaType: 'text/csv', hasHeader: true, dateCreated: '2024-01-02T00:00:00Z', disabled: false]
        ]])
    }

    // ---- newInputStream - latest version ----

    def "newInputStream resolves latest version and downloads correct content"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson())
        tc.sendApiRequest("${ENDPOINT}/datasets/ds-1/versions?workspaceId=10") >> ok(versionsJson())
        final csvContent = 'col1,col2\n1,2\n3,4\n'
        tc.sendStreamingRequest("${ENDPOINT}/datasets/ds-1/v/2/n/samples_v2.csv?workspaceId=10") >> new ByteArrayInputStream(csvContent.getBytes('UTF-8'))

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        final text = fs.provider().newInputStream(path).text

        then:
        text == csvContent
    }

    // ---- newInputStream - pinned version ----

    def "newInputStream uses pinned version when @ver suffix given"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson())
        tc.sendApiRequest("${ENDPOINT}/datasets/ds-1/versions?workspaceId=10") >> ok(versionsJson())
        final csvContent = 'col1,col2\n1,2\n'
        tc.sendStreamingRequest("${ENDPOINT}/datasets/ds-1/v/1/n/samples.csv?workspaceId=10") >> new ByteArrayInputStream(csvContent.getBytes('UTF-8'))

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@1')

        when:
        final text = fs.provider().newInputStream(path).text

        then:
        text == csvContent
    }

    // ---- newInputStream - missing dataset ----

    def "newInputStream throws NoSuchFileException for unknown dataset"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >>
            ok(JsonOutput.toJson([datasets: [], totalSize: 0]))

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/missing-dataset')

        when:
        fs.provider().newInputStream(path)

        then:
        thrown(NoSuchFileException)
    }

    // ---- newInputStream - pinned version not found ----

    def "newInputStream throws NoSuchFileException for unknown pinned version"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson())
        tc.sendApiRequest("${ENDPOINT}/datasets/ds-1/versions?workspaceId=10") >> ok(versionsJson())

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@99')

        when:
        fs.provider().newInputStream(path)

        then:
        thrown(NoSuchFileException)
    }

    // ---- readAttributes ----

    def "readAttributes returns directory attributes for depth < 4"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research')

        when:
        final attrs = fs.provider().readAttributes(path, java.nio.file.attribute.BasicFileAttributes)

        then:
        attrs.isDirectory()
        !attrs.isRegularFile()
    }

    def "readAttributes returns file attributes for dataset path"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson())

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        final attrs = fs.provider().readAttributes(path, BasicFileAttributes)

        then:
        !attrs.isDirectory()
        attrs.isRegularFile()
    }

    // ---- newDirectoryStream (T023) ----

    def "newDirectoryStream on root returns org names"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())

        final fs = buildFs(tc)
        final root = new SeqeraPath(fs, 'seqera://')

        when:
        def entries = fs.provider().newDirectoryStream(root, null).toList()

        then:
        entries.size() == 1
        entries[0].toString() == 'seqera://acme'
    }

    def "newDirectoryStream on org returns workspace names"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())

        final fs = buildFs(tc)
        final orgPath = new SeqeraPath(fs, 'seqera://acme')

        when:
        def entries = fs.provider().newDirectoryStream(orgPath, null).toList()

        then:
        entries.size() == 1
        entries[0].toString() == 'seqera://acme/research'
    }

    def "newDirectoryStream on workspace returns datasets resource type"() {
        given:
        def tc = spyTower()
        final fs = buildFs(tc)
        final wsPath = new SeqeraPath(fs, 'seqera://acme/research')

        when:
        def entries = fs.provider().newDirectoryStream(wsPath, null).toList()

        then:
        entries.size() == 1
        entries[0].toString() == 'seqera://acme/research/datasets'
    }

    def "newDirectoryStream on datasets dir returns dataset names"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(datasetsJson())

        final fs = buildFs(tc)
        final dsDir = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def entries = fs.provider().newDirectoryStream(dsDir, null).toList()

        then:
        entries.size() == 1
        entries[0].toString() == 'seqera://acme/research/datasets/samples'
    }

    def "newDirectoryStream on datasets dir with empty workspace returns empty stream"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >>
            ok(JsonOutput.toJson([datasets: [], totalSize: 0]))

        final fs = buildFs(tc)
        final dsDir = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def entries = fs.provider().newDirectoryStream(dsDir, null).toList()

        then:
        entries.isEmpty()
    }

    def "newDirectoryStream filter is applied to entries"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >> ok(JsonOutput.toJson([datasets: [
            [id: 'ds-1', name: 'samples', version: 1L, mediaType: 'text/csv', workspaceId: 10L,
             dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z'],
            [id: 'ds-2', name: 'results', version: 1L, mediaType: 'text/csv', workspaceId: 10L,
             dateCreated: '2024-01-01T00:00:00Z', lastUpdated: '2024-01-02T00:00:00Z']
        ], totalSize: 2]))

        final fs = buildFs(tc)
        final dsDir = new SeqeraPath(fs, 'seqera://acme/research/datasets')
        final filter = { java.nio.file.Path p -> p.toString().contains('results') } as DirectoryStream.Filter

        when:
        def entries = fs.provider().newDirectoryStream(dsDir, filter).toList()

        then:
        entries.size() == 1
        entries[0].toString() == 'seqera://acme/research/datasets/results'
    }

    // ---- error scenarios (T028) ----

    def "readAttributes throws NoSuchFileException for unknown org"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://unknown-org/research')

        when:
        fs.provider().readAttributes(path, BasicFileAttributes)

        then:
        thrown(NoSuchFileException)
    }

    def "newInputStream throws NoSuchFileException containing dataset name for missing dataset"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        tc.sendApiRequest("${ENDPOINT}/datasets?workspaceId=10") >>
            ok(JsonOutput.toJson([datasets: [], totalSize: 0]))

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/missing-dataset')

        when:
        fs.provider().newInputStream(path)

        then:
        def e = thrown(NoSuchFileException)
        e.file?.contains('missing-dataset')
    }

    def "getUserInfo 401 propagates as AbortOperationException"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> new TowerClient.Response(401, 'Unauthorized')

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        fs.provider().newInputStream(path)

        then:
        thrown(AbortOperationException)
    }

    def "getUserInfo 403 propagates as AccessDeniedException"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> new TowerClient.Response(403, 'Forbidden')

        final fs = buildFs(tc)
        final path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        fs.provider().newInputStream(path)

        then:
        thrown(AccessDeniedException)
    }

    def "SeqeraPath constructor throws InvalidPathException for path with empty workspace segment"() {
        given:
        def tc = spyTower()
        final fs = buildFs(tc)

        when:
        new SeqeraPath(fs, 'seqera://acme//datasets/samples')

        then:
        thrown(InvalidPathException)
    }

    // ---- newFileSystem contract ----

    def "newFileSystem throws FileSystemAlreadyExistsException when filesystem exists"() {
        given: 'a provider with an existing filesystem'
        def tc = spyTower()
        def provider = new SeqeraFileSystemProvider()
        def fs = new SeqeraFileSystem(provider, new SeqeraDatasetClient(tc))
        provider.@fileSystem = fs

        when:
        provider.newFileSystem(new URI('seqera://test'), [:])

        then:
        thrown(FileSystemAlreadyExistsException)
    }
}
