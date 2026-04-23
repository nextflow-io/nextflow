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

package io.seqera.tower.plugin.fs.handler

import java.nio.file.AccessMode
import java.nio.file.NoSuchFileException

import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.DatasetVersionDto
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath
import spock.lang.Specification

class DatasetsResourceHandlerTest extends Specification {

    def fs = Mock(SeqeraFileSystem)
    def client = Mock(SeqeraDatasetClient)
    def handler = new DatasetsResourceHandler(fs, client)

    private static DatasetDto ds(String id, String name, long wsId = 10L) {
        def d = new DatasetDto()
        d.id = id; d.name = name; d.workspaceId = wsId
        return d
    }

    private static DatasetVersionDto ver(String dsId, long v, String file, boolean disabled = false) {
        def dv = new DatasetVersionDto()
        dv.datasetId = dsId; dv.version = v; dv.fileName = file; dv.disabled = disabled
        return dv
    }

    def "getResourceType returns 'datasets'"() {
        expect:
        handler.resourceType == 'datasets'
    }

    def "list at depth 3 returns one path per dataset"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def paths = handler.list(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDatasets(10L) >> [ds('d1', 'one'), ds('d2', 'two')]
        paths*.toString() == [
                'seqera://acme/research/datasets/one',
                'seqera://acme/research/datasets/two'
        ]
    }

    def "list result is cached across calls"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        handler.list(path)
        handler.list(path)

        then:
        2 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDatasets(10L) >> [ds('d1', 'one')]
    }

    def "newInputStream resolves latest non-disabled version when no pin"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')
        def dataset = ds('d1', 'samples')
        def v1 = ver('d1', 1, 'a.csv')
        def v2 = ver('d1', 2, 'b.csv')
        def content = new ByteArrayInputStream('x'.bytes)

        when:
        def stream = handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDatasets(10L) >> [dataset]
        1 * client.listVersions('d1', 10L) >> [v1, v2]
        1 * client.downloadDataset('d1', '2', 'b.csv', 10L) >> content
        stream === content
    }

    def "newInputStream honors @version pin"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@1')
        def dataset = ds('d1', 'samples')
        def v1 = ver('d1', 1, 'a.csv')
        def v2 = ver('d1', 2, 'b.csv')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDatasets(10L) >> [dataset]
        1 * client.listVersions('d1', 10L) >> [v1, v2]
        1 * client.downloadDataset('d1', '1', 'a.csv', 10L) >> new ByteArrayInputStream('x'.bytes)
    }

    def "newInputStream throws NoSuchFileException when dataset is missing"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/ghost')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDatasets(10L) >> [ds('d1', 'samples')]
        thrown(NoSuchFileException)
    }

    def "newInputStream throws when pinned version is unknown"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@99')
        def dataset = ds('d1', 'samples')
        def v1 = ver('d1', 1, 'a.csv')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDatasets(10L) >> [dataset]
        1 * client.listVersions('d1', 10L) >> [v1]
        thrown(NoSuchFileException)
    }

    def "newInputStream falls back to latest when pinned version is omitted"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')
        def dataset = ds('d1', 'samples')
        def enabled = ver('d1', 3, 'c.csv', false)
        def disabled = ver('d1', 4, 'd.csv', true)

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDatasets(10L) >> [dataset]
        1 * client.listVersions('d1', 10L) >> [enabled, disabled]
        1 * client.downloadDataset('d1', '3', 'c.csv', 10L) >> new ByteArrayInputStream('x'.bytes)
    }

    def "readAttributes at depth 3 reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        attr.directory
        !attr.regularFile
    }

    def "readAttributes at depth 4 returns file attributes"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDatasets(10L) >> [ds('d1', 'samples')]
        attr.regularFile
        !attr.directory
        attr.fileKey() == 'd1'
    }

    def "checkAccess rejects WRITE"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        handler.checkAccess(path, AccessMode.WRITE)

        then:
        thrown(java.nio.file.AccessDeniedException)
    }
}
