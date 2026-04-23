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

import java.net.http.HttpClient
import java.net.http.HttpResponse
import java.nio.file.AccessMode
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.model.DataLinkProvider
import io.seqera.tower.model.DataLinkDownloadUrlResponse
import io.seqera.tower.plugin.datalink.PagedDataLinkContent
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath
import spock.lang.Specification

class DataLinksResourceHandlerTest extends Specification {

    private SeqeraFileSystem fs = Mock(SeqeraFileSystem)
    private SeqeraDataLinkClient client = Mock(SeqeraDataLinkClient)
    private HttpClient http = Mock(HttpClient)
    private DataLinksResourceHandler handler = new DataLinksResourceHandler(fs, client, http)

    private static DataLinkDto dl(String id, String name, DataLinkProvider p) {
        def d = new DataLinkDto(); d.id = id; d.name = name; d.provider = p; return d
    }

    private static DataLinkItem item(String name, DataLinkItemType t, long size) {
        def i = new DataLinkItem(); i.name = name; i.type = t; i.size = size; return i
    }

    /** Single-page {@link PagedDataLinkContent} for tests. */
    private static PagedDataLinkContent pagedContent(List<DataLinkItem> items, String originalPath = null) {
        return new PagedDataLinkContent(originalPath, items, null, new PagedDataLinkContent.PageFetcher() {
            Map<String, Object> fetch(String t) { throw new IllegalStateException('no more pages') }
        })
    }

    /** Iterator over a fixed list — what {@code client.listDataLinks(...)} is expected to return. */
    private static Iterator<DataLinkDto> iter(List<DataLinkDto> list) { list.iterator() }

    private static List<Path> asList(Iterable<Path> iterable) {
        final out = new ArrayList<Path>()
        for (Path p : iterable) out.add(p)
        return out
    }

    // =====================================================
    // newInputStream — MVP
    // =====================================================

    def "newInputStream resolves (provider,name,subPath) and streams the signed URL"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/a.fq')
        def signedBody = new ByteArrayInputStream('data'.bytes)
        def httpResp = Mock(HttpResponse) {
            statusCode() >> 200
            body() >> signedBody
        }
        def urlResp = new DataLinkDownloadUrlResponse(); urlResp.url = 'https://signed/a'

        when:
        def stream = handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L) >> urlResp
        1 * http.send(_, _) >> httpResp
        stream === signedBody
    }

    def "newInputStream throws NoSuchFileException when data-link unknown"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/unknown/reads/a.fq')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        thrown(NoSuchFileException)
    }

    def "newInputStream requires trail.size >= 3 (file path, not the data-link root itself)"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        handler.newInputStream(path)

        then:
        thrown(IllegalArgumentException)
    }

    def "newInputStream surfaces signed-URL HTTP failure as IOException"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/a.fq')
        def urlResp = new DataLinkDownloadUrlResponse(); urlResp.url = 'https://signed/a'
        def httpResp = Mock(HttpResponse) {
            statusCode() >> 403
            body() >> new ByteArrayInputStream(new byte[0])
        }

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L) >> urlResp
        1 * http.send(_, _) >> httpResp
        thrown(IOException)
    }

    // =====================================================
    // list — US2 browse
    // =====================================================

    def "list at data-links/ returns distinct providers in use"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links')

        when:
        def paths = asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDataLinks(10L) >> iter([
                dl('dl-1', 'a', DataLinkProvider.AWS),
                dl('dl-2', 'b', DataLinkProvider.GOOGLE),
                dl('dl-3', 'c', DataLinkProvider.AWS)
        ])
        paths*.toString().sort() == [
                'seqera://acme/research/data-links/aws',
                'seqera://acme/research/data-links/google'
        ]
    }

    def "list at data-links/<provider>/ returns data-link names for that provider"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws')

        when:
        def paths = asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([
                dl('dl-1', 'inputs', DataLinkProvider.AWS),
                dl('dl-2', 'archive', DataLinkProvider.AWS),
                dl('dl-3', 'onGcs', DataLinkProvider.GOOGLE)
        ])
        paths*.toString() == [
                'seqera://acme/research/data-links/aws/archive',
                'seqera://acme/research/data-links/aws/inputs'
        ]
    }

    def "list at data-link root returns top-level objects"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        def paths = asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getContent('dl-1', '', 10L) >> pagedContent([
                item('reads', DataLinkItemType.FOLDER, 0),
                item('samplesheet.csv', DataLinkItemType.FILE, 42)
        ])
        paths*.toString() == [
                'seqera://acme/research/data-links/aws/inputs/reads',
                'seqera://acme/research/data-links/aws/inputs/samplesheet.csv'
        ]
    }

    def "list at deep sub-path browses the correct sub-path"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads')

        when:
        def paths = asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getContent('dl-1', 'reads', 10L) >> pagedContent([
                item('a.fq', DataLinkItemType.FILE, 1),
                item('b.fq', DataLinkItemType.FILE, 2)
        ])
        paths*.toString() == [
                'seqera://acme/research/data-links/aws/inputs/reads/a.fq',
                'seqera://acme/research/data-links/aws/inputs/reads/b.fq'
        ]
    }

    // =====================================================
    // readAttributes
    // =====================================================

    def "readAttributes at data-links/ resource-type dir reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        attr.directory
        !attr.regularFile
    }

    def "readAttributes at data-link root reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        attr.directory
    }

    def "readAttributes on a file sub-path reports file with size"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/a.fq')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getContent('dl-1', 'reads/a.fq', 10L) >> pagedContent([
                item('a.fq', DataLinkItemType.FILE, 123)
        ])
        attr.regularFile
        attr.size() == 123L
    }

    def "readAttributes on a directory sub-path reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getContent('dl-1', 'reads', 10L) >> pagedContent(
                [item('a.fq', DataLinkItemType.FILE, 1), item('b.fq', DataLinkItemType.FILE, 2)],
                'reads/')
        attr.directory
    }

    // =====================================================
    // error paths — US3
    // =====================================================

    def "list at data-links/<provider>/ throws when no data-links for that provider"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/azure')

        when:
        asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'x', DataLinkProvider.AWS)])
        def ex = thrown(NoSuchFileException)
        ex.reason?.toLowerCase()?.contains('no data-links')
    }

    def "unknown data-link under a known provider throws"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/ghost/a.fq')

        when:
        asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        thrown(NoSuchFileException)
    }

    def "missing sub-path inside a data-link surfaces as NoSuchFileException"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/does/not/exist')

        when:
        handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> iter([dl('dl-1', 'inputs', DataLinkProvider.AWS)])
        1 * client.getContent('dl-1', 'does/not/exist', 10L) >> pagedContent([])
        thrown(NoSuchFileException)
    }

    def "checkAccess with WRITE is rejected"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/a.fq')

        when:
        handler.checkAccess(path, AccessMode.WRITE)

        then:
        thrown(java.nio.file.AccessDeniedException)
    }
}
