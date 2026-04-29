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
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import io.seqera.tower.model.DataLinkCredentials
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.model.DataLinkProvider
import io.seqera.tower.model.DataLinkDownloadUrlResponse
import io.seqera.tower.plugin.datalink.PagedIterable
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath
import spock.lang.Specification

class DataLinksResourceHandlerTest extends Specification {

    private SeqeraFileSystem fs = Mock(SeqeraFileSystem)
    private SeqeraDataLinkClient client = Mock(SeqeraDataLinkClient)
    private HttpClient http = Mock(HttpClient)
    private DataLinksResourceHandler handler = new DataLinksResourceHandler(fs, client, http)

    private static DataLinkDto dl(String id, String name, DataLinkProvider p, String credId = null) {
        def d = new DataLinkDto()
        d.id = id; d.name = name; d.provider = p
        if (credId) {
            def c = new DataLinkCredentials(); c.id = credId
            d.credentials = [c]
        }
        return d
    }

    private static DataLinkItem item(String name, DataLinkItemType t, long size) {
        def i = new DataLinkItem(); i.name = name; i.type = t; i.size = size; return i
    }

    private static PagedIterable<DataLinkItem> pagedContent(List<DataLinkItem> items) {
        // Single-page test fixture: firstPage already loaded, no further pages.
        return new PagedIterable<DataLinkItem>(items, true, new PagedIterable.NextPageFetcher<DataLinkItem>() {
            PagedIterable.Page<DataLinkItem> fetch() {
                throw new IllegalStateException('no more pages')
            }
        })
    }

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
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L, null) >> urlResp
        1 * http.send(_, _) >> httpResp
        stream === signedBody
    }

    def "newInputStream forwards credentialsId from the data-link's credentials"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/a.fq')
        def urlResp = new DataLinkDownloadUrlResponse(); urlResp.url = 'https://signed/a'
        def httpResp = Mock(HttpResponse) {
            statusCode() >> 200
            body() >> new ByteArrayInputStream('x'.bytes)
        }

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS, 'cred-42')
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L, 'cred-42') >> urlResp
        1 * http.send(_, _) >> httpResp
    }

    def "newInputStream throws NoSuchFileException when data-link unknown"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/unknown/reads/a.fq')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.getDataLink(10L, 'aws', 'unknown') >> null
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
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L, null) >> urlResp
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
        1 * client.getDataLinkProviders(10L) >> (['aws', 'google'] as Set)
        paths*.toString() == [
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

    def "list at data-link root returns top-level objects with cached attributes"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        def paths = asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', '', 10L, null) >> pagedContent([
                item('reads', DataLinkItemType.FOLDER, 0),
                item('samplesheet.csv', DataLinkItemType.FILE, 42)
        ])
        paths*.toString() == [
                'seqera://acme/research/data-links/aws/inputs/reads',
                'seqera://acme/research/data-links/aws/inputs/samplesheet.csv'
        ]
        // Attributes attached without follow-up API calls
        (paths[0] as SeqeraPath).cachedAttributes.directory
        (paths[1] as SeqeraPath).cachedAttributes.regularFile
        (paths[1] as SeqeraPath).cachedAttributes.size() == 42L
    }

    def "list forwards credentialsId to getContent when data-link has credentials"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS, 'cred-42')
        1 * client.getContent('dl-1', '', 10L, 'cred-42') >> pagedContent([])
    }

    def "list at deep sub-path browses the correct sub-path"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads')

        when:
        def paths = asList(handler.list(path))

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', 'reads', 10L, null) >> pagedContent([
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

    def "readAttributes at data-links/<provider>/ reports directory when provider exists"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLinkProviders(10L) >> (['aws', 'google'] as Set)
        attr.directory
    }

    def "readAttributes at data-links/<provider>/ throws when the provider has no data-links"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/azure')

        when:
        handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLinkProviders(10L) >> (['aws'] as Set)
        def ex = thrown(NoSuchFileException)
        ex.reason?.contains("No data-links for provider 'azure'")
    }

    def "readAttributes at data-link root reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        attr.directory
    }

    def "readAttributes on a deep file sub-path queries the parent dir and finds the entry"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/a.fq')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', 'reads', 10L, null, 'a.fq') >> pagedContent([
                item('a.fq', DataLinkItemType.FILE, 123),
                item('b.fq', DataLinkItemType.FILE, 456)
        ])
        attr.regularFile
        attr.size() == 123L
    }

    def "readAttributes on a top-level file queries the data-link root and finds the entry"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/samplesheet.csv')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', '', 10L, null, 'samplesheet.csv') >> pagedContent([
                item('reads', DataLinkItemType.FOLDER, 0),
                item('samplesheet.csv', DataLinkItemType.FILE, 999)
        ])
        attr.regularFile
        attr.size() == 999L
    }

    def "readAttributes on a directory sub-path queries the parent dir and reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', '', 10L, null, 'reads') >> pagedContent([
                item('reads', DataLinkItemType.FOLDER, 0),
                item('samplesheet.csv', DataLinkItemType.FILE, 1)
        ])
        attr.directory
    }

    def "readAttributes on a non-existent sub-path throws NoSuchFileException"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/missing.fq')

        when:
        handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', 'reads', 10L, null, 'missing.fq') >> pagedContent([
                item('a.fq', DataLinkItemType.FILE, 1)
        ])
        thrown(NoSuchFileException)
    }

    def "readAttributes distinguishes a directory and a same-named file inside it (deferred #1 case)"() {
        // Directory 'foo' contains a file 'foo'. Path .../bucket/foo refers to the dir.
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/bucket/foo')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'bucket') >> dl('dl-1', 'bucket', DataLinkProvider.AWS)
        // root listing shows 'foo' as a folder
        1 * client.getContent('dl-1', '', 10L, null, 'foo') >> pagedContent([
                item('foo', DataLinkItemType.FOLDER, 0)
        ])
        attr.directory

        when: 'now query the inner file at .../bucket/foo/foo'
        def innerPath = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/bucket/foo/foo')
        def innerAttr = handler.readAttributes(innerPath)

        then: 'parent dir is foo, lastSeg is foo, found as FILE'
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'bucket') >> dl('dl-1', 'bucket', DataLinkProvider.AWS)
        1 * client.getContent('dl-1', 'foo', 10L, null, 'foo') >> pagedContent([
                item('foo', DataLinkItemType.FILE, 42)
        ])
        innerAttr.regularFile
        innerAttr.size() == 42L
    }

    def "readAttributes returns NoSuchFileException when the data-link is missing"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/ghost/a.fq')

        when:
        handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'ghost') >> null
        def ex = thrown(NoSuchFileException)
        ex.reason?.contains("Data-link 'ghost' not found")
    }

    def "readAttributes short-circuits when path has cached attributes (no API call)"() {
        given:
        def attrs = new io.seqera.tower.plugin.fs.SeqeraFileAttributes(99L, java.time.Instant.EPOCH, java.time.Instant.EPOCH, 'key')
        def parent = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads')
        def path = parent.resolveWithAttributes('a.fq', attrs)

        when:
        def got = handler.readAttributes(path)

        then:
        0 * fs.resolveWorkspaceId(_, _)
        0 * client.getDataLink(_, _, _)
        0 * client.getContent(_, _, _, _)
        got === attrs
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
        1 * client.getDataLink(10L, 'aws', 'ghost') >> null
        thrown(NoSuchFileException)
    }

    def "missing sub-path inside a data-link surfaces as NoSuchFileException"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/does/not/exist')

        when:
        handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.getDataLink(10L, 'aws', 'inputs') >> dl('dl-1', 'inputs', DataLinkProvider.AWS)
        // parent listing of 'does/not' has no 'exist' entry
        1 * client.getContent('dl-1', 'does/not', 10L, null, 'exist') >> pagedContent([])
        thrown(NoSuchFileException)
    }
}
