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

package io.seqera.tower.plugin.datalink

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import groovy.json.JsonOutput
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.plugin.TowerClient
import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 * Tests for {@link SeqeraDataLinkClient} using a spy {@link TowerClient}.
 */
class SeqeraDataLinkClientTest extends Specification {

    private static final String EP = 'https://api.example.com'

    private TowerClient tower() {
        def tc = Spy(TowerClient)
        tc.@endpoint = EP
        return tc
    }

    private static TowerClient.Response ok(String body) { new TowerClient.Response(200, body) }
    private static TowerClient.Response err(int code)    { new TowerClient.Response(code, "error $code") }

    private static List<DataLinkDto> drain(Iterator<DataLinkDto> it) {
        final out = new ArrayList<DataLinkDto>()
        while (it.hasNext()) out.add(it.next())
        return out
    }

    // ---- listDataLinks ----

    def "listDataLinks yields DTOs lazily for a single page"() {
        given:
        def body = JsonOutput.toJson([dataLinks: [
                [id: 'dl-1', name: 'inputs',  provider: 'aws',    resourceRef: 's3://bucket/'],
                [id: 'dl-2', name: 'archive', provider: 'google', resourceRef: 'gs://bucket/']
        ], totalSize: 2])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=0") >> ok(body)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def list = drain(client.listDataLinks(10L))

        then:
        list.size() == 2
        list[0].id == 'dl-1'
        list[1].provider.toString() == 'google'
    }

    def "listDataLinks exhausts pagination across multiple pages"() {
        given:
        def p1 = JsonOutput.toJson([dataLinks: [[id: 'dl-1', name: 'a', provider: 'aws']], totalSize: 3])
        def p2 = JsonOutput.toJson([dataLinks: [[id: 'dl-2', name: 'b', provider: 'aws']], totalSize: 3])
        def p3 = JsonOutput.toJson([dataLinks: [[id: 'dl-3', name: 'c', provider: 'aws']], totalSize: 3])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=0") >> ok(p1)
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=1") >> ok(p2)
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=2") >> ok(p3)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def list = drain(client.listDataLinks(10L))

        then:
        list*.id == ['dl-1', 'dl-2', 'dl-3']
    }

    def "listDataLinks short-circuits — only fetches enough pages to satisfy the consumer"() {
        given:
        def p1 = JsonOutput.toJson([dataLinks: [[id: 'dl-1', name: 'a', provider: 'aws']], totalSize: 5])
        def tc = tower()
        def client = new SeqeraDataLinkClient(tc)

        when:
        def it = client.listDataLinks(10L)
        def first = it.next()

        then:
        1 * tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=0") >> ok(p1)
        0 * tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=1")
        first.id == 'dl-1'
    }

    def "listDataLinks returns empty iterator when workspace has none"() {
        given:
        def body = JsonOutput.toJson([dataLinks: [], totalSize: 0])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=0") >> ok(body)
        def client = new SeqeraDataLinkClient(tc)

        expect:
        !client.listDataLinks(10L).hasNext()
    }

    // ---- getContent ----

    def "getContent on a sub-path uses /browse/{path}"() {
        given:
        def body = JsonOutput.toJson([
                originalPath: 'reads/',
                objects: [
                        [name: 'a.fq', type: 'FILE', size: 123, mimeType: 'application/gzip'],
                        [name: 'b.fq', type: 'FILE', size: 456, mimeType: 'application/gzip']
                ]])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links/dl-1/browse/reads/?workspaceId=10") >> ok(body)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def resp = client.getContent('dl-1', 'reads/', 10L)

        then:
        resp.originalPath == 'reads/'
        resp.firstPage.size() == 2
        resp.firstPage[0].name == 'a.fq'
        resp.firstPage[0].size == 123L
    }

    def "getContent at the data-link root uses /browse (no path)"() {
        given:
        def body = JsonOutput.toJson([originalPath: '', objects: [[name: 'a', type: 'FILE', size: 1]]])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links/dl-1/browse?workspaceId=10") >> ok(body)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def resp = client.getContent('dl-1', null, 10L)

        then:
        resp.firstPage*.name == ['a']
    }

    def "getContent iterator lazily paginates across pages"() {
        given:
        def p1 = JsonOutput.toJson([originalPath: '', objects: [[name: 'a', type: 'FILE', size: 1]], nextPageToken: 'T2'])
        def p2 = JsonOutput.toJson([originalPath: '', objects: [[name: 'b', type: 'FILE', size: 2]]])
        def tc = tower()
        def client = new SeqeraDataLinkClient(tc)

        when: 'the caller iterates the full stream'
        def resp = client.getContent('dl-1', null, 10L)
        def names = resp.collect { it.name }

        then: 'first page fetched eagerly; second page fetched only when iterator advances past page 1'
        1 * tc.sendApiRequest("${EP}/data-links/dl-1/browse?workspaceId=10") >> ok(p1)
        1 * tc.sendApiRequest("${EP}/data-links/dl-1/browse?workspaceId=10&nextPageToken=T2") >> ok(p2)
        names == ['a', 'b']
    }

    def "getContent does not fetch page 2 if the caller only consumes the first page"() {
        given:
        def p1 = JsonOutput.toJson([originalPath: '', objects: [[name: 'a', type: 'FILE', size: 1]], nextPageToken: 'T2'])
        def tc = tower()
        def client = new SeqeraDataLinkClient(tc)

        when: 'caller only reads firstPage metadata without iterating'
        def resp = client.getContent('dl-1', null, 10L)
        def first = resp.firstPage

        then:
        1 * tc.sendApiRequest("${EP}/data-links/dl-1/browse?workspaceId=10") >> ok(p1)
        0 * tc.sendApiRequest("${EP}/data-links/dl-1/browse?workspaceId=10&nextPageToken=T2")
        first*.name == ['a']
    }

    // ---- getDownloadUrl ----

    def "getDownloadUrl returns the signed URL from /generate-download-url"() {
        given:
        def tc = tower()
        def expectedUrl = "${EP}/data-links/dl-1/generate-download-url?workspaceId=10&filePath=" + URLEncoder.encode('reads/a.fq', 'UTF-8')
        tc.sendApiRequest(expectedUrl) >> ok(JsonOutput.toJson([url: 'https://signed']))
        def client = new SeqeraDataLinkClient(tc)

        when:
        def resp = client.getDownloadUrl('dl-1', 'reads/a.fq', 10L)

        then:
        resp.url == 'https://signed'
    }

    // ---- error mapping ----

    def "401 throws AbortOperationException"() {
        given:
        def tc = tower()
        tc.sendApiRequest(_) >> err(401)
        def client = new SeqeraDataLinkClient(tc)

        when:
        drain(client.listDataLinks(10L))

        then:
        thrown(AbortOperationException)
    }

    def "403 throws AccessDeniedException"() {
        given:
        def tc = tower()
        tc.sendApiRequest(_) >> err(403)
        def client = new SeqeraDataLinkClient(tc)

        when:
        client.getContent('dl-1', '', 10L)

        then:
        thrown(AccessDeniedException)
    }

    def "404 throws NoSuchFileException"() {
        given:
        def tc = tower()
        tc.sendApiRequest(_) >> err(404)
        def client = new SeqeraDataLinkClient(tc)

        when:
        client.getDownloadUrl('dl-1', 'missing', 10L)

        then:
        thrown(NoSuchFileException)
    }
}
