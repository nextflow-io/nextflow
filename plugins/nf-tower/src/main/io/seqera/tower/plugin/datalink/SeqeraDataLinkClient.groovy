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

import groovy.transform.Memoized

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DataLinkCredentials
import io.seqera.tower.model.DataLinkDownloadUrlResponse
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.model.DataLinkProvider
import io.seqera.tower.plugin.TowerClient
import nextflow.exception.AbortOperationException

/**
 * Typed client for Seqera Platform data-link API endpoints.
 *
 * Paginated endpoints return lazy iterators so callers don't materialize the
 * full result set in memory — only the current page is buffered at any time.
 */
@Slf4j
@CompileStatic
class SeqeraDataLinkClient {

    private static final int LIST_PAGE_SIZE = 100

    private final TowerClient towerClient

    SeqeraDataLinkClient(TowerClient towerClient) {
        this.towerClient = towerClient
    }

    private String getEndpoint() { towerClient.endpoint }

    /**
     * Lazy iterator over every data-link in the workspace.
     * The first page is fetched eagerly from {@code GET /data-links?workspaceId=<ws>&max=<n>&offset=<o>}
     * (so any IOException surfaces here, not on the first {@code hasNext()}); subsequent
     * pages are fetched on demand as the iterator advances.
     */
    Iterator<DataLinkDto> listDataLinks(long workspaceId) throws IOException {
        return PagedIterable.start(new DataLinkListFetcher(towerClient, endpoint, workspaceId, LIST_PAGE_SIZE, null)).iterator()
    }

    /**
     * Lazy iterator over the data-links of a single provider, filtered server-side via the
     * Platform's {@code search=provider:<provider>} keyword so the workspace's other
     * providers are never paged over. Callers should still verify {@code provider} equality
     * on each result, as the keyword search is not guaranteed to be an exact match.
     */
    Iterator<DataLinkDto> listDataLinksByProvider(long workspaceId, String provider) throws IOException {
        final search = "provider:${provider}".toString()
        return PagedIterable.start(new DataLinkListFetcher(towerClient, endpoint, workspaceId, LIST_PAGE_SIZE, search)).iterator()
    }

    /**
     * Distinct provider identifiers present in the workspace, sorted.
     * The returned set is unmodifiable; memoized per workspace.
     */
    @Memoized
    Set<String> getDataLinkProviders(long workspaceId) {
        final providers = new TreeSet<String>()
        final Iterator<DataLinkDto> it = listDataLinks(workspaceId)
        while (it.hasNext()) {
            final p = it.next().provider?.toString()
            if (p)
                providers.add(p)
        }
        return Collections.unmodifiableSet(providers)
    }

    /**
     * Resolve a data-link by {@code (provider, name)} in the given workspace.
     * Filters server-side using the Platform's keyword-search syntax
     * ({@code search=<name> provider:<provider>}) so the response contains at most
     * the matching data-link; this method returns the first result or {@code null}.
     *
     * Memoized per {@code (workspaceId, provider, name)} tuple, including {@code null}
     * misses — a path that repeatedly references a non-existent data-link does not
     * re-issue the search. Caveat: a data-link created on the Platform after a miss is
     * cached will not be visible until a new {@link SeqeraDataLinkClient} is constructed
     * (i.e. for the lifetime of a pipeline run, misses are sticky).
     */
    @Memoized
    DataLinkDto getDataLink(long workspaceId, String provider, String name) throws IOException {
        final search = "${name} provider:${provider}".toString()
        final it = PagedIterable.start(new DataLinkListFetcher(towerClient, endpoint, workspaceId, LIST_PAGE_SIZE, search)).iterator()
        return it.hasNext() ? it.next() : null
    }

    /**
     * Browse the content of a data-link.
     * The first page is fetched eagerly (so any IOException surfaces here, not at the
     * first iterator call); subsequent pages are fetched on demand as the returned
     * {@link PagedIterable} is iterated.
     *
     * Endpoints: {@code GET /data-links/{id}/browse} (root) and
     * {@code GET /data-links/{id}/browse/{path}} (sub-path).
     *
     * @param credentialsId optional data-link credentials identifier (from
     *     {@code DataLinkDto.credentials[0].id}); forwarded as a query param when set.
     * @param search optional server-side prefix filter on entry names
     */
    PagedIterable<DataLinkItem> getContent(String dataLinkId, String subPath, long workspaceId, String credentialsId = null, String search = null) throws IOException {
        log.debug("Getting content for data-link: $dataLinkId, path: $subPath, workspace: $workspaceId, credentialsId: $credentialsId")
        final pathSegment = subPath ? '/' + encodePath(subPath) : ''
        final baseUrl = "${endpoint}/data-links/${encodePath(dataLinkId)}/browse${pathSegment}"
        return PagedIterable.start(new DataLinkContentFetcher(towerClient, baseUrl, workspaceId, credentialsId, search))
    }

    /** {@code GET /data-links/{id}/generate-download-url?workspaceId=<ws>&filePath=<sub>[&credentialsId=<c>]} */
    DataLinkDownloadUrlResponse getDownloadUrl(String dataLinkId, String subPath, long workspaceId, String credentialsId = null) {
        String url = "${endpoint}/data-links/${encodePath(dataLinkId)}/generate-download-url?workspaceId=${workspaceId}&filePath=${encodeQuery(subPath ?: '')}"
        if (credentialsId)
            url += "&credentialsId=${encodeQuery(credentialsId)}"
        log.debug "Getting downloadURL: GET $url"
        final resp = towerClient.sendApiRequest(url)
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        final out = new DataLinkDownloadUrlResponse()
        out.url = json.url as String
        return out
    }

    // ---- page fetchers ----

    /**
     * {@link io.seqera.tower.plugin.datalink.PagedIterable.NextPageFetcher} for the
     * {@code /data-links} list endpoint (offset pagination).
     * Cursor state (offset + server-reported total) lives in instance fields.
     */
    @CompileStatic
    private static class DataLinkListFetcher implements PagedIterable.NextPageFetcher<DataLinkDto> {
        private final TowerClient towerClient
        private final String endpoint
        private final long workspaceId
        private final int pageSize
        private final String search

        private int offset = 0
        private long total = -1L   // unknown until the server reports totalSize

        DataLinkListFetcher(TowerClient towerClient, String endpoint, long workspaceId, int pageSize, String search) {
            this.towerClient = towerClient
            this.endpoint = endpoint
            this.workspaceId = workspaceId
            this.pageSize = pageSize
            this.search = search
        }

        @Override
        PagedIterable.Page<DataLinkDto> fetch() throws IOException {
            final url = "${endpoint}/data-links?workspaceId=${workspaceId}&max=${pageSize}&offset=${offset}${search ? '&search=' + encodeQuery(search) : ''}"
            log.debug "Fetching next list of datalinks: GET $url"
            final resp = towerClient.sendApiRequest(url)
            checkFsResponse(resp, url)
            final json = new JsonSlurper().parseText(resp.message) as Map
            final items = (json.dataLinks as List<Map>)?.collect { Map m -> mapDataLink(m) } ?: Collections.<DataLinkDto>emptyList()
            offset += items.size()
            if (total < 0 && json.totalSize != null)
                total = json.totalSize as Long
            final isLast = items.isEmpty() || (total >= 0 && offset >= total)
            return new PagedIterable.Page<DataLinkDto>(items, isLast)
        }
    }

    /**
     * {@link io.seqera.tower.plugin.datalink.PagedIterable.NextPageFetcher} for a
     * data-link's {@code /browse[/path]} endpoint (token pagination).
     * The next-page cursor lives in the {@code nextToken} instance field.
     */
    @CompileStatic
    private static class DataLinkContentFetcher implements PagedIterable.NextPageFetcher<DataLinkItem> {
        private final TowerClient towerClient
        private final String baseUrl
        private final long workspaceId
        private final String credentialsId
        private final String search

        private String nextToken = null

        DataLinkContentFetcher(TowerClient towerClient, String baseUrl, long workspaceId, String credentialsId, String search) {
            this.towerClient = towerClient
            this.baseUrl = baseUrl
            this.workspaceId = workspaceId
            this.credentialsId = credentialsId
            this.search = search
        }

        @Override
        PagedIterable.Page<DataLinkItem> fetch() throws IOException {
            String url = "${baseUrl}?workspaceId=${workspaceId}"
            if (search)
                url += "&search=${encodeQuery(search)}"
            if (credentialsId)
                url += "&credentialsId=${encodeQuery(credentialsId)}"
            if (nextToken)
                url += "&nextPageToken=${encodeQuery(nextToken)}"
            log.debug "Fetching Browse page GET $url"
            final resp = towerClient.sendApiRequest(url)
            checkFsResponse(resp, url)
            final json = new JsonSlurper().parseText(resp.message) as Map
            final items = (json.objects as List<Map>)?.collect { Map m -> mapItem(m) } ?: Collections.<DataLinkItem>emptyList()
            nextToken = json.nextPageToken as String
            return new PagedIterable.Page<DataLinkItem>(items, nextToken == null)
        }
    }

    // ---- encoding / error mapping ----

    /** URL-encode a path value while preserving {@code /} as path separators. */
    private static String encodePath(String s) {
        new URI(null, null, s ?: '', null).rawPath ?: ''
    }

    /** URL-encode a value intended for use as a query-string value. */
    private static String encodeQuery(String s) {
        URLEncoder.encode(s ?: '', 'UTF-8')
    }

    private static void checkFsResponse(TowerClient.Response resp, String url) {
        if (!resp.error)
            return
        final code = resp.code
        if (code == 401)
            throw new AbortOperationException("Seqera authentication failed — check tower.accessToken or TOWER_ACCESS_TOKEN")
        if (code == 403)
            throw new AccessDeniedException(url, null, "Forbidden — check workspace permissions")
        if (code == 404)
            throw new NoSuchFileException(url)
        throw new IOException("Seqera API error: HTTP ${code} for ${url}")
    }

    private static DataLinkDto mapDataLink(Map m) {
        final dto = new DataLinkDto()
        dto.id = m.id as String
        dto.name = m.name as String
        dto.description = m.description as String
        dto.resourceRef = m.resourceRef as String
        if (m.provider)
            dto.provider = parseProvider(m.provider as String)
        dto.region = m.region as String
        final credList = m.credentials as List<Map>
        if (credList)
            dto.credentials = credList.collect { Map c -> mapCredentials(c) }
        return dto
    }

    private static DataLinkCredentials mapCredentials(Map m) {
        final c = new DataLinkCredentials()
        c.id = m.id as String
        c.name = m.name as String
        if (m.provider)
            c.provider = parseProvider(m.provider as String)
        return c
    }

    private static DataLinkItem mapItem(Map m) {
        final it = new DataLinkItem()
        it.name = m.name as String
        if (m.type)
            it.type = parseItemType(m.type as String)
        it.size = (m.size as Long) ?: 0L
        it.mimeType = m.mimeType as String
        return it
    }

    private static DataLinkProvider parseProvider(String value) {
        try {
            return DataLinkProvider.fromValue(value)
        } catch (Throwable ignored) {
            return DataLinkProvider.values().find { DataLinkProvider p -> p.toString() == value }
        }
    }

    private static DataLinkItemType parseItemType(String value) {
        try {
            return DataLinkItemType.fromValue(value)
        } catch (Throwable ignored) {
            return DataLinkItemType.values().find { DataLinkItemType t -> t.toString() == value }
        }
    }
}
