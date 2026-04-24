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

import java.nio.file.Path

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
     * Pages are fetched from {@code GET /data-links?workspaceId=<ws>&max=<n>&offset=<o>}
     * on demand as the iterator advances.
     */
    Iterator<DataLinkDto> listDataLinks(long workspaceId) {
        return new DataLinkListIterator(towerClient, endpoint, workspaceId, LIST_PAGE_SIZE)
    }

    /**
     * Resolve a data-link providers in the given workspace.
     */
    @Memoized
    Set<String> getDataLinkProviders(long workspaceId) {
        final providers = new TreeSet<String>()
        final Iterator<DataLinkDto> it = listDataLinks(workspaceId)
        while (it.hasNext()) {
            final p = it.next().provider?.toString()
            if (p) providers.add(p)
        }
        return providers
    }

    /**
     * Resolve a data-link by {@code (provider, name)} in the given workspace.
     * Iterates the API's list endpoint lazily and short-circuits on first match.
     */
    @Memoized
    DataLinkDto getDataLink(long workspaceId, String provider, String name) {
        final Iterator<DataLinkDto> it = new DataLinkListIterator(towerClient, endpoint, workspaceId, LIST_PAGE_SIZE, name)
        while( it.hasNext() ) {
            final dl = it.next()
            if( dl.provider?.toString() == provider )
                return dl
        }
        throw new NoSuchFileException(
            "seqera://.../data-links/${provider}/${name}",
            null,
            "Data-link '${name}' not found for provider '${provider}' in workspace '$workspaceId'")
    }

    /**
     * Browse the content of a data-link.
     * The first page is fetched eagerly to populate metadata ({@code originalPath},
     * first-page items). Subsequent pages are fetched on demand as the returned
     * {@link PagedDataLinkContent} is iterated.
     *
     * Endpoints: {@code GET /data-links/{id}/browse} (root) and
     * {@code GET /data-links/{id}/browse/{path}} (sub-path).
     *
     * @param credentialsId optional data-link credentials identifier (from
     *     {@code DataLinkDto.credentials[0].id}); forwarded as a query param when set.
     */
    PagedDataLinkContent getContent(String dataLinkId, String subPath, long workspaceId, String credentialsId = null) {
        final pathSegment = subPath ? '/' + encodePath(subPath) : ''
        final baseUrl = "${endpoint}/data-links/${encodePath(dataLinkId)}/browse${pathSegment}"
        final page = fetchBrowsePage(baseUrl, workspaceId, credentialsId, null)
        final firstItems = page.objects
        final firstToken = page.nextPageToken
        final originalPath = page.originalPath
        final fetcher = new PagedDataLinkContent.PageFetcher() {
            @Override
            Map<String, Object> fetch(String token) throws IOException {
                final next = fetchBrowsePage(baseUrl, workspaceId, credentialsId, token)
                return [objects: next.objects, nextPageToken: next.nextPageToken] as Map<String, Object>
            }
        }
        return new PagedDataLinkContent(originalPath, firstItems, firstToken, fetcher)
    }

    /** {@code GET /data-links/{id}/generate-download-url?workspaceId=<ws>&filePath=<sub>[&credentialsId=<c>]} */
    DataLinkDownloadUrlResponse getDownloadUrl(String dataLinkId, String subPath, long workspaceId, String credentialsId = null) {
        String url = "${endpoint}/data-links/${encodePath(dataLinkId)}/generate-download-url?workspaceId=${workspaceId}&filePath=${encodeQuery(subPath ?: '')}"
        if (credentialsId) url += "&credentialsId=${encodeQuery(credentialsId)}"
        log.debug "Getting downloadURL: GET $url"
        final resp = towerClient.sendApiRequest(url)
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        final out = new DataLinkDownloadUrlResponse()
        out.url = json.url as String
        return out
    }

    // ---- page-fetching helpers ----

    /** Fetch one browse page and normalize it into a {@link BrowsePage}. */
    private BrowsePage fetchBrowsePage(String baseUrl, long workspaceId, String credentialsId, String nextPageToken) {
        String url = "${baseUrl}?workspaceId=${workspaceId}"
        if (credentialsId) url += "&credentialsId=${encodeQuery(credentialsId)}"
        if (nextPageToken) url += "&nextPageToken=${encodeQuery(nextPageToken)}"
        log.debug "Fetching Browse page GET $url"
        final resp = towerClient.sendApiRequest(url)
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        final items = (json.objects as List<Map>)?.collect { Map m -> mapItem(m) } ?: Collections.<DataLinkItem>emptyList()
        return new BrowsePage(json.originalPath as String, items, json.nextPageToken as String)
    }

    @CompileStatic
    private static class BrowsePage {
        final String originalPath
        final List<DataLinkItem> objects
        final String nextPageToken

        BrowsePage(String originalPath, List<DataLinkItem> objects, String nextPageToken) {
            this.originalPath = originalPath
            this.objects = objects
            this.nextPageToken = nextPageToken
        }
    }

    /** Lazy iterator for the {@code /data-links} list endpoint (offset pagination). */
    @CompileStatic
    private static class DataLinkListIterator implements Iterator<DataLinkDto> {
        private final TowerClient towerClient
        private final String endpoint
        private final long workspaceId
        private final int pageSize
        private final String search

        private Iterator<DataLinkDto> current = Collections.<DataLinkDto>emptyIterator()
        private int offset = 0
        private long total = -1L   // unknown until first fetch
        private boolean exhausted = false

        DataLinkListIterator(TowerClient towerClient, String endpoint, long workspaceId, int pageSize, String search = null) {
            this.towerClient = towerClient
            this.endpoint = endpoint
            this.workspaceId = workspaceId
            this.pageSize = pageSize
            this.search = search
        }

        @Override
        boolean hasNext() {
            while (!current.hasNext()) {
                if (exhausted) return false
                fetchNextPage()
            }
            return true
        }

        @Override
        DataLinkDto next() {
            if (!hasNext()) throw new NoSuchElementException()
            return current.next()
        }

        private void fetchNextPage() {
            final url = "${endpoint}/data-links?workspaceId=${workspaceId}&max=${pageSize}&offset=${offset}${search ? '&search='+ encodeQuery(search) :''}"
            log.debug "Fetching next list of datalinks: GET $url"
            final resp = towerClient.sendApiRequest(url)
            checkFsResponse(resp, url)
            final json = new JsonSlurper().parseText(resp.message) as Map
            final items = (json.dataLinks as List<Map>)?.collect { Map m -> mapDataLink(m) } ?: Collections.<DataLinkDto>emptyList()
            current = items.iterator()
            offset += items.size()
            if (total < 0) total = (json.totalSize as Long) ?: 0L
            if (items.isEmpty() || offset >= total) exhausted = true
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
        if (!resp.error) return
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
        if (m.provider) dto.provider = parseProvider(m.provider as String)
        dto.region = m.region as String
        final credList = m.credentials as List<Map>
        if (credList) dto.credentials = credList.collect { Map c -> mapCredentials(c) }
        return dto
    }

    private static DataLinkCredentials mapCredentials(Map m) {
        final c = new DataLinkCredentials()
        c.id = m.id as String
        c.name = m.name as String
        if (m.provider) c.provider = parseProvider(m.provider as String)
        return c
    }

    private static DataLinkItem mapItem(Map m) {
        final it = new DataLinkItem()
        it.name = m.name as String
        if (m.type) it.type = parseItemType(m.type as String)
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
