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

import groovy.transform.CompileStatic
import io.seqera.tower.model.DataLinkItem

/**
 * Lazy, paginated view over a data-link's content.
 *
 * The first page is fetched eagerly by the producer so callers can inspect
 * {@link #getOriginalPath()} and {@link #getFirstPage()} without triggering
 * additional HTTP calls. Iterating yields items from the first page followed
 * by subsequent pages fetched on demand via the injected page fetcher.
 */
@CompileStatic
class PagedDataLinkContent implements Iterable<DataLinkItem> {

    /**
     * Opaque page fetcher. Given a {@code nextPageToken}, returns the next page
     * as a map with keys {@code objects} ({@code List<DataLinkItem>}) and
     * {@code nextPageToken} ({@code String}, null if no more pages).
     */
    static interface PageFetcher {
        Map<String, Object> fetch(String nextPageToken) throws IOException
    }

    private final String originalPath
    private final List<DataLinkItem> firstPage
    private final String firstPageNextToken
    private final PageFetcher pageFetcher

    PagedDataLinkContent(String originalPath,
                         List<DataLinkItem> firstPage,
                         String firstPageNextToken,
                         PageFetcher pageFetcher) {
        this.originalPath = originalPath
        this.firstPage = firstPage ?: Collections.<DataLinkItem>emptyList()
        this.firstPageNextToken = firstPageNextToken
        this.pageFetcher = pageFetcher
    }

    String getOriginalPath() { originalPath }

    /** First page, loaded eagerly — bounded in size by the server's page size. */
    List<DataLinkItem> getFirstPage() { Collections.unmodifiableList(firstPage) }

    boolean isEmpty() { firstPage.isEmpty() && !firstPageNextToken }

    @Override
    Iterator<DataLinkItem> iterator() {
        return new PagedIterator(firstPage, firstPageNextToken, pageFetcher)
    }

    /** Lazy iterator that paginates on demand. */
    @CompileStatic
    private static class PagedIterator implements Iterator<DataLinkItem> {
        private Iterator<DataLinkItem> current
        private String nextToken
        private final PageFetcher fetcher

        PagedIterator(List<DataLinkItem> firstPage, String firstPageNextToken, PageFetcher fetcher) {
            this.current = firstPage.iterator()
            this.nextToken = firstPageNextToken
            this.fetcher = fetcher
        }

        @Override
        boolean hasNext() {
            while (!current.hasNext()) {
                if (!nextToken) return false
                try {
                    final page = fetcher.fetch(nextToken)
                    final items = (page?.objects ?: []) as List<DataLinkItem>
                    current = items.iterator()
                    nextToken = page?.nextPageToken as String
                } catch (IOException e) {
                    throw new RuntimeException(e)
                }
            }
            return true
        }

        @Override
        DataLinkItem next() {
            if (!hasNext()) throw new NoSuchElementException()
            return current.next()
        }
    }
}
