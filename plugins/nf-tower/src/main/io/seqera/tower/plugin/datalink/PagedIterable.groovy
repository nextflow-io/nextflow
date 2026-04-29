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

/**
 * Generic lazy paginated view over a sequence of {@code T}.
 *
 * The first page is fetched eagerly (via {@link #start}) so any {@link IOException} from
 * the underlying source surfaces at the construction call site, not later at the first
 * {@link Iterator#hasNext()} call. Subsequent pages are fetched on demand as the
 * iterator advances.
 *
 * Pagination cursor state is captured by the {@link NextPage} implementation (closure
 * variables / fields) — this base class only knows "fetch returns more" vs "exhausted".
 */
@CompileStatic
class PagedIterable<T> implements Iterable<T> {

    /**
     * Stateful "give me the next page" callback. Implementations track their own
     * cursor (offset, token, etc.) across invocations.
     */
    static interface NextPageFetcher<T> {
        /** @return next page, never {@code null} (use empty items + {@code isLast=true} for end) */
        Page<T> fetch() throws IOException
    }

    @CompileStatic
    static class Page<T> {
        final List<T> items
        final boolean isLast
        Page(List<T> items, boolean isLast) {
            this.items = items ?: Collections.<T>emptyList()
            this.isLast = isLast
        }
    }

    protected final List<T> firstPage
    protected final boolean firstPageIsLast
    protected final NextPageFetcher<T> fetcher

    PagedIterable(List<T> firstPage, boolean firstPageIsLast, NextPageFetcher<T> fetcher) {
        this.firstPage = firstPage ?: Collections.<T>emptyList()
        this.firstPageIsLast = firstPageIsLast
        this.fetcher = fetcher
    }

    /** Eagerly fetch the first page; later pages on demand. Throws IOException at the call site on failure. */
    static <T> PagedIterable<T> start(NextPageFetcher<T> fetcher) throws IOException {
        final p = fetcher.fetch()
        if (p == null) return new PagedIterable<T>(Collections.<T>emptyList(), true, fetcher)
        return new PagedIterable<T>(p.items, p.isLast, fetcher)
    }

    /** First page (eagerly loaded). */
    List<T> getFirstPage() { Collections.unmodifiableList(firstPage) }

    boolean isEmpty() { firstPage.isEmpty() && firstPageIsLast }

    @Override
    Iterator<T> iterator() {
        return new PagedIterator()
    }

    /**
     * Lazy iterator that yields first-page items, then advances pages on demand.
     * Fetch failures are wrapped in {@link UncheckedIOException} (the {@link Iterator}
     * contract does not declare {@link IOException}).
     */
    @CompileStatic
    private class PagedIterator implements Iterator<T> {
        private Iterator<T> current = firstPage.iterator()
        private boolean exhausted = firstPageIsLast

        @Override
        boolean hasNext() {
            while (!current.hasNext()) {
                if (exhausted) return false
                try {
                    final p = fetcher.fetch()
                    final items = p?.items ?: Collections.<T>emptyList()
                    current = items.iterator()
                    if (p == null || p.isLast) exhausted = true
                } catch (IOException e) {
                    throw new UncheckedIOException(e)
                }
            }
            return true
        }

        @Override
        T next() {
            if (!hasNext()) throw new NoSuchElementException()
            return current.next()
        }
    }
}
