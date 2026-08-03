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

import spock.lang.Specification

/**
 * Tests for {@link PagedIterable} — eager-first-page, lazy subsequent pages, and the
 * single-use {@link PagedIterable#iterator()} contract.
 */
class PagedIterableTest extends Specification {

    /** Fetcher that hands out the given pages in order, then reports exhaustion. */
    private static class ListFetcher<T> implements PagedIterable.NextPageFetcher<T> {
        private final List<List<T>> pages
        int fetchCount = 0
        ListFetcher(List<List<T>> pages) { this.pages = pages }
        @Override
        PagedIterable.Page<T> fetch() throws IOException {
            final idx = fetchCount++
            if (idx >= pages.size())
                return new PagedIterable.Page<T>([], true)
            final items = pages[idx]
            return new PagedIterable.Page<T>(items, idx == pages.size() - 1)
        }
    }

    private static <T> List<T> drain(Iterator<T> it) {
        final out = new ArrayList<T>()
        while (it.hasNext()) out.add(it.next())
        return out
    }

    def "yields items from a single page"() {
        given:
        def fetcher = new ListFetcher<String>([['a', 'b', 'c']])

        when:
        def paged = PagedIterable.start(fetcher)

        then:
        paged.firstPage == ['a', 'b', 'c']
        !paged.isEmpty()

        and:
        drain(paged.iterator()) == ['a', 'b', 'c']
    }

    def "iterates lazily across multiple pages"() {
        given:
        def fetcher = new ListFetcher<String>([['a'], ['b'], ['c']])

        when:
        def paged = PagedIterable.start(fetcher)

        then: 'only the first page has been fetched eagerly'
        fetcher.fetchCount == 1

        when:
        def all = drain(paged.iterator())

        then:
        all == ['a', 'b', 'c']
    }

    def "does not fetch page 2 if only the first page is consumed"() {
        given:
        def fetcher = new ListFetcher<String>([['a'], ['b']])

        when:
        def paged = PagedIterable.start(fetcher)
        def first = paged.iterator().next()

        then:
        first == 'a'
        fetcher.fetchCount == 1
    }

    def "reports empty when the first page is empty and last"() {
        given:
        def fetcher = new ListFetcher<String>([[]])

        when:
        def paged = PagedIterable.start(fetcher)

        then:
        paged.isEmpty()
        !paged.iterator().hasNext()
    }

    def "surfaces IOException eagerly at the start() call site"() {
        given:
        def fetcher = Stub(PagedIterable.NextPageFetcher) {
            fetch() >> { throw new IOException("boom") }
        }

        when:
        PagedIterable.start(fetcher)

        then:
        thrown(IOException)
    }

    def "wraps a later-page IOException as UncheckedIOException during iteration"() {
        given: 'first page ok, second page throws'
        int n = 0
        def fetcher = new PagedIterable.NextPageFetcher<String>() {
            @Override
            PagedIterable.Page<String> fetch() throws IOException {
                if (n++ == 0) return new PagedIterable.Page<String>(['a'], false)
                throw new IOException("page 2 failed")
            }
        }
        def paged = PagedIterable.start(fetcher)

        when:
        drain(paged.iterator())

        then:
        thrown(UncheckedIOException)
    }

    def "getFirstPage is immutable"() {
        given:
        def paged = PagedIterable.start(new ListFetcher<String>([['a']]))

        when:
        paged.firstPage.add('b')

        then:
        thrown(UnsupportedOperationException)
    }

    def "iterator() is single-use — a second call throws IllegalStateException"() {
        given:
        def paged = PagedIterable.start(new ListFetcher<String>([['a', 'b']]))

        when: 'first iterator is fine'
        paged.iterator()

        and: 'second iterator is rejected'
        paged.iterator()

        then:
        thrown(IllegalStateException)
    }

    def "getFirstPage does not consume the single iterator"() {
        given:
        def paged = PagedIterable.start(new ListFetcher<String>([['a', 'b']]))

        when: 'reading firstPage does not count as iterating'
        paged.firstPage
        def items = drain(paged.iterator())

        then:
        items == ['a', 'b']
    }
}