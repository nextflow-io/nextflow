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
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.time.Duration
import java.time.Instant

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.plugin.datalink.PagedIterable
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.ResourceTypeHandler
import io.seqera.tower.plugin.fs.SeqeraFileAttributes
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath

/**
 * {@link ResourceTypeHandler} for the {@code data-links} resource type.
 *
 * Listings and attribute queries go through the Seqera Platform API; file reads
 * use a pre-signed URL obtained from {@code /generate-download-url} and fetched
 * with a plain JDK {@link HttpClient} — the Seqera {@code Authorization} header
 * must not be sent to the cloud-backed URL.
 *
 * Data-link list and directory content are streamed lazily to avoid materializing
 * potentially large result sets in memory.
 */
@Slf4j
@CompileStatic
class DataLinksResourceHandler implements ResourceTypeHandler {

    public static final String TYPE = 'data-links'

    private final SeqeraFileSystem fs
    private final SeqeraDataLinkClient client
    private final HttpClient httpClient

    DataLinksResourceHandler(SeqeraFileSystem fs, SeqeraDataLinkClient client) {
        this(fs, client, HttpClient.newBuilder()
                .connectTimeout(Duration.ofSeconds(10))
                .followRedirects(HttpClient.Redirect.NORMAL)
                .build())
    }

    /** Test-only constructor to inject a mock {@link HttpClient}. */
    DataLinksResourceHandler(SeqeraFileSystem fs, SeqeraDataLinkClient client, HttpClient httpClient) {
        this.fs = fs
        this.client = client
        this.httpClient = httpClient
    }

    @Override
    String getResourceType() { TYPE }

    @Override
    Iterable<Path> list(SeqeraPath dir) throws IOException {
        final workspaceId = fs.resolveWorkspaceId(dir.org, dir.workspace)
        final trail = dir.trail
        if (trail.isEmpty()) {
            // data-links/ → distinct providers in use (sorted). Iterate the stream,
            // collect distinct provider names — small output.
            final providers = client.getDataLinkProviders(workspaceId)
            return providers.collect { String p -> dir.resolve(p) as Path }
        }
        if (trail.size() == 1) {
            // data-links/<provider>/ → sorted data-link names for that provider
            final prov = trail[0]
            final names = new TreeSet<String>()
            final Iterator<DataLinkDto> it = client.listDataLinks(workspaceId)
            while (it.hasNext()) {
                final dl = it.next()
                if (dl.provider?.toString() == prov)
                    names.add(dl.name)
            }
            if (names.isEmpty())
                throw new NoSuchFileException(dir.toString(), null, "No data-links for provider '$prov' in workspace '${dir.workspace}'")
            return names.collect { String n -> dir.resolve(n) as Path }
        }
        // trail.size() >= 2 — browse inside a specific data-link.
        // Content can be very large, so we stream it lazily.
        final dl = requireDataLink(workspaceId, trail[0], trail[1], dir)
        final subPath = trail.size() > 2 ? trail.subList(2, trail.size()).join('/') : ''
        log.debug("Listing files for $dl.name path $subPath")
        final content = client.getContent(dl.id, subPath, workspaceId, credentialsIdOf(dl))
        return new PathMappingIterable(content, dir)
    }

    @Override
    SeqeraFileAttributes readAttributes(SeqeraPath p) throws IOException {
        // Short-circuit: attributes attached when this path was produced by a listing
        if (p.cachedAttributes)
            return p.cachedAttributes
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final trail = p.trail
        if (trail.isEmpty()) {
            // data-links/ — always a directory
            return new SeqeraFileAttributes(true)
        }
        if (trail.size() == 1) {
            // data-links/<provider> — validate the provider has at least one data-link
            final providers = client.getDataLinkProviders(workspaceId)
            if (!providers.contains(trail[0]))
                throw new NoSuchFileException(p.toString(), null, "No data-links for provider '${trail[0]}' in workspace '${p.workspace}'")
            return new SeqeraFileAttributes(true)
        }
        final dl = requireDataLink(workspaceId, trail[0], trail[1], p)
        if (trail.size() == 2)
            return new SeqeraFileAttributes(true) // data-link root
        final subPath = trail.subList(2, trail.size()).join('/')
        log.debug("Reading attributes for $p")
        return resolveAttrsViaParent(dl, subPath, workspaceId, p)
    }

    @Override
    InputStream newInputStream(SeqeraPath p) throws IOException {
        if (p.trail.size() < 3)
            throw new IllegalArgumentException("newInputStream requires a file path inside a data-link: $p")
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final dl = requireDataLink(workspaceId, p.trail[0], p.trail[1], p)
        final subPath = p.trail.subList(2, p.trail.size()).join('/')
        final urlResp = client.getDownloadUrl(dl.id, subPath, workspaceId, credentialsIdOf(dl))
        if (!urlResp.url)
            throw new NoSuchFileException(p.toString(), null, "Platform returned no download URL")
        return fetchSignedUrl(urlResp.url)
    }

    /** First credentials entry on the data-link (or null if none). */
    private static String credentialsIdOf(DataLinkDto dl) {
        final creds = dl?.credentials
        return (creds && !creds.isEmpty()) ? creds[0].id : null
    }

    /**
     * Resolve a data-link by (provider, name) and throw {@link NoSuchFileException}
     * with a uniform error message when missing. Wraps {@link SeqeraDataLinkClient#getDataLink}
     * (which returns {@code null} on miss).
     */
    private DataLinkDto requireDataLink(long workspaceId, String provider, String name, SeqeraPath pathForErrors) throws NoSuchFileException {
        final dl = client.getDataLink(workspaceId, provider, name)
        if (!dl)
            throw new NoSuchFileException(pathForErrors.toString(), null, "Data-link '${name}' not found for provider '${provider}' in workspace '$workspaceId'")
        return dl
    }

    // ---- private helpers ----

    /**
     * Determine attributes for a path inside a data-link by listing the path's parent
     * directory and finding the entry by name. The entry's {@code type} (FILE/FOLDER)
     * is the authoritative file-vs-directory signal; absence of the entry means the
     * path does not exist.
     *
     * The Platform's {@code /browse/{path}} response does not reliably distinguish
     * "file path", "directory path", and "missing path" by itself ({@code originalPath}
     * is populated in all three), so we always query the parent.
     *
     * Cost: one API call (or more if the parent listing paginates and the entry isn't
     * on the first page). Iteration is short-circuited as soon as the entry is found.
     */
    private SeqeraFileAttributes resolveAttrsViaParent(DataLinkDto dl, String subPath, long workspaceId, SeqeraPath pathForErrors) throws IOException {
        final lastSlash = subPath.lastIndexOf('/')
        final parentSubPath = lastSlash > 0 ? subPath.substring(0, lastSlash) : ''
        final lastSeg = lastSlash > 0 ? subPath.substring(lastSlash + 1) : subPath
        log.debug("Looking for $lastSeg in data-link: ${dl.id} path: $parentSubPath")
        final parent = client.getContent(dl.id, parentSubPath, workspaceId, credentialsIdOf(dl), lastSeg)
        DataLinkItem found = null
        for (DataLinkItem it : parent) {
            log.trace("Item: $it")
            if (it.name == lastSeg || it.name == lastSeg + '/') {
                found = it
                break
            }
        }
        if (found == null)
            throw new NoSuchFileException(pathForErrors.toString(), null, "Path '${subPath}' not found inside data-link '${dl.name}'")
        if (found.type == DataLinkItemType.FILE) {
            final long size = (found.size != null) ? found.size.longValue() : 0L
            return new SeqeraFileAttributes(size, Instant.EPOCH, Instant.EPOCH, pathForErrors.toString())
        }
        return new SeqeraFileAttributes(true)
    }

    private InputStream fetchSignedUrl(String url) throws IOException {
        final req = HttpRequest.newBuilder(URI.create(url))
                .timeout(Duration.ofMinutes(5))
                .GET()
                .build()
        try {
            final HttpResponse<InputStream> resp = httpClient.send(req, HttpResponse.BodyHandlers.ofInputStream())
            final status = resp.statusCode()
            if (status >= 200 && status < 300)
                return resp.body()
            try {
                resp.body()?.close()
            }
            catch (Throwable ignored) {
                // best-effort close of the error-response body
            }
            throw new IOException("Signed URL fetch failed: HTTP $status for ${redactUrl(url)}")
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("Interrupted while fetching signed URL", e)
        }
    }

    /**
     * Strip the query string from a pre-signed URL before it is logged or surfaced in an
     * exception — the query carries temporary credentials/signatures (e.g. {@code X-Amz-Signature},
     * SAS tokens, GCS {@code Signature}) that must not leak into logs or error reports.
     */
    private static String redactUrl(String url) {
        try {
            final u = URI.create(url)
            return new URI(u.scheme, u.authority, u.path, null, null).toString()
        } catch (Exception ignored) {
            return '(redacted)'
        }
    }

    /**
     * Lazy {@link Iterable} that maps each {@link DataLinkItem} from a
     * {@link PagedIterable<DataLinkItem>} to a child {@link SeqeraPath} under
     * {@code parent}. Each produced path carries cached attributes built from the
     * item, so a follow-up {@code readAttributes()} call does not re-browse the
     * Platform. Pages are fetched on demand as the iterator advances.
     */
    @CompileStatic
    private static class PathMappingIterable implements Iterable<Path> {
        private final PagedIterable<DataLinkItem> content
        private final SeqeraPath parent

        PathMappingIterable(PagedIterable<DataLinkItem> content, SeqeraPath parent) {
            this.content = content
            this.parent = parent
        }

        @Override
        Iterator<Path> iterator() {
            final Iterator<DataLinkItem> inner = content.iterator()
            return new Iterator<Path>() {
                @Override
                boolean hasNext() { inner.hasNext() }

                @Override
                Path next() {
                    final item = inner.next()
                    return parent.resolveWithAttributes(item.name, attributesFor(item)) as Path
                }
            }
        }

        private static SeqeraFileAttributes attributesFor(DataLinkItem item) {
            if (item.type == DataLinkItemType.FILE)
                return new SeqeraFileAttributes(item.size ?: 0L, Instant.EPOCH, Instant.EPOCH, item.name)
            return new SeqeraFileAttributes(true)
        }
    }
}
