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

package io.seqera.tower.plugin.fs

import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.WatchService
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.OrgAndWorkspaceDto
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient

/**
 * FileSystem instance for the {@code seqera://} scheme.
 * One instance per {@link SeqeraFileSystemProvider}.
 *
 * Resource-type-agnostic: the filesystem owns the org/workspace cache (shared across
 * resource types) and a registry of {@link ResourceTypeHandler}s. Each handler owns
 * its own API client and resource-specific caches.
 */
@Slf4j
@CompileStatic
class SeqeraFileSystem extends FileSystem {

    private final SeqeraFileSystemProvider provider0
    private SeqeraDatasetClient orgWorkspaceClient

    /** orgName → orgId */
    private final Map<String, Long> orgCache = new LinkedHashMap<>()
    /** "orgName/workspaceName" → workspaceId */
    private final Map<String, Long> workspaceCache = new LinkedHashMap<>()

    /** resourceType → handler */
    private final Map<String, ResourceTypeHandler> handlers = new LinkedHashMap<>()

    private volatile boolean orgWorkspaceCacheLoaded = false

    SeqeraFileSystem(SeqeraFileSystemProvider provider) {
        this.provider0 = provider
    }

    /**
     * Attach the dataset client used for user-info / workspaces lookup. The org/workspace
     * listing uses dataset endpoints today ({@code GET /user-info}, {@code GET /user/{id}/workspaces});
     * keeping the client on the filesystem avoids duplicating it across handlers.
     */
    void setOrgWorkspaceClient(SeqeraDatasetClient client) {
        this.orgWorkspaceClient = client
    }

    @Override
    FileSystemProvider provider() { provider0 }

    @Override
    void close() { /* no-op */ }

    @Override
    boolean isOpen() { true }

    @Override
    boolean isReadOnly() { true }

    @Override
    String getSeparator() { '/' }

    @Override
    Iterable<Path> getRootDirectories() {
        return [getPath('seqera://')] as Iterable<Path>
    }

    @Override
    Iterable<FileStore> getFileStores() { Collections.emptyList() }

    @Override
    Set<String> supportedFileAttributeViews() { Collections.singleton('basic') }

    @Override
    Path getPath(String first, String... more) {
        final full = more ? ([first] + more.toList()).join(getSeparator()) : first
        return new SeqeraPath(this, full)
    }

    @Override
    PathMatcher getPathMatcher(String syntaxAndPattern) {
        throw new UnsupportedOperationException("PathMatcher not supported by seqera:// filesystem")
    }

    @Override
    UserPrincipalLookupService getUserPrincipalLookupService() {
        throw new UnsupportedOperationException("UserPrincipalLookupService not supported by seqera:// filesystem")
    }

    @Override
    WatchService newWatchService() {
        throw new UnsupportedOperationException("WatchService not supported by seqera:// filesystem")
    }

    // ---- org/workspace cache (shared across handlers) ----

    synchronized void loadOrgWorkspaceCache() {
        if (orgWorkspaceCacheLoaded) return
        if (!orgWorkspaceClient)
            throw new IllegalStateException("SeqeraFileSystem has no orgWorkspaceClient attached")
        log.debug "Loading Seqera org/workspace cache"
        final entries = orgWorkspaceClient.listUserWorkspacesAndOrgs(orgWorkspaceClient.getUserId())
        for (OrgAndWorkspaceDto entry : entries) {
            if (entry.orgName)
                orgCache.put(entry.orgName, entry.orgId)
            if (entry.orgName && entry.workspaceName && entry.workspaceId)
                workspaceCache.put("${entry.orgName}/${entry.workspaceName}" as String, entry.workspaceId)
        }
        orgWorkspaceCacheLoaded = true
    }

    synchronized Set<String> listOrgNames() {
        loadOrgWorkspaceCache()
        return Collections.unmodifiableSet(orgCache.keySet())
    }

    synchronized List<String> listWorkspaceNames(String org) {
        loadOrgWorkspaceCache()
        return workspaceCache.keySet()
                .findAll { String k -> k.startsWith("${org}/") }
                .collect { String k -> k.substring(org.length() + 1) }
    }

    synchronized long resolveWorkspaceId(String org, String workspace) throws NoSuchFileException {
        loadOrgWorkspaceCache()
        final key = "${org}/${workspace}" as String
        final id = workspaceCache.get(key)
        if (id == null)
            throw new NoSuchFileException("seqera://${key}", null, "Org or workspace not found or not accessible")
        return id
    }

    // ---- handler registry ----

    synchronized void registerHandler(ResourceTypeHandler handler) {
        handlers.put(handler.resourceType, handler)
    }

    synchronized ResourceTypeHandler getHandler(String resourceType) {
        handlers.get(resourceType)
    }

    synchronized Set<String> getResourceTypes() {
        Collections.unmodifiableSet(new LinkedHashSet<String>(handlers.keySet()))
    }
}
