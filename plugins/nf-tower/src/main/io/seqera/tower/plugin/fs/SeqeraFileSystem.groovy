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
import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.OrgAndWorkspaceDto
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient

/**
 * FileSystem instance for the {@code seqera://} scheme.
 * One instance per (endpoint + credentials) pair, cached by {@link SeqeraFileSystemProvider}.
 *
 * Lazily populates org/workspace/dataset caches on first access.
 * Cache is invalidated on dataset write operations.
 *
 * @author Seqera Labs
 */
@Slf4j
@CompileStatic
class SeqeraFileSystem extends FileSystem {

    private final SeqeraFileSystemProvider provider0
    final SeqeraDatasetClient client

    /** orgName → orgId */
    private final Map<String, Long> orgCache = new LinkedHashMap<>()
    /** "orgName/workspaceName" → workspaceId */
    private final Map<String, Long> workspaceCache = new LinkedHashMap<>()
    /** workspaceId → list of DatasetDto */
    private final Map<Long, List<DatasetDto>> datasetCache = new LinkedHashMap<>()

    private volatile boolean orgWorkspaceCacheLoaded = false

    SeqeraFileSystem(SeqeraFileSystemProvider provider, SeqeraDatasetClient client) {
        this.provider0 = provider
        this.client = client
    }

    @Override
    FileSystemProvider provider() { provider0 }

    @Override
    void close() { /* no-op: platform API connection is stateless */ }

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

    // ---- cache management ----

    /**
     * Ensure the org/workspace cache is populated. Thread-safe: loads at most once.
     * Calls GET /user-info then GET /user/{userId}/workspaces.
     */
    synchronized void loadOrgWorkspaceCache() {
        if (orgWorkspaceCacheLoaded) return
        log.debug "Loading Seqera org/workspace cache"
        final userInfo = client.getUserInfo()
        final userId = (userInfo.id as long)
        final entries = client.listUserWorkspacesAndOrgs(userId)
        for (OrgAndWorkspaceDto entry : entries) {
            if (entry.orgName)
                orgCache.put(entry.orgName, entry.orgId)
            if (entry.orgName && entry.workspaceName && entry.workspaceId)
                workspaceCache.put("${entry.orgName}/${entry.workspaceName}" as String, entry.workspaceId)
        }
        orgWorkspaceCacheLoaded = true
    }

    /**
     * @return distinct org names visible to the authenticated user
     */
    synchronized Set<String> listOrgNames() {
        loadOrgWorkspaceCache()
        return Collections.unmodifiableSet(orgCache.keySet())
    }

    /**
     * @return workspace names for the given org
     */
    synchronized List<String> listWorkspaceNames(String org) {
        loadOrgWorkspaceCache()
        return workspaceCache.keySet()
            .findAll { String k -> k.startsWith("${org}/") }
            .collect { String k -> k.substring(org.length() + 1) }
    }

    /**
     * Resolve a workspace ID by org and workspace name.
     * @throws NoSuchFileException if the org or workspace is not in the cache
     */
    synchronized long resolveWorkspaceId(String org, String workspace) throws NoSuchFileException {
        loadOrgWorkspaceCache()
        final key = "${org}/${workspace}" as String
        final id = workspaceCache.get(key)
        if (id == null)
            throw new NoSuchFileException("seqera://${key}", null, "Org or workspace not found or not accessible")
        return id
    }

    /**
     * Return datasets for the given workspace, populating the cache on first access.
     */
    synchronized List<DatasetDto> resolveDatasets(long workspaceId) {
        List<DatasetDto> cached = datasetCache.get(workspaceId)
        if (cached == null) {
            cached = client.listDatasets(workspaceId)
            datasetCache.put(workspaceId, cached)
        }
        return cached
    }

    /**
     * Invalidate the dataset cache for a workspace (call after a write operation).
     */
    synchronized void invalidateDatasetCache(long workspaceId) {
        datasetCache.remove(workspaceId)
    }

    /**
     * Resolve a DatasetDto by name within a workspace.
     * @throws NoSuchFileException if no dataset with the given name exists
     */
    DatasetDto resolveDataset(long workspaceId, String name) throws NoSuchFileException {
        final datasets = resolveDatasets(workspaceId)
        final found = datasets.find { DatasetDto d -> d.name == name }
        if (!found)
            throw new NoSuchFileException("seqera:///datasets/${name}", null, "Dataset '${name}' not found in workspace")
        return found
    }
}
