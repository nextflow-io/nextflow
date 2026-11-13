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

import java.nio.file.AccessDeniedException
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
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.exception.ForbiddenException
import io.seqera.tower.plugin.exception.NotFoundException
import io.seqera.tower.plugin.exception.UnauthorizedException
import nextflow.exception.AbortOperationException

/**
 * FileSystem instance for the {@code seqera://} scheme.
 * One instance per {@link SeqeraFileSystemProvider}.
 *
 * Resource-type-agnostic: the filesystem owns the user-id and org/workspace caches
 * (shared across resource types) and a registry of {@link ResourceTypeHandler}s.
 * Each handler owns its own API client and resource-specific caches.
 */
@Slf4j
@CompileStatic
class SeqeraFileSystem extends FileSystem {

    private final SeqeraFileSystemProvider provider0
    private final TowerClient towerClient

    /** Cached current-user id; the user is fixed by the {@code TowerClient}'s access token. */
    private volatile Long cachedUserId

    /** orgName → orgId */
    private final Map<String, Long> orgCache = new LinkedHashMap<>()
    /** "orgName/workspaceName" → workspaceId */
    private final Map<String, Long> workspaceCache = new LinkedHashMap<>()

    /** resourceType → handler */
    private final Map<String, ResourceTypeHandler> handlers = new LinkedHashMap<>()

    private volatile boolean orgWorkspaceCacheLoaded = false

    SeqeraFileSystem(SeqeraFileSystemProvider provider, TowerClient towerClient) {
        this.provider0 = provider
        this.towerClient = towerClient
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

    // ---- user-id / org / workspace caches (shared across handlers) ----

    /**
     * Resolve the current user's numeric ID via {@code GET /user-info}.
     * Cached for the lifetime of this filesystem — the token does not change
     * during a pipeline run, so neither does the resolved user.
     */
    synchronized long getUserId() throws IOException {
        if (cachedUserId != null) return cachedUserId
        try {
            final info = towerClient.getUserInfo()
            if (info?.id == null)
                throw new AbortOperationException("Unable to retrieve user ID from Seqera Platform — check your access token")
            cachedUserId = info.id as Long
            return cachedUserId
        } catch (UnauthorizedException e) {
            throw new AbortOperationException(e.getMessage())
        } catch (ForbiddenException e) {
            throw new AccessDeniedException("${towerClient.endpoint}/user-info", null, e.message)
        } catch (NotFoundException e) {
            throw new NoSuchFileException("${towerClient.endpoint}/user-info")
        }
    }

    /** {@code GET /user/{userId}/workspaces} — reachable orgs and workspaces. */
    private List<OrgAndWorkspaceDto> fetchUserWorkspacesAndOrgs(long userId) throws IOException {
        try {
            final list = towerClient.listUserWorkspacesAndOrgs(userId as String)
            return list.collect { Map m -> mapOrgAndWorkspace(m) }
        } catch (UnauthorizedException e) {
            throw new AbortOperationException(e.getMessage())
        } catch (ForbiddenException e) {
            throw new AccessDeniedException("${towerClient.endpoint}/user/$userId/workspaces", null, e.message)
        } catch (NotFoundException e) {
            throw new NoSuchFileException("${towerClient.endpoint}/user/$userId/workspaces")
        }
    }

    private static OrgAndWorkspaceDto mapOrgAndWorkspace(Map m) {
        final dto = new OrgAndWorkspaceDto()
        dto.orgId = (m.orgId as Long) ?: 0L
        dto.orgName = m.orgName as String
        dto.workspaceId = (m.workspaceId as Long) ?: 0L
        dto.workspaceName = m.workspaceName as String
        dto.workspaceFullName = m.workspaceFullName as String
        return dto
    }

    /**
     * Ensure the org/workspace cache is populated. Thread-safe: loads at most once.
     * Calls GET /user-info then GET /user/{userId}/workspaces.
     */
    synchronized void loadOrgWorkspaceCache() {
        if (orgWorkspaceCacheLoaded) return
        log.debug "Loading Seqera org/workspace cache"
        final entries = fetchUserWorkspacesAndOrgs(getUserId())
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
