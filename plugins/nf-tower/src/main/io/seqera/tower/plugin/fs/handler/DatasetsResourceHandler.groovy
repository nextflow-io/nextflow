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

import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.time.Instant

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.DatasetVersionDto
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import io.seqera.tower.plugin.fs.ResourceTypeHandler
import io.seqera.tower.plugin.fs.SeqeraFileAttributes
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath

/**
 * {@link ResourceTypeHandler} for the {@code datasets} resource type.
 * Owns its own dataset/version caches and {@code @version} suffix parsing.
 */
@Slf4j
@CompileStatic
class DatasetsResourceHandler implements ResourceTypeHandler {

    public static final String TYPE = 'datasets'

    private final SeqeraFileSystem fs
    private final SeqeraDatasetClient client

    /** workspaceId → list of DatasetDto */
    private final Map<Long, List<DatasetDto>> datasetCache = new LinkedHashMap<>()
    /** datasetId → list of DatasetVersionDto */
    private final Map<String, List<DatasetVersionDto>> versionCache = new LinkedHashMap<>()

    DatasetsResourceHandler(SeqeraFileSystem fs, SeqeraDatasetClient client) {
        this.fs = fs
        this.client = client
    }

    @Override
    String getResourceType() { TYPE }

    @Override
    Iterable<Path> list(SeqeraPath dir) throws IOException {
        final d = dir.depth()
        if (d == 3) {
            final workspaceId = fs.resolveWorkspaceId(dir.org, dir.workspace)
            return resolveDatasets(workspaceId).collect { DatasetDto ds -> dir.resolve(ds.name) as Path }
        }
        throw new IllegalArgumentException("datasets handler does not list depth $d paths: $dir")
    }

    @Override
    SeqeraFileAttributes readAttributes(SeqeraPath p) throws IOException {
        final d = p.depth()
        if (d == 3) {
            fs.resolveWorkspaceId(p.org, p.workspace) // validates
            return new SeqeraFileAttributes(true)
        }
        if (d != 4)
            throw new NoSuchFileException(p.toString(), null, "Invalid dataset path depth: $d")
        final names = parseNameAndVersion(p.trail[0])
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final dataset = resolveDataset(workspaceId, names[0])
        if (!dataset)
            throw new NoSuchFileException(p.toString(), null, "Dataset '${names[0]}' not found in workspace ${p.workspace}")
        return new SeqeraFileAttributes(
                0L,
                dataset.lastUpdated?.toInstant() ?: Instant.EPOCH,
                dataset.dateCreated?.toInstant() ?: Instant.EPOCH,
                dataset.id)
    }

    @Override
    InputStream newInputStream(SeqeraPath p) throws IOException {
        if (p.depth() != 4)
            throw new IllegalArgumentException("Operation `newInputStream` requires a dataset path (depth 4): $p")
        final names = parseNameAndVersion(p.trail[0])
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final dataset = resolveDataset(workspaceId, names[0])
        if (!dataset)
            throw new NoSuchFileException(p.toString(), null, "Dataset '${names[0]}' not found in workspace ${p.workspace}")
        final version = resolveVersion(dataset, names[1], p)
        log.debug "Downloading dataset '${names[0]}' version ${version.version} (${version.fileName}) from workspace $workspaceId"
        return client.downloadDataset(dataset.id, String.valueOf(version.version), version.fileName, dataset.workspaceId)
    }

    @Override
    void checkAccess(SeqeraPath p, AccessMode... modes) throws IOException {
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE || m == AccessMode.EXECUTE)
                throw new AccessDeniedException(p.toString(), null, "seqera:// datasets are read-only")
        }
        readAttributes(p)
    }

    // ---- helpers ----

    /**
     * Split a trail segment into {@code [name, version]}. Version is {@code null} when
     * the segment does not contain an {@code @}.
     */
    private static String[] parseNameAndVersion(String segment) {
        final atIdx = segment.lastIndexOf('@')
        if (atIdx > 0)
            return [segment.substring(0, atIdx), segment.substring(atIdx + 1)] as String[]
        return [segment, null] as String[]
    }

    private synchronized List<DatasetDto> resolveDatasets(long workspaceId) {
        def cached = datasetCache.get(workspaceId)
        if (cached == null) {
            cached = client.listDatasets(workspaceId)
            datasetCache.put(workspaceId, cached)
        }
        return cached
    }

    private synchronized DatasetDto resolveDataset(long workspaceId, String name) {
        return resolveDatasets(workspaceId).find { DatasetDto d -> d.name == name }
    }

    private synchronized List<DatasetVersionDto> resolveVersions(String datasetId, long workspaceId) {
        def cached = versionCache.get(datasetId)
        if (cached == null) {
            cached = client.listVersions(datasetId, workspaceId)
            versionCache.put(datasetId, cached)
        }
        return cached
    }

    private DatasetVersionDto resolveVersion(DatasetDto dataset, String pinnedVersion, SeqeraPath p) throws IOException {
        final versions = resolveVersions(dataset.id, dataset.workspaceId)
        if (versions.isEmpty())
            throw new NoSuchFileException(p.toString(), null, "No versions available for dataset '${dataset.name}'")
        if (pinnedVersion) {
            final found = versions.find { DatasetVersionDto v -> String.valueOf(v.version) == pinnedVersion }
            if (!found)
                throw new NoSuchFileException(p.toString(), null, "Version '$pinnedVersion' not found for dataset '${dataset.name}'")
            return found
        }
        final latest = versions.findAll { DatasetVersionDto v -> !v.disabled }.max { DatasetVersionDto v -> v.version }
        if (!latest)
            throw new NoSuchFileException(p.toString(), null, "No enabled versions for dataset '${dataset.name}'")
        return latest
    }
}
