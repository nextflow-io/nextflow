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

import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.DirectoryIteratorException
import java.nio.file.NoSuchFileException
import java.nio.file.NotDirectoryException
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.DatasetVersionDto
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.TowerFactory
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient

/**
 * NIO {@link FileSystemProvider} for the {@code seqera://} scheme.
 * Registered via {@code META-INF/services/java.nio.file.spi.FileSystemProvider}.
 *
 * Enables Nextflow pipelines to read Seqera Platform datasets as ordinary file paths:
 * {@code seqera://<org>/<workspace>/datasets/<dataset-name>}
 *
 * Follows the {@code LinFileSystemProvider} pattern for structure.
 * Write support follows the {@code AzFileSystemProvider} buffered-upload pattern.
 *
 * @author Seqera Labs
 */
@Slf4j
@CompileStatic
class SeqeraFileSystemProvider extends FileSystemProvider {

    public static final String SCHEME = 'seqera'

    /** Single filesystem instance — TowerClient is a singleton per session */
    private volatile SeqeraFileSystem fileSystem

    @Override
    String getScheme() { SCHEME }

    // ---- FileSystem lifecycle ----

    @Override
    synchronized FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        checkScheme(uri)
        if (fileSystem)
            throw new FileSystemAlreadyExistsException("File system `seqera://` already exists")
        final TowerClient towerClient = TowerFactory.client()
        if (!towerClient)
            throw new IllegalStateException("File system `seqera://` requires the Seqera Platform access token to be provided - use `tower.accessToken` config option or TOWER_ACCESS_TOKEN env variable")
        final client = new SeqeraDatasetClient(towerClient)
        fileSystem = new SeqeraFileSystem(this, client)
        return fileSystem
    }

    @Override
    synchronized FileSystem getFileSystem(URI uri) {
        checkScheme(uri)
        if (!fileSystem) throw new FileSystemNotFoundException("No seqera:// filesystem has been created yet")
        return fileSystem
    }

    synchronized SeqeraFileSystem getOrCreateFileSystem(URI uri, Map<String, ?> env) {
        checkScheme(uri)
        if (!fileSystem) {
            final envMap = env ?: Collections.<String, Object>emptyMap()
            newFileSystem(uri, envMap as Map<String, ?>)
        }
        return fileSystem
    }

    @Override
    SeqeraPath getPath(URI uri) {
        final fs = getOrCreateFileSystem(uri, Collections.emptyMap())
        return new SeqeraPath(fs, uri.toString())
    }

    // ---- Read operations ----

    @Override
    InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final sp = toSeqeraPath(path)
        if (sp.depth() != 4)
            throw new IllegalArgumentException("Operation `newInputStream` requires a dataset path (depth 4): $path")
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final workspaceId = fs.resolveWorkspaceId(sp.org, sp.workspace)
        final dataset = fs.resolveDataset(workspaceId, sp.datasetName)
        if (!dataset)
            throw new NoSuchFileException(sp.toString(), null, "Dataset '${sp.datasetName}' not found in workspace $sp.workspace")
        final version = resolveVersion(fs, dataset, sp)
        log.debug "Downloading dataset '${sp.datasetName}' version ${version.version} (${version.fileName}) from workspace $workspaceId"
        return fs.client.downloadDataset(dataset.id, String.valueOf(version.version), version.fileName, dataset.workspaceId)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        if (options?.contains(StandardOpenOption.WRITE) || options?.contains(StandardOpenOption.APPEND))
            throw new UnsupportedOperationException("File system `seqera://` is read-only")
        final inputStream = newInputStream(path)
        return new DatasetInputStream(inputStream)
    }

    // ---- Metadata ----

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        if (!BasicFileAttributes.isAssignableFrom(type))
            throw new UnsupportedOperationException("Attribute type not supported: $type")
        final sp = toSeqeraPath(path)
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final d = sp.depth()
        if (d < 4) {
            // Virtual directory — validate the path exists (throws NoSuchFileException if not)
            validateDirectoryExists(fs, sp)
            return (A) new SeqeraFileAttributes(true)
        }
        // Dataset file
        final workspaceId = fs.resolveWorkspaceId(sp.org, sp.workspace)
        final dataset = fs.resolveDataset(workspaceId, sp.datasetName)
        if (!dataset)
            throw new NoSuchFileException(sp.toString(), null, "Dataset '${sp.datasetName}' not found in workspace $sp.workspace")
        return (A) new SeqeraFileAttributes(dataset)
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Operation `readAttributes(String)` not supported by `seqera://` file system")
    }

    // ---- Access check ----

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        final sp = toSeqeraPath(path)
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE || m == AccessMode.EXECUTE)
                throw new AccessDeniedException(path.toString(), null, "seqera:// filesystem is read-only")
        }
        // For READ, verify the path resolves without throwing NoSuchFileException
        if (sp.depth() >= 1) {
            final fs = sp.getFileSystem() as SeqeraFileSystem
            if (sp.depth() == 1) {
                fs.loadOrgWorkspaceCache()
                if (!fs.listOrgNames().contains(sp.org))
                    throw new NoSuchFileException(path.toString(), null, "Organisation not found")
            } else {
                fs.resolveWorkspaceId(sp.org, sp.workspace)
            }
        }
    }

    // ---- Directory stream ----

    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final sp = toSeqeraPath(dir)
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final d = sp.depth()
        List<Path> entries
        if (d == 0) {
            // Root: list distinct org names
            fs.loadOrgWorkspaceCache()
            entries = fs.listOrgNames().collect { String org -> sp.resolve(org) as Path }
        } else if (d == 1) {
            // Org: list workspace names
            fs.loadOrgWorkspaceCache()
            entries = fs.listWorkspaceNames(sp.org).collect { String ws -> sp.resolve(ws) as Path }
        } else if (d == 2) {
            // Workspace: static resource types
            entries = ['datasets'].collect { String rt -> sp.resolve(rt) as Path }
        } else if (d == 3) {
            // Resource type directory: list dataset names
            final workspaceId = fs.resolveWorkspaceId(sp.org, sp.workspace)
            entries = fs.resolveDatasets(workspaceId).collect { DatasetDto ds ->
                sp.resolve(ds.name) as Path
            }
        } else {
            throw new NotDirectoryException(dir.toString())
        }

        final filtered = filter ? entries.findAll { Path p ->
            try { filter.accept(p) }
            catch (IOException e) { throw new DirectoryIteratorException(e) }
        } : entries

        return new DirectoryStream<Path>() {
            @Override Iterator<Path> iterator() { filtered.iterator() }
            @Override void close() {}
        }
    }

    // ---- Copy ----

    @Override
    void copy(Path source, Path target, CopyOption... options) throws IOException {
        toSeqeraPath(source)
        if (target instanceof SeqeraPath)
            throw new UnsupportedOperationException("seqera:// filesystem is read-only")
        // cross-provider (seqera → local): stream to target
        try (final InputStream is = newInputStream(source)) {
            Files.copy(is, target, options)
        }
    }

    // ---- Unsupported mutations ----

    @Override
    void move(Path source, Path target, CopyOption... options) {
        throw new UnsupportedOperationException("move() not supported by seqera:// filesystem")
    }

    @Override
    void delete(Path path) {
        throw new UnsupportedOperationException("delete() not supported by seqera:// filesystem")
    }

    @Override
    void createDirectory(Path dir, FileAttribute<?>... attrs) {
        throw new UnsupportedOperationException("createDirectory() not supported by seqera:// filesystem")
    }

    // ---- Misc ----

    @Override
    boolean isSameFile(Path path, Path path2) throws IOException {
        return path == path2
    }

    @Override
    boolean isHidden(Path path) { false }

    @Override
    FileStore getFileStore(Path path) {
        throw new UnsupportedOperationException("getFileStore() not supported by seqera:// filesystem")
    }

    @Override
    <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        return null
    }

    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) {
        throw new UnsupportedOperationException("setAttribute() not supported by seqera:// filesystem")
    }

    // ---- private helpers ----

    private static SeqeraPath toSeqeraPath(Path path) {
        if (path !instanceof SeqeraPath)
            throw new ProviderMismatchException()
        return (SeqeraPath) path
    }

    private static void checkScheme(URI uri) {
        if (uri.scheme?.toLowerCase() != SCHEME)
            throw new IllegalArgumentException("Not a seqera:// URI: $uri")
    }

    private static void validateDirectoryExists(SeqeraFileSystem fs, SeqeraPath sp) throws NoSuchFileException {
        final d = sp.depth()
        if (d == 0) return
        // Depth 1+: ensure org/workspace cache is loaded
        fs.loadOrgWorkspaceCache()
        if (d >= 1 && !fs.listOrgNames().contains(sp.org))
            throw new NoSuchFileException("seqera://${sp.org}", null, "Organisation not found")
        if (d >= 2)
            fs.resolveWorkspaceId(sp.org, sp.workspace)
        if (d >= 3 && sp.resourceType != 'datasets')
            throw new NoSuchFileException("seqera://${sp.org}/${sp.workspace}/${sp.resourceType}", null, "Unsupported resource type")
    }

    private static DatasetVersionDto resolveVersion(SeqeraFileSystem fs, DatasetDto dataset, SeqeraPath sp) throws IOException {
        final pinnedVersion = sp.version
        final versions = fs.resolveVersions(dataset.id, dataset.workspaceId)
        if (versions.isEmpty())
            throw new NoSuchFileException(sp.toString(), null, "No versions available for dataset '${dataset.name}'")
        if (pinnedVersion) {
            final found = versions.find { DatasetVersionDto v -> String.valueOf(v.version) == pinnedVersion }
            if (!found)
                throw new NoSuchFileException(sp.toString(), null, "Version '${pinnedVersion}' not found for dataset '${dataset.name}'")
            return found
        }
        // Latest non-disabled version
        final latest = versions.findAll { DatasetVersionDto v -> !v.disabled }
                               .max { DatasetVersionDto v -> v.version }
        if (!latest)
            throw new NoSuchFileException(sp.toString(), null, "No enabled versions for dataset '${dataset.name}'")
        return latest
    }

}
