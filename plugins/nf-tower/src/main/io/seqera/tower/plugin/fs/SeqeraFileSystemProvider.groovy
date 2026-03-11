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

import nextflow.Global
import nextflow.Session

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.StandardCopyOption
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.TowerFactory
import io.seqera.tower.plugin.dataset.DatasetDto
import io.seqera.tower.plugin.dataset.DatasetVersionDto
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

    static final String SCHEME = 'seqera'

    static final long MAX_UPLOAD_BYTES = 10L * 1024 * 1024  // 10 MB

    /** One filesystem per (endpoint + accessToken) key */
    private final Map<String, SeqeraFileSystem> fileSystems = new LinkedHashMap<>()

    @Override
    String getScheme() { SCHEME }

    // ---- FileSystem lifecycle ----

    @Override
    synchronized FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        checkScheme(uri)
        final key = fsKey(uri, env)
        if (fileSystems.containsKey(key))
            return fileSystems.get(key)
        checkTowerEnabled()
        final TowerClient towerClient = TowerFactory.client()
        if (!towerClient)
            throw new IllegalStateException("seqera:// paths require Seqera Platform to be enabled (tower.accessToken and tower.endpoint)")
        final client = new SeqeraDatasetClient(towerClient)
        final seqeraFs = new SeqeraFileSystem(this, client)
        fileSystems.put(key, seqeraFs)
        return seqeraFs
    }

    @Override
    FileSystem getFileSystem(URI uri) {
        checkScheme(uri)
        final fs = fileSystems.values().find { true }
        if (!fs) throw new FileSystemNotFoundException("No seqera:// filesystem has been created yet")
        return fs
    }

    synchronized SeqeraFileSystem getOrCreateFileSystem(URI uri, Map<String, ?> env) {
        checkScheme(uri)
        final key = fsKey(uri, env)
        SeqeraFileSystem fs = fileSystems.get(key)
        if (!fs) {
            final envMap = env ?: Collections.<String, Object>emptyMap()
            fs = (SeqeraFileSystem) newFileSystem(uri, envMap as Map<String, ?>)
        }
        return fs
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
            throw new IllegalArgumentException("newInputStream requires a dataset path (depth 4): $path")
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final workspaceId = fs.resolveWorkspaceId(sp.org, sp.workspace)
        final dataset = fs.resolveDataset(workspaceId, sp.datasetName)
        final version = resolveVersion(fs, dataset, sp.version)
        log.debug "Downloading dataset '${sp.datasetName}' version ${version.version} (${version.fileName}) from workspace $workspaceId"
        return fs.client.downloadDataset(dataset.id, String.valueOf(version.version), version.fileName, dataset.workspaceId)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        if (options?.contains(StandardOpenOption.WRITE) || options?.contains(StandardOpenOption.APPEND))
            throw new UnsupportedOperationException("WRITE channel not yet supported for seqera:// paths — use newOutputStream")
        final inputStream = newInputStream(path)
        return new InputStreamSeekableByteChannel(inputStream)
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
        return (A) new SeqeraFileAttributes(dataset)
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("readAttributes(String) not supported by seqera:// filesystem")
    }

    // ---- Access check ----

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        final sp = toSeqeraPath(path)
        for (AccessMode m : modes) {
            if (m == AccessMode.EXECUTE)
                throw new AccessDeniedException(path.toString(), null, "EXECUTE not supported on seqera:// paths")
        }
        // For READ and WRITE, verify the path resolves without throwing NoSuchFileException
        if (sp.depth() >= 2) {
            final fs = sp.getFileSystem() as SeqeraFileSystem
            fs.resolveWorkspaceId(sp.org, sp.workspace)
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
            try { filter.accept(p) } catch (IOException e) { false }
        } : entries

        return new DirectoryStream<Path>() {
            @Override Iterator<Path> iterator() { filtered.iterator() }
            @Override void close() {}
        }
    }

    // ---- Write operations ----

    @Override
    OutputStream newOutputStream(Path path, OpenOption... options) throws IOException {
        final sp = toSeqeraPath(path)
        if (sp.depth() != 4)
            throw new IllegalArgumentException("newOutputStream requires a dataset path (depth 4): $path")
        return new DatasetOutputStream(sp)
    }

    // ---- Copy ----

    @Override
    void copy(Path source, Path target, CopyOption... options) throws IOException {
        final sp = toSeqeraPath(source)
        if (target instanceof SeqeraPath) {
            // within-provider: download then upload
            final bytes = newInputStream(source).bytes
            final out = newOutputStream(target)
            out.write(bytes)
            out.close()
        } else {
            // cross-provider (seqera → local): stream to target
            final inputStream = newInputStream(source)
            Files.copy(inputStream, target, StandardCopyOption.REPLACE_EXISTING)
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

    private static String fsKey(URI uri, Map<String, ?> env) {
        // Key by scheme — one filesystem per JVM (TowerClient is also a singleton per session)
        return uri.scheme
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

    private static DatasetVersionDto resolveVersion(SeqeraFileSystem fs, DatasetDto dataset, String pinnedVersion) throws IOException {
        final versions = fs.client.listVersions(dataset.id, dataset.workspaceId)
        if (versions.isEmpty())
            throw new NoSuchFileException("seqera:///datasets/${dataset.name}", null, "No versions available for dataset '${dataset.name}'")
        if (pinnedVersion) {
            final found = versions.find { DatasetVersionDto v -> String.valueOf(v.version) == pinnedVersion }
            if (!found)
                throw new NoSuchFileException("seqera:///datasets/${dataset.name}@${pinnedVersion}", null, "Version '${pinnedVersion}' not found for dataset '${dataset.name}'")
            return found
        }
        // Latest non-disabled version
        final latest = versions.findAll { DatasetVersionDto v -> !v.disabled }
                               .max { DatasetVersionDto v -> v.version }
        if (!latest)
            throw new NoSuchFileException("seqera:///datasets/${dataset.name}", null, "No enabled versions for dataset '${dataset.name}'")
        return latest
    }

    /**
     * Check if tower is enabled and enables it if not
     */
    private void checkTowerEnabled() {
        Session session = Global.session as Session
        if( !session )
             throw new IllegalStateException("Session not found creating a Seqera file system provider")
        if( !session.config.tower ) {
            session.config.tower = [enabled: true]
            return
        }
        final tower = session.config.tower as Map
        if ( !tower.enabled )
            tower.enabled = true

    }

    // ---- inner classes ----

    /**
     * Minimal {@link SeekableByteChannel} backed by an {@link InputStream}.
     * Supports sequential reads only (no seek/position).
     */
    @CompileStatic
    private static class InputStreamSeekableByteChannel implements SeekableByteChannel {
        private final InputStream inputStream
        private long position0 = 0L
        private boolean open = true

        InputStreamSeekableByteChannel(InputStream inputStream) {
            this.inputStream = inputStream
        }

        @Override
        int read(ByteBuffer dst) throws IOException {
            final bytes = new byte[dst.remaining()]
            final n = inputStream.read(bytes)
            if (n > 0) {
                dst.put(bytes, 0, n)
                position0 += n
            }
            return n
        }

        @Override
        int write(ByteBuffer src) { throw new UnsupportedOperationException() }

        @Override
        long position() { position0 }

        @Override
        SeekableByteChannel position(long newPosition) { throw new UnsupportedOperationException("seek not supported") }

        @Override
        long size() { -1L }

        @Override
        SeekableByteChannel truncate(long size) { throw new UnsupportedOperationException() }

        @Override
        boolean isOpen() { open }

        @Override
        void close() throws IOException {
            open = false
            inputStream.close()
        }
    }

    /**
     * Buffered {@link OutputStream} that uploads to the Seqera Platform on {@code close()}.
     * Enforces the 10 MB platform size limit before making the API call.
     */
    @CompileStatic
    private static class DatasetOutputStream extends OutputStream {
        private final SeqeraPath path
        private final ByteArrayOutputStream buffer = new ByteArrayOutputStream()
        private boolean closed = false

        DatasetOutputStream(SeqeraPath path) {
            this.path = path
        }

        @Override
        void write(int b) throws IOException {
            buffer.write(b)
        }

        @Override
        void write(byte[] b, int off, int len) throws IOException {
            buffer.write(b, off, len)
        }

        @Override
        void close() throws IOException {
            if (closed) return
            closed = true
            final bytes = buffer.toByteArray()
            if (bytes.length > MAX_UPLOAD_BYTES)
                throw new IOException("Dataset '${path.datasetName}' exceeds the 10 MB Seqera Platform limit (${bytes.length} bytes)")

            final fs = path.getFileSystem() as SeqeraFileSystem
            final workspaceId = fs.resolveWorkspaceId(path.org, path.workspace)

            // Resolve or create the dataset
            String datasetId
            try {
                final existing = fs.resolveDataset(workspaceId, path.datasetName)
                datasetId = existing.id
            } catch (NoSuchFileException ignored) {
                final created = fs.client.createDataset(workspaceId, path.datasetName)
                datasetId = created.id
            }

            final fileName = path.datasetName
            fs.client.uploadDataset(datasetId, bytes, fileName, false)
            fs.invalidateDatasetCache(workspaceId)
        }
    }

    /** Thrown by newDirectoryStream for non-directory paths */
    @CompileStatic
    private static class NotDirectoryException extends IOException {
        NotDirectoryException(String path) { super("Not a directory: $path") }
    }
}
