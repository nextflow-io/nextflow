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
import java.nio.file.DirectoryIteratorException
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
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
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.TowerFactory
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import io.seqera.tower.plugin.fs.handler.DataLinksResourceHandler
import io.seqera.tower.plugin.fs.handler.DatasetsResourceHandler

/**
 * NIO {@link FileSystemProvider} for the {@code seqera://} scheme. Registered via
 * {@code META-INF/services/java.nio.file.spi.FileSystemProvider}.
 *
 * Generic for depth &ge; 3: dispatches to a {@link ResourceTypeHandler} selected by
 * {@code SeqeraPath.resourceType}. The handlers own all resource-specific logic.
 */
@Slf4j
@CompileStatic
class SeqeraFileSystemProvider extends FileSystemProvider {

    public static final String SCHEME = 'seqera'

    private volatile SeqeraFileSystem fileSystem

    @Override
    String getScheme() { SCHEME }

    // ---- lifecycle ----

    @Override
    synchronized FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        checkScheme(uri)
        if (fileSystem)
            throw new FileSystemAlreadyExistsException("File system `seqera://` already exists")
        final TowerClient tc = TowerFactory.client()
        if (!tc)
            throw new IllegalStateException("File system `seqera://` requires the Seqera Platform access token — use `tower.accessToken` config option or TOWER_ACCESS_TOKEN env variable")
        final datasetClient = new SeqeraDatasetClient(tc)
        fileSystem = new SeqeraFileSystem(this)
        fileSystem.setOrgWorkspaceClient(datasetClient)
        fileSystem.registerHandler(new DatasetsResourceHandler(fileSystem, datasetClient))
        fileSystem.registerHandler(new DataLinksResourceHandler(fileSystem, new SeqeraDataLinkClient(tc)))
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
        if (!fileSystem) newFileSystem(uri, env ?: Collections.<String, Object>emptyMap())
        return fileSystem
    }

    @Override
    SeqeraPath getPath(URI uri) {
        final fs = getOrCreateFileSystem(uri, Collections.emptyMap())
        return new SeqeraPath(fs, uri.toString())
    }

    // ---- read ----

    @Override
    InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final sp = toSeqeraPath(path)
        if (sp.depth() < 3)
            throw new IllegalArgumentException("newInputStream requires a leaf path: $path")
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final h = fs.getHandler(sp.resourceType)
        if (!h)
            throw new NoSuchFileException(path.toString(), null, "Unsupported resource type: ${sp.resourceType}")
        return h.newInputStream(sp)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        if (options?.contains(StandardOpenOption.WRITE) || options?.contains(StandardOpenOption.APPEND))
            throw new UnsupportedOperationException("seqera:// filesystem is read-only")
        return new DatasetInputStream(newInputStream(path))
    }

    // ---- attributes ----

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        if (!BasicFileAttributes.isAssignableFrom(type))
            throw new UnsupportedOperationException("Attribute type not supported: $type")
        final sp = toSeqeraPath(path)
        if (sp.cachedAttributes)
            return (A) sp.cachedAttributes
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final d = sp.depth()
        if (d < 3) {
            validateSharedDirectoryExists(fs, sp)
            return (A) new SeqeraFileAttributes(true)
        }
        final h = fs.getHandler(sp.resourceType)
        if (!h)
            throw new NoSuchFileException(path.toString(), null, "Unsupported resource type: ${sp.resourceType}")
        return (A) h.readAttributes(sp)
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Operation `readAttributes(String)` not supported by `seqera://` file system")
    }

    // ---- access ----

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        final sp = toSeqeraPath(path)
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE || m == AccessMode.EXECUTE)
                throw new AccessDeniedException(path.toString(), null, "seqera:// filesystem is read-only")
        }
        final d = sp.depth()
        if (d == 0) return
        if (d < 3) {
            validateSharedDirectoryExists(sp.getFileSystem() as SeqeraFileSystem, sp)
            return
        }
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final h = fs.getHandler(sp.resourceType)
        if (!h)
            throw new NoSuchFileException(path.toString(), null, "Unsupported resource type: ${sp.resourceType}")
        h.checkAccess(sp, modes)
    }

    // ---- directory stream ----

    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final sp = toSeqeraPath(dir)
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final d = sp.depth()
        Iterable<Path> entries
        if (d == 0) {
            fs.loadOrgWorkspaceCache()
            entries = fs.listOrgNames().collect { String org -> sp.resolve(org) as Path }
        } else if (d == 1) {
            fs.loadOrgWorkspaceCache()
            entries = fs.listWorkspaceNames(sp.org).collect { String ws -> sp.resolve(ws) as Path }
        } else if (d == 2) {
            fs.resolveWorkspaceId(sp.org, sp.workspace)
            entries = fs.getResourceTypes().collect { String rt -> sp.resolve(rt) as Path }
        } else {
            final h = fs.getHandler(sp.resourceType)
            if (!h)
                throw new NoSuchFileException(dir.toString(), null, "Unsupported resource type: ${sp.resourceType}")
            entries = h.list(sp)
        }

        final source = entries
        return new DirectoryStream<Path>() {
            @Override
            Iterator<Path> iterator() {
                final inner = source.iterator()
                if (!filter) return inner
                return new FilteredIterator<Path>(inner, filter)
            }
            @Override void close() {}
        }
    }

    /** Lazy filtering iterator: calls the filter as each element is consumed. */
    @CompileStatic
    private static class FilteredIterator<T> implements Iterator<T> {
        private final Iterator<T> inner
        private final DirectoryStream.Filter<? super T> filter
        private T buffered
        private boolean hasBuffered = false

        FilteredIterator(Iterator<T> inner, DirectoryStream.Filter<? super T> filter) {
            this.inner = inner
            this.filter = filter
        }

        @Override
        boolean hasNext() {
            while (!hasBuffered && inner.hasNext()) {
                final candidate = inner.next()
                try {
                    if (filter.accept(candidate)) {
                        buffered = candidate
                        hasBuffered = true
                    }
                } catch (IOException e) {
                    throw new DirectoryIteratorException(e)
                }
            }
            return hasBuffered
        }

        @Override
        T next() {
            if (!hasNext()) throw new NoSuchElementException()
            final out = buffered
            buffered = null
            hasBuffered = false
            return out
        }
    }

    // ---- copy ----

    @Override
    void copy(Path source, Path target, CopyOption... options) throws IOException {
        toSeqeraPath(source)
        if (target instanceof SeqeraPath)
            throw new UnsupportedOperationException("seqera:// filesystem is read-only")
        try (final InputStream is = newInputStream(source)) {
            Files.copy(is, target, options)
        }
    }

    // ---- unsupported mutations ----

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

    // ---- misc ----

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

    // ---- helpers ----

    private static SeqeraPath toSeqeraPath(Path path) {
        if (path !instanceof SeqeraPath) throw new ProviderMismatchException()
        return (SeqeraPath) path
    }

    private static void checkScheme(URI uri) {
        if (uri.scheme?.toLowerCase() != SCHEME)
            throw new IllegalArgumentException("Not a seqera:// URI: $uri")
    }

    private static void validateSharedDirectoryExists(SeqeraFileSystem fs, SeqeraPath sp) throws NoSuchFileException {
        final d = sp.depth()
        if (d == 0) return
        fs.loadOrgWorkspaceCache()
        if (d >= 1 && !fs.listOrgNames().contains(sp.org))
            throw new NoSuchFileException("seqera://${sp.org}", null, "Organisation not found")
        if (d >= 2)
            fs.resolveWorkspaceId(sp.org, sp.workspace)
    }
}
