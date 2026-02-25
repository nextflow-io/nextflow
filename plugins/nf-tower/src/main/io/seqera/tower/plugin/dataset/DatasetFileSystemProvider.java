/*
 * Copyright 2013-2024, Seqera Labs
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

package io.seqera.tower.plugin.dataset;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.AccessMode;
import java.nio.file.CopyOption;
import java.nio.file.DirectoryStream;
import java.nio.file.FileStore;
import java.nio.file.FileSystem;
import java.nio.file.FileSystemAlreadyExistsException;
import java.nio.file.LinkOption;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.ReadOnlyFileSystemException;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.spi.FileSystemProvider;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * NIO FileSystemProvider for the {@code dataset://} scheme.
 * <p>
 * Resolves dataset URIs of the form {@code dataset://name} or
 * {@code dataset://name?version=N} to their backing cloud storage
 * path via the Seqera Platform API, then delegates all I/O
 * operations to the resolved path's provider.
 * <p>
 * Phase 1: read-only — write operations throw {@link ReadOnlyFileSystemException}.
 *
 * @author Edmund Miller
 */
public class DatasetFileSystemProvider extends FileSystemProvider {

    private static final Logger log = LoggerFactory.getLogger(DatasetFileSystemProvider.class);

    private volatile DatasetFileSystem fileSystem;

    @Override
    public String getScheme() {
        return "dataset";
    }

    @Override
    public FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        if (fileSystem != null) {
            throw new FileSystemAlreadyExistsException("Dataset filesystem already exists");
        }
        synchronized (this) {
            if (fileSystem != null) {
                throw new FileSystemAlreadyExistsException("Dataset filesystem already exists");
            }
            fileSystem = new DatasetFileSystem(this, env);
            return fileSystem;
        }
    }

    @Override
    public FileSystem getFileSystem(URI uri) {
        if (fileSystem == null) {
            throw new java.nio.file.FileSystemNotFoundException("Dataset filesystem not yet created. Use newFileSystem() first");
        }
        return fileSystem;
    }

    @Override
    public Path getPath(URI uri) {
        if (!"dataset".equals(uri.getScheme())) {
            throw new IllegalArgumentException("URI scheme must be 'dataset': " + uri);
        }
        return new DatasetPath(getOrCreateFileSystem(), uri);
    }

    // -- read operations: delegate to resolved cloud path --

    @Override
    public InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final Path resolved = resolvedPath(path);
        log.debug("newInputStream for dataset path {} -> {}", path, resolved);
        return resolved.getFileSystem().provider().newInputStream(resolved, options);
    }

    @Override
    public SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        final Path resolved = resolvedPath(path);
        return resolved.getFileSystem().provider().newByteChannel(resolved, options, attrs);
    }

    @Override
    public DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        throw new UnsupportedOperationException("Dataset paths do not support directory listing");
    }

    @Override
    @SuppressWarnings("unchecked")
    public <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        final Path resolved = resolvedPath(path);
        return resolved.getFileSystem().provider().readAttributes(resolved, type, options);
    }

    @Override
    public Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        final Path resolved = resolvedPath(path);
        return resolved.getFileSystem().provider().readAttributes(resolved, attributes, options);
    }

    @Override
    public <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        try {
            final Path resolved = resolvedPath(path);
            return resolved.getFileSystem().provider().getFileAttributeView(resolved, type, options);
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to resolve dataset path: " + path, e);
        }
    }

    @Override
    public void checkAccess(Path path, AccessMode... modes) throws IOException {
        for (AccessMode mode : modes) {
            if (mode == AccessMode.WRITE) {
                throw new ReadOnlyFileSystemException();
            }
        }
        final Path resolved = resolvedPath(path);
        resolved.getFileSystem().provider().checkAccess(resolved, modes);
    }

    // -- write operations: read-only in phase 1 --

    @Override
    public void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        throw new ReadOnlyFileSystemException();
    }

    @Override
    public void delete(Path path) throws IOException {
        throw new ReadOnlyFileSystemException();
    }

    @Override
    public void copy(Path source, Path target, CopyOption... options) throws IOException {
        throw new ReadOnlyFileSystemException();
    }

    @Override
    public void move(Path source, Path target, CopyOption... options) throws IOException {
        throw new ReadOnlyFileSystemException();
    }

    @Override
    public boolean isSameFile(Path path1, Path path2) throws IOException {
        if (path1 instanceof DatasetPath && path2 instanceof DatasetPath) {
            return path1.equals(path2);
        }
        return false;
    }

    @Override
    public boolean isHidden(Path path) throws IOException {
        return false;
    }

    @Override
    public FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new ReadOnlyFileSystemException();
    }

    // -- internal helpers --

    private DatasetFileSystem getOrCreateFileSystem() {
        if (fileSystem == null) {
            synchronized (this) {
                if (fileSystem == null) {
                    fileSystem = new DatasetFileSystem(this, null);
                }
            }
        }
        return fileSystem;
    }

    /**
     * Resolve a dataset path to its backing cloud storage path.
     * Caches the resolved path on the DatasetPath instance.
     */
    private Path resolvedPath(Path path) throws IOException {
        if (!(path instanceof DatasetPath)) {
            throw new IllegalArgumentException("Path must be a DatasetPath: " + path);
        }
        return ((DatasetPath) path).getResolvedPath();
    }
}
