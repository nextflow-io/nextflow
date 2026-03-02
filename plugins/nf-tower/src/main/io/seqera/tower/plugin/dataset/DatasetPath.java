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
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.Iterator;
import java.util.Objects;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A {@link Path} representing a Seqera Platform dataset reference.
 * <p>
 * URI format: {@code dataset://name} or {@code dataset://name?version=N}
 * <p>
 * The path lazily resolves to the backing cloud storage path
 * (S3/GCS/Azure) via the Platform API on first I/O access.
 *
 * @author Edmund Miller
 */
public class DatasetPath implements Path {

    private static final Logger log = LoggerFactory.getLogger(DatasetPath.class);

    private final DatasetFileSystem fileSystem;
    private final String datasetName;
    private final String version; // null = latest
    private final URI uri;

    /** Cached resolved cloud path — populated lazily on first I/O */
    private volatile Path resolvedPath;

    /**
     * Construct from a URI (e.g. from provider.getPath(URI))
     */
    DatasetPath(DatasetFileSystem fileSystem, URI uri) {
        this.fileSystem = fileSystem;
        this.uri = uri;
        // dataset://my-samplesheet or dataset:///my-samplesheet
        // host = dataset name, or if host is null, first path segment is the name
        String name = uri.getHost();
        if (name == null || name.isEmpty()) {
            // handle dataset:///name form
            String path = uri.getPath();
            if (path != null && path.startsWith("/")) {
                path = path.substring(1);
            }
            name = path;
        }
        this.datasetName = name;
        // parse ?version=N from query string
        this.version = parseVersion(uri.getQuery());
    }

    /**
     * Construct from string path (e.g. from fileSystem.getPath())
     */
    DatasetPath(DatasetFileSystem fileSystem, String path) {
        this.fileSystem = fileSystem;
        // strip leading slash if present
        if (path.startsWith("/")) {
            path = path.substring(1);
        }
        // check for version suffix: name@version
        int atIdx = path.indexOf('@');
        if (atIdx > 0) {
            this.datasetName = path.substring(0, atIdx);
            this.version = path.substring(atIdx + 1);
        }
        else {
            this.datasetName = path;
            this.version = null;
        }
        this.uri = URI.create("dataset://" + datasetName + (version != null ? "?version=" + version : ""));
    }

    public String getDatasetName() {
        return datasetName;
    }

    public String getVersion() {
        return version;
    }

    /**
     * Resolve this dataset reference to the backing cloud storage path.
     * Lazily initialized and cached.
     */
    Path getResolvedPath() throws IOException {
        if (resolvedPath == null) {
            synchronized (this) {
                if (resolvedPath == null) {
                    log.debug("Resolving dataset '{}' version={}", datasetName, version != null ? version : "latest");
                    resolvedPath = DatasetResolver.resolve(datasetName, version);
                    log.debug("Resolved dataset '{}' -> {}", datasetName, resolvedPath);
                }
            }
        }
        return resolvedPath;
    }

    // -- Path interface --

    @Override
    public FileSystem getFileSystem() {
        return fileSystem;
    }

    @Override
    public boolean isAbsolute() {
        return true;
    }

    @Override
    public Path getRoot() {
        return null;
    }

    @Override
    public Path getFileName() {
        // The dataset name is the "file name"
        return new DatasetPath(fileSystem, datasetName);
    }

    @Override
    public Path getParent() {
        return null;
    }

    @Override
    public int getNameCount() {
        return 1;
    }

    @Override
    public Path getName(int index) {
        if (index != 0) {
            throw new IllegalArgumentException("Invalid name index: " + index);
        }
        return this;
    }

    @Override
    public Path subpath(int beginIndex, int endIndex) {
        if (beginIndex != 0 || endIndex != 1) {
            throw new IllegalArgumentException("Invalid subpath range");
        }
        return this;
    }

    @Override
    public boolean startsWith(Path other) {
        return equals(other);
    }

    @Override
    public boolean endsWith(Path other) {
        return equals(other);
    }

    @Override
    public Path normalize() {
        return this;
    }

    @Override
    public Path resolve(Path other) {
        // dataset paths are leaf nodes, cannot resolve children
        throw new UnsupportedOperationException("Cannot resolve against a dataset path");
    }

    @Override
    public Path relativize(Path other) {
        throw new UnsupportedOperationException("Cannot relativize dataset paths");
    }

    @Override
    public URI toUri() {
        return uri;
    }

    @Override
    public Path toAbsolutePath() {
        return this;
    }

    @Override
    public Path toRealPath(LinkOption... options) throws IOException {
        // Return the resolved cloud path as the "real" path
        return getResolvedPath();
    }

    @Override
    public WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public int compareTo(Path other) {
        if (other instanceof DatasetPath) {
            DatasetPath o = (DatasetPath) other;
            int cmp = datasetName.compareTo(o.datasetName);
            if (cmp != 0) return cmp;
            if (version == null && o.version == null) return 0;
            if (version == null) return -1;
            if (o.version == null) return 1;
            return version.compareTo(o.version);
        }
        return toString().compareTo(other.toString());
    }

    @Override
    public Iterator<Path> iterator() {
        return java.util.List.<Path>of(this).iterator();
    }

    @Override
    public String toString() {
        return "dataset://" + datasetName + (version != null ? "?version=" + version : "");
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DatasetPath)) return false;
        DatasetPath that = (DatasetPath) o;
        return Objects.equals(datasetName, that.datasetName) && Objects.equals(version, that.version);
    }

    @Override
    public int hashCode() {
        return Objects.hash(datasetName, version);
    }

    // -- helpers --

    private static String parseVersion(String query) {
        if (query == null || query.isEmpty()) return null;
        for (String param : query.split("&")) {
            String[] kv = param.split("=", 2);
            if (kv.length == 2 && "version".equals(kv[0])) {
                return kv[1];
            }
        }
        return null;
    }
}
