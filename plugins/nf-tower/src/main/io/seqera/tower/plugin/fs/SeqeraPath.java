/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.tower.plugin.fs;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.ProviderMismatchException;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.Iterator;
import java.util.Objects;

/**
 * Seqera Platform Data-Link file system path
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraPath implements Path {

    private static final Logger log = LoggerFactory.getLogger(SeqeraPath.class);

    public static final String SEPARATOR = "/";
    public static final String SCHEME = "seqera";
    public static final String SEQERA_PROT = SCHEME + "://";

    private static final String[] EMPTY = new String[]{};

    private final SeqeraFileSystem fileSystem;
    private final String path;

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected SeqeraPath() {
        this.fileSystem = null;
        this.path = null;
    }

    public SeqeraPath(SeqeraFileSystem fs, URI uri) {
        if (!SCHEME.equals(uri.getScheme())) {
            throw new IllegalArgumentException("Invalid Seqera URI - scheme must be '" + SCHEME + "'");
        }
        this.fileSystem = fs;
        // URI format: seqera://datalink-name/path/to/file
        this.path = normalizePath(uri.getPath() != null ? uri.getPath() : "");
    }

    public SeqeraPath(SeqeraFileSystem fs, String path) {
        this.fileSystem = fs;
        this.path = normalizePath(path);
    }

    private static String normalizePath(String path) {
        if (path == null || path.isEmpty() || path.equals(SEPARATOR)) {
            return "";
        }
        // Normalize and remove leading/trailing separators
        path = path.trim();
        if (path.startsWith(SEPARATOR)) {
            path = path.substring(1);
        }
        if (path.endsWith(SEPARATOR) && path.length() > 1) {
            path = path.substring(0, path.length() - 1);
        }
        return path;
    }

    public static boolean isSeqeraUri(String path) {
        return path != null && path.startsWith(SEQERA_PROT);
    }

    public static URI asUri(String path) {
        if (path == null) {
            throw new IllegalArgumentException("Missing 'path' argument");
        }
        if (!path.startsWith(SEQERA_PROT)) {
            throw new IllegalArgumentException("Invalid Seqera file system path URI - it must start with '" + SEQERA_PROT + "' prefix");
        }
        return URI.create(path);
    }

    @Override
    public FileSystem getFileSystem() {
        return fileSystem;
    }

    @Override
    public boolean isAbsolute() {
        return fileSystem != null;
    }

    @Override
    public Path getRoot() {
        return new SeqeraPath(fileSystem, "");
    }

    @Override
    public Path getFileName() {
        if (path == null || path.isEmpty()) {
            return null;
        }
        final int idx = path.lastIndexOf(SEPARATOR);
        return idx == -1
            ? new SeqeraPath(null, path)
            : new SeqeraPath(null, path.substring(idx + 1));
    }

    @Override
    public Path getParent() {
        if (path == null || path.isEmpty()) {
            return null;
        }
        final int idx = path.lastIndexOf(SEPARATOR);
        if (idx == -1) {
            return new SeqeraPath(fileSystem, "");
        }
        return new SeqeraPath(fileSystem, path.substring(0, idx));
    }

    @Override
    public int getNameCount() {
        if (path == null || path.isEmpty()) {
            return 0;
        }
        return path.split(SEPARATOR).length;
    }

    @Override
    public Path getName(int index) {
        if (index < 0) {
            throw new IllegalArgumentException("Path name index cannot be less than zero - offending value: " + index);
        }
        final String[] parts = path.split(SEPARATOR);
        if (index >= parts.length) {
            throw new IllegalArgumentException("Index out of bounds: " + index);
        }
        return new SeqeraPath(null, parts[index]);
    }

    @Override
    public Path subpath(int beginIndex, int endIndex) {
        if (beginIndex < 0) {
            throw new IllegalArgumentException("subpath begin index cannot be less than zero - offending value: " + beginIndex);
        }
        if (beginIndex >= endIndex) {
            throw new IllegalArgumentException("begin index must be less than end index");
        }
        final String[] parts = path.split(SEPARATOR);
        if (endIndex > parts.length) {
            throw new IllegalArgumentException("end index out of bounds: " + endIndex);
        }
        final StringBuilder sb = new StringBuilder();
        for (int i = beginIndex; i < endIndex; i++) {
            if (i > beginIndex) {
                sb.append(SEPARATOR);
            }
            sb.append(parts[i]);
        }
        return new SeqeraPath(beginIndex == 0 ? fileSystem : null, sb.toString());
    }

    @Override
    public Path normalize() {
        return new SeqeraPath(fileSystem, normalizePath(path));
    }

    @Override
    public Path resolve(Path other) {
        if (!(other instanceof SeqeraPath)) {
            throw new ProviderMismatchException();
        }

        final SeqeraPath that = (SeqeraPath) other;

        if (that.fileSystem != null && this.fileSystem != that.fileSystem) {
            return other;
        }
        if (that.isAbsolute()) {
            return that;
        }
        if (this.path == null || this.path.isEmpty()) {
            return that;
        }
        return new SeqeraPath(fileSystem, path + SEPARATOR + that.path);
    }

    @Override
    public Path resolve(String other) {
        if (other == null || other.isEmpty()) {
            return this;
        }
        // If it's a Seqera URI, parse it
        if (other.startsWith(SEQERA_PROT)) {
            final Path that = fileSystem.provider().getPath(asUri(other));
            return resolve(that);
        }
        // Otherwise, treat as relative path
        return resolve(new SeqeraPath(null, other));
    }

    @Override
    public Path relativize(Path other) {
        if (!(other instanceof SeqeraPath)) {
            throw new ProviderMismatchException();
        }
        final SeqeraPath that = (SeqeraPath) other;
        if (this.isAbsolute() != that.isAbsolute()) {
            throw new IllegalArgumentException("Cannot compare absolute with relative paths");
        }

        final String[] thisParts = this.path != null && !this.path.isEmpty() ? this.path.split(SEPARATOR) : new String[0];
        final String[] thatParts = that.path != null && !that.path.isEmpty() ? that.path.split(SEPARATOR) : new String[0];

        // Find common prefix
        int common = 0;
        while (common < thisParts.length && common < thatParts.length && thisParts[common].equals(thatParts[common])) {
            common++;
        }

        // Build relative path
        final StringBuilder result = new StringBuilder();
        final int upLevels = thisParts.length - common;
        for (int i = 0; i < upLevels; i++) {
            if (i > 0) {
                result.append(SEPARATOR);
            }
            result.append("..");
        }
        for (int i = common; i < thatParts.length; i++) {
            if (result.length() > 0) {
                result.append(SEPARATOR);
            }
            result.append(thatParts[i]);
        }

        return new SeqeraPath(null, result.toString());
    }

    @Override
    public URI toUri() {
        final String dataLinkName = fileSystem != null ? fileSystem.getDataLinkName() : "unknown";
        final String uriPath = path != null && !path.isEmpty() ? SEPARATOR + path : "";
        return URI.create(SEQERA_PROT + dataLinkName + uriPath);
    }

    public String toUriString() {
        return toUri().toString();
    }

    @Override
    public Path toAbsolutePath() {
        return this;
    }

    @Override
    public Path toRealPath(LinkOption... options) throws IOException {
        // For Seqera paths, real path is the same as the path itself
        return this;
    }

    @Override
    public File toFile() {
        throw new UnsupportedOperationException("toFile not supported by SeqeraPath");
    }

    @Override
    public WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException("Register not supported by SeqeraPath");
    }

    @Override
    public WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        throw new UnsupportedOperationException("Register not supported by SeqeraPath");
    }

    @Override
    public Iterator<Path> iterator() {
        return new Iterator<Path>() {
            private int index = 0;
            private final String[] parts = path != null && !path.isEmpty() ? path.split(SEPARATOR) : new String[0];

            @Override
            public boolean hasNext() {
                return index < parts.length;
            }

            @Override
            public Path next() {
                return new SeqeraPath(null, parts[index++]);
            }
        };
    }

    @Override
    public int compareTo(Path other) {
        return toString().compareTo(other.toString());
    }

    @Override
    public boolean startsWith(Path other) {
        return startsWith(other.toString());
    }

    @Override
    public boolean startsWith(String other) {
        return path != null && path.startsWith(normalizePath(other));
    }

    @Override
    public boolean endsWith(Path other) {
        return endsWith(other.toString());
    }

    @Override
    public boolean endsWith(String other) {
        return path != null && path.endsWith(normalizePath(other));
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof SeqeraPath)) {
            return false;
        }
        final SeqeraPath that = (SeqeraPath) other;
        return Objects.equals(this.fileSystem, that.fileSystem) && Objects.equals(this.path, that.path);
    }

    @Override
    public int hashCode() {
        return Objects.hash(fileSystem, path);
    }

    @Override
    public String toString() {
        return path != null ? path : "";
    }

    /**
     * Get the path string for API calls
     */
    public String getPathForApi() {
        return path != null ? path : "";
    }
}
