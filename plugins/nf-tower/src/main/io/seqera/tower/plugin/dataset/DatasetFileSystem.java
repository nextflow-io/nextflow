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
import java.nio.file.FileStore;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.WatchService;
import java.nio.file.attribute.UserPrincipalLookupService;
import java.nio.file.spi.FileSystemProvider;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

/**
 * Minimal read-only FileSystem for the {@code dataset://} scheme.
 * <p>
 * Datasets are single files resolved via the Seqera Platform API,
 * so most filesystem operations (roots, stores, path matching)
 * are either trivial or unsupported.
 *
 * @author Edmund Miller
 */
public class DatasetFileSystem extends FileSystem {

    static final String PATH_SEPARATOR = "/";

    private final DatasetFileSystemProvider provider;
    private final Map<String, ?> env;
    private volatile boolean open = true;

    DatasetFileSystem(DatasetFileSystemProvider provider, Map<String, ?> env) {
        this.provider = provider;
        this.env = env != null ? env : Collections.emptyMap();
    }

    @Override
    public FileSystemProvider provider() {
        return provider;
    }

    @Override
    public void close() throws IOException {
        open = false;
    }

    @Override
    public boolean isOpen() {
        return open;
    }

    @Override
    public boolean isReadOnly() {
        return true;
    }

    @Override
    public String getSeparator() {
        return PATH_SEPARATOR;
    }

    @Override
    public Iterable<Path> getRootDirectories() {
        return Collections.emptyList();
    }

    @Override
    public Iterable<FileStore> getFileStores() {
        return Collections.emptyList();
    }

    @Override
    public Set<String> supportedFileAttributeViews() {
        return Set.of("basic");
    }

    @Override
    public Path getPath(String first, String... more) {
        // Build a dataset path from string components
        StringBuilder sb = new StringBuilder(first);
        for (String part : more) {
            if (!part.isEmpty()) {
                sb.append(PATH_SEPARATOR).append(part);
            }
        }
        return new DatasetPath(this, sb.toString());
    }

    @Override
    public PathMatcher getPathMatcher(String syntaxAndPattern) {
        throw new UnsupportedOperationException("Dataset filesystem does not support path matching");
    }

    @Override
    public UserPrincipalLookupService getUserPrincipalLookupService() {
        throw new UnsupportedOperationException();
    }

    @Override
    public WatchService newWatchService() throws IOException {
        throw new UnsupportedOperationException();
    }

    Map<String, ?> getEnv() {
        return env;
    }
}
