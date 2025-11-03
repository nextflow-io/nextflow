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

import com.google.common.collect.ImmutableSet;
import io.seqera.http.HxClient;
import io.seqera.tower.plugin.datalink.DataLink;

import java.io.IOException;
import java.net.URI;
import java.nio.file.FileStore;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.WatchService;
import java.nio.file.attribute.UserPrincipalLookupService;
import java.nio.file.spi.FileSystemProvider;
import java.util.Collections;
import java.util.Objects;
import java.util.Set;

/**
 * File system for Seqera Platform Data-Link Paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraFileSystem extends FileSystem {

    private final SeqeraFileSystemProvider provider;
    private final String endpoint;
    private final String workspaceId;
    private final HxClient httpClient;
    private final DataLink dataLink;

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected SeqeraFileSystem() {
        this.provider = null;
        this.endpoint = null;
        this.workspaceId = null;
        this.httpClient = null;
        this.dataLink = null;
    }

    public SeqeraFileSystem(SeqeraFileSystemProvider provider, DataLink dataLink, String endpoint, String workspaceId, HxClient httpClient) {
        this.provider = provider;
        this.endpoint = endpoint;
        this.workspaceId = workspaceId;
        this.httpClient = httpClient;
        this.dataLink = dataLink;
    }

    public String getEndpoint() {
        return endpoint;
    }

    public String getWorkspaceId() {
        return workspaceId;
    }

    public DataLink getDataLink() {
        return dataLink;
    }

    public String getDataLinkName() {
        return dataLink != null ? dataLink.getName() : null;
    }

    public HxClient getHttpClient() {
        return httpClient;
    }

    @Override
    public boolean equals(Object other) {
        if (this.getClass() != other.getClass()) {
            return false;
        }
        final SeqeraFileSystem that = (SeqeraFileSystem) other;
        return Objects.equals(this.provider, that.provider) &&
               Objects.equals(this.dataLink.getId(), that.dataLink.getId()) &&
               Objects.equals(this.endpoint, that.endpoint);
    }

    @Override
    public int hashCode() {
        return Objects.hash(provider, dataLink.getId(), endpoint);
    }

    @Override
    public FileSystemProvider provider() {
        return provider;
    }

    @Override
    public void close() throws IOException {
        // Cleanup if needed
    }

    @Override
    public boolean isOpen() {
        return true;
    }

    @Override
    public boolean isReadOnly() {
        return false;  // Allow both read and write operations
    }

    @Override
    public String getSeparator() {
        return SeqeraPath.SEPARATOR;
    }

    @Override
    public Iterable<Path> getRootDirectories() {
        return Collections.singletonList(getPath("/"));
    }

    @Override
    public Iterable<FileStore> getFileStores() {
        return null;
    }

    @Override
    public Set<String> supportedFileAttributeViews() {
        return ImmutableSet.of("basic");
    }

    @Override
    public Path getPath(String first, String... more) {
        final StringBuilder pathBuilder = new StringBuilder(first);
        if (more != null && more.length > 0) {
            for (String segment : more) {
                pathBuilder.append(SeqeraPath.SEPARATOR).append(segment);
            }
        }
        return new SeqeraPath(this, pathBuilder.toString());
    }

    public Path getPath(URI uri) {
        return new SeqeraPath(this, uri);
    }

    @Override
    public PathMatcher getPathMatcher(String syntaxAndPattern) {
        throw new UnsupportedOperationException("Path matcher not supported");
    }

    @Override
    public UserPrincipalLookupService getUserPrincipalLookupService() {
        throw new UnsupportedOperationException("User Principal Lookup Service not supported");
    }

    @Override
    public WatchService newWatchService() throws IOException {
        throw new UnsupportedOperationException("Watch Service not supported");
    }
}
