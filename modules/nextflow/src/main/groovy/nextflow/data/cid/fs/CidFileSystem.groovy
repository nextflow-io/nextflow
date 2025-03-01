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
 *
 */

package nextflow.data.cid.fs

import nextflow.data.cid.DefaultCidStore

import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.WatchService
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import nextflow.data.config.DataConfig

/**
 * File system for CID Paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class CidFileSystem extends FileSystem {

    private CidFileSystemProvider provider

    private Path basePath

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected CidFileSystem(){}

    CidFileSystem(CidFileSystemProvider provider, DataConfig config) {
        this.provider = provider
        this.basePath = DefaultCidStore.getMetadataPath(config)
    }

    Path getBasePath() {
        return basePath
    }

    @Override
    boolean equals( Object other ) {
        if( this.class != other.class ) return false
        final that = (CidFileSystem)other
        this.provider == that.provider && this.basePath == that.basePath
    }

    @Override
    int hashCode() {
        Objects.hash(provider,basePath)
    }

    @Override
    FileSystemProvider provider() {
        return provider
    }

    @Override
    void close() throws IOException {

    }

    @Override
    boolean isOpen() {
        return false
    }

    @Override
    boolean isReadOnly() {
        return true
    }

    @Override
    String getSeparator() {
        return CidPath.SEPARATOR
    }

    @Override
    Iterable<Path> getRootDirectories() {
        return null
    }

    @Override
    Iterable<FileStore> getFileStores() {
        return null
    }

    @Override
    Set<String> supportedFileAttributeViews() {
        return null
    }

    @Override
    Path getPath(String first, String... more) {
        return new CidPath(this,first,more)
    }

    @Override
    PathMatcher getPathMatcher(String syntaxAndPattern) {
        throw new UnsupportedOperationException();
    }

    @Override
    UserPrincipalLookupService getUserPrincipalLookupService() {
        throw new UnsupportedOperationException('User Principal Lookup Service not supported')
    }

    @Override
    WatchService newWatchService() throws IOException {
        throw new UnsupportedOperationException('Watch Service not supported')
    }
}
