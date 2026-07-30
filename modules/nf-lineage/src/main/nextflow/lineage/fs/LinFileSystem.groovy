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

package nextflow.lineage.fs

import com.google.common.collect.ImmutableSet
import nextflow.lineage.LinStore
import nextflow.lineage.LinStoreFactory

import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.WatchService
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import nextflow.lineage.LinStore
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.config.LineageConfig

/**
 * File system for LID Paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class LinFileSystem extends FileSystem {

    private LinFileSystemProvider provider

    private LinStore store

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected LinFileSystem(){}

    LinFileSystem(LinFileSystemProvider provider, LineageConfig config) {
        this.provider = provider
        this.store = LinStoreFactory.create(config)
    }

    LinStore getStore() {
        return store
    }

    @Override
    boolean equals( Object other ) {
        if( this.class != other.class ) return false
        final that = (LinFileSystem)other
        this.provider == that.provider && this.store == that.store
    }

    @Override
    int hashCode() {
        Objects.hash(provider,store)
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
        return LinPath.SEPARATOR
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
        return ImmutableSet.of("basic")
    }

    @Override
    Path getPath(String first, String... more) {
        final path = more ? LinPath.SEPARATOR + more.join(LinPath. SEPARATOR) : ''
        return getPath(LinPath.asUri(LinPath.LID_PROT + first + path))
    }

    Path getPath(URI uri){
        return new LinPath(this, uri)
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
