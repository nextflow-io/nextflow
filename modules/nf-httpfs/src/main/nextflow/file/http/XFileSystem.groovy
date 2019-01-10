/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.file.http

import groovy.transform.CompileStatic

import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.WatchService
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider

import groovy.transform.PackageScope

/**
 * Implements a read-only JSR-203 complaint file system for http/ftp protocols
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 */
@PackageScope
@CompileStatic
class XFileSystem extends FileSystem {

    static private String PATH_SEPARATOR = '/'

    private XFileSystemProvider provider

    private URI base

    XFileSystem(XFileSystemProvider provider, URI base) {
        this.provider = provider
        this.base = base
    }

    @Override
    boolean equals( Object other ) {
        if( this.class != other.class ) return false
        final that = (XFileSystem)other
        this.provider == that.provider && this.base == that.base
    }

    @Override
    int hashCode() {
        Objects.hash(provider,base)
    }

    @Override
    FileSystemProvider provider() {
        return provider
    }

    URI getBaseUri() {
        return base
    }

    @Override
    void close() throws IOException {

    }

    @Override
    boolean isOpen() {
        return true
    }

    @Override
    boolean isReadOnly() {
        return true
    }

    @Override
    String getSeparator() {
        return PATH_SEPARATOR
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
        return new XPath(this,first,more)
    }

    @Override
    PathMatcher getPathMatcher(String syntaxAndPattern) {
        return null
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
