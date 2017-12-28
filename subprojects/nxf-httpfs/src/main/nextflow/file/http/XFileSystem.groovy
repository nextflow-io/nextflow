/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
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
