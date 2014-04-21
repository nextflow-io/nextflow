/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.file.ggfs
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.WatchService
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider

import com.google.common.collect.ImmutableList
import com.google.common.collect.ImmutableSet
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import org.gridgain.grid.ggfs.GridGgfs
/**
 * Implements a {@code FileSystem} for the GridGain provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GgFileSystem extends FileSystem {

    private GgFileSystemProvider provider

    /**
     * The underlying {@code GridGgfs} instance
     */
    private GridGgfs ggfs

    GgFileSystem( GgFileSystemProvider provider, GridGgfs ggfs ) {
        this.provider = provider
        this.ggfs = ggfs
    }

    /**
     * @return The underlying {@code GridGgfs} instance providing GridGain platform in-memory storage
     */
    @PackageScope
    GridGgfs getGgfs() { ggfs }


    /**
     * @Inheritdoc
     */
    @Override
    FileSystemProvider provider() {
        return provider
    }

    /**
     * @Inheritdoc
     */
    @Override
    void close() throws IOException {
        ggfs.stop()
    }

    /**
     * @Inheritdoc
     */
    @Override
    boolean isOpen() {
        return ggfs != null
    }

    /**
     * @return {@code false}
     */
    @Override
    boolean isReadOnly() {
        return false
    }

    /**
     * @return {@code GgPath#PATH_SEPARATOR}
     */
    @Override
    String getSeparator() {
        return GgPath.PATH_SEPARATOR
    }

    @Override
    Iterable<Path> getRootDirectories() {
        ImmutableList.of( new GgPath(this, GgPath.PATH_SEPARATOR) )
    }

    @Override
    public Iterable<FileStore> getFileStores() {
        ImmutableList.of()
    }

    /**
     * @Inheritdoc
     */
    @Override
    public Set<String> supportedFileAttributeViews() {
        ImmutableSet.of("basic");
    }

    /**
     * @Inheritdoc
     */
    @Override
    public GgPath getPath(String first, String... more) {
        if (more.length == 0) {
            return new GgPath(this, first);
        }

        new GgPath(this, first, more);
    }

    /**
     * TODO
     * @param syntaxAndPattern
     * @return
     */
    @Override
    PathMatcher getPathMatcher(String syntaxAndPattern) {
        throw new UnsupportedOperationException("Method 'getUserPrincipalLookupService' not implemented yet")
    }

    /**
     * TODO
     * @return
     */
    @Override
    UserPrincipalLookupService getUserPrincipalLookupService() {
        throw new UnsupportedOperationException("Method 'getUserPrincipalLookupService' not supported")
    }

    /**
     * TODO
     * @return
     * @throws IOException
     */
    @Override
    WatchService newWatchService() throws IOException {
        throw new UnsupportedOperationException("Method 'newWatchService' not supported")
    }
}
