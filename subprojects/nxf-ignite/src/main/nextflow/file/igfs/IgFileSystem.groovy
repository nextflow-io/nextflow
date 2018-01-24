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

package nextflow.file.igfs
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
import nextflow.file.FileHelper
import org.apache.ignite.IgniteFileSystem
/**
 * Implements a {@code FileSystem} for the Ignite provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class IgFileSystem extends FileSystem {

    private IgFileSystemProvider provider

    /**
     * The underlying {@link IgniteFileSystem} instance
     */
    private IgniteFileSystem igfs

    IgFileSystem( IgFileSystemProvider provider, IgniteFileSystem igfs ) {
        this.provider = provider
        this.igfs = igfs
    }

    /**
     * @return The underlying {@link IgniteFileSystem} instance providing Ignite platform in-memory storage
     */
    @PackageScope
    IgniteFileSystem getIgfs() { igfs }


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

    }

    /**
     * @Inheritdoc
     */
    @Override
    boolean isOpen() {
        return igfs != null
    }

    /**
     * @return {@code false}
     */
    @Override
    boolean isReadOnly() {
        return false
    }

    /**
     * @return {@link IgPath#PATH_SEPARATOR}
     */
    @Override
    String getSeparator() {
        return IgPath.PATH_SEPARATOR
    }

    @Override
    Iterable<Path> getRootDirectories() {
        ImmutableList.of((Path) new IgPath(this, IgPath.PATH_SEPARATOR) )
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
    public IgPath getPath(String first, String... more) {
        if (more.length == 0) {
            return new IgPath(this, first);
        }

        new IgPath(this, first, more);
    }

    /**
     *
     * @Inheritdoc
     */
    @Override
    PathMatcher getPathMatcher(String syntaxAndInput) {
        FileHelper.getDefaultPathMatcher(syntaxAndInput)
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
