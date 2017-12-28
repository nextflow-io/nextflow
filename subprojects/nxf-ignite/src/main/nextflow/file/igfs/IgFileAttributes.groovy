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
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import org.apache.ignite.igfs.IgfsFile
import org.apache.ignite.internal.processors.igfs.IgfsFileImpl
/**
 * Implements the basic attributes object for Ignite file system provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@EqualsAndHashCode
@ToString
class IgFileAttributes implements BasicFileAttributes {

    private final FileTime lastModifiedTime;
    private final long size;
    private final boolean directory;
    private final boolean regularFile;
    private final key;

    protected IgFileAttributes( key, FileTime lastModifiedTime, long size, boolean isDirectory, boolean isRegularFile) {
        this.key = key;
        this.lastModifiedTime = lastModifiedTime;
        this.size = size;
        directory = isDirectory;
        regularFile = isRegularFile;
    }

    IgFileAttributes( IgfsFile fileInfo ) {
        this.lastModifiedTime = FileTime.fromMillis(fileInfo.modificationTime() )
        this.size = fileInfo.length()
        this.directory = fileInfo.isDirectory()
        this.regularFile = fileInfo.isFile()

        if( fileInfo instanceof IgfsFileImpl )
            this.key = (fileInfo as IgfsFileImpl).fileId()
    }

    /**
     * @Inheritdoc
     */
    @Override
    FileTime lastModifiedTime() {
        return lastModifiedTime;
    }

    /**
     * @Inheritdoc
     */
    @Override
    FileTime lastAccessTime() {
        return lastModifiedTime;
    }

    /**
     * Creation time not supported by Ignite file system provider
     */
    @Override
    FileTime creationTime() {
        null
    }

    /**
     * @Inheritdoc
     */
    @Override
    boolean isRegularFile() {
        regularFile;
    }

    /**
     * @Inheritdoc
     */
    @Override
    public boolean isDirectory() {
        directory;
    }

    /**
     * Sym link are not supported by Ignite file system provider
     *
     * @return {@code false} by definition
     */
    @Override
    public boolean isSymbolicLink() {
        false;
    }

    /**
     * Other group is not supported by Ignite file system provider
     *
     * @return {@code false} by definition
     */

    @Override
    public boolean isOther() {
        false;
    }

    /**
     * @Inheritdoc
     */
    @Override
    public long size() {
        size;
    }

    /**
     * @Inheritdoc
     */
    @Override
    public Object fileKey() {
        key;
    }

}
