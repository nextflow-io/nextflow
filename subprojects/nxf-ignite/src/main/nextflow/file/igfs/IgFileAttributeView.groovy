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

import java.nio.file.NoSuchFileException
import java.nio.file.attribute.BasicFileAttributeView
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

import groovy.transform.CompileStatic

/**
 * A file attribute view that provides a view of a <em>basic set</em> of file
 * attributes common to many file systems. The basic set of file attributes
 * consist of <em>mandatory</em> and <em>optional</em> file attributes as
 * defined by the {@link BasicFileAttributes} interface.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class IgFileAttributeView implements BasicFileAttributeView {

    static String NAME = "basic";

    private IgPath path;

    public IgFileAttributeView(IgPath path) {
        this.path = path
    }

    /**
     * @Inheritdoc
     */
    @Override
    public String name() {
         return NAME;
    }

    /**
     * @Inheritdoc
     */
    @Override
    public BasicFileAttributes readAttributes() throws IOException {
        def attrs = path.nativeReadAttributes()
        if( attrs )
            return new IgFileAttributes(attrs)
        else
            throw new NoSuchFileException("Cannot read attributes for file: $path")
    }


    /**
     * @Inheritdoc
     *
     * NOTE: Creation time is not supported by the Ignite file system, this attribute is just ignored
     *
     * @param lastModifiedTime
     * @param lastAccessTime
     * @param createTime
     * @throws IOException
     */
    @Override
    public void setTimes(FileTime lastModifiedTime, FileTime lastAccessTime, FileTime createTime) throws IOException {
        long lastAccess = lastAccessTime?.toMillis() ?: 0
        long lastModified = lastModifiedTime?.toMillis() ?: 0
        path.nativeSetTime( lastAccess, lastModified )
    }

}