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

import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

import groovy.transform.CompileStatic
import groovy.transform.PackageScope

/**
 * Implements a {@link BasicFileAttributes} for http/ftp
 * file system providers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
@PackageScope
@CompileStatic
class XFileAttributes implements BasicFileAttributes {

    private FileTime lastModifiedTime
    private long size

    XFileAttributes(FileTime lastModifiedTime, long size) {
        this.lastModifiedTime = lastModifiedTime
        this.size = size
    }

    @Override
    FileTime lastModifiedTime() {
        return lastModifiedTime
    }

    @Override
    FileTime lastAccessTime() {
        return null
    }

    @Override
    FileTime creationTime() {
        return null
    }

    @Override
    boolean isRegularFile() {
        return true
    }

    @Override
    boolean isDirectory() {
        return false
    }

    @Override
    boolean isSymbolicLink() {
        return false
    }

    @Override
    boolean isOther() {
        return false
    }

    @Override
    long size() {
        return size
    }

    @Override
    Object fileKey() {
        return null
    }
}
