/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

import groovy.transform.CompileStatic
import groovy.transform.PackageScope

/**
 * Implements a {@link BasicFileAttributes} for http/ftp
 * file system providers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
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
