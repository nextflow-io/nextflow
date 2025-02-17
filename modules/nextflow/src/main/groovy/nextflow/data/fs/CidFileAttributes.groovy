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

package nextflow.data.fs

import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

/**
 * Model CID file attributes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CidFileAttributes implements BasicFileAttributes {

    private final CidFileId fileId

    private final boolean directory

    private final long size

    private final FileTime creationTime

    private final FileTime lastModifiedTime

    private final FileTime lastAccessTime

    CidFileAttributes(CidFileId fileId, BasicFileAttributes attrs) {
        this.fileId = fileId
        this.directory = attrs.isDirectory()
        this.size = attrs.size()
        this.creationTime = attrs.creationTime()
        this.lastModifiedTime = attrs.lastModifiedTime()
        this.lastAccessTime = attrs.lastAccessTime()
    }

    @Override
    FileTime lastModifiedTime() {
        return lastModifiedTime
    }

    @Override
    FileTime lastAccessTime() {
        return lastAccessTime
    }

    @Override
    FileTime creationTime() {
        return creationTime
    }

    @Override
    boolean isRegularFile() {
        return !directory
    }

    @Override
    boolean isDirectory() {
        return directory
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
        return fileId
    }
}
