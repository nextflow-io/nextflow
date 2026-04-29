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

package io.seqera.tower.plugin.fs

import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime
import java.time.Instant

import groovy.transform.CompileStatic

/**
 * {@link BasicFileAttributes} for {@code seqera://} paths.
 *
 * Resource-type agnostic: virtual directories use the {@code (boolean isDir)}
 * constructor; file-like entries use the explicit {@code (size, lastMod, created, key)}
 * constructor. Handlers build instances using whatever metadata the underlying
 * resource exposes.
 */
@CompileStatic
class SeqeraFileAttributes implements BasicFileAttributes {

    private final boolean directory
    private final long size
    private final Instant lastModified
    private final Instant created
    private final Object fileKey

    /** Construct attributes for a virtual directory. */
    SeqeraFileAttributes(boolean isDir) {
        this.directory = isDir
        this.size = 0L
        this.lastModified = Instant.EPOCH
        this.created = Instant.EPOCH
        this.fileKey = null
    }

    /** Construct attributes for a regular file with explicit metadata. */
    SeqeraFileAttributes(long size, Instant lastModified, Instant created, Object fileKey) {
        this.directory = false
        this.size = size >= 0 ? size : 0L
        this.lastModified = lastModified ?: Instant.EPOCH
        this.created = created ?: Instant.EPOCH
        this.fileKey = fileKey
    }

    @Override FileTime lastModifiedTime() { FileTime.from(lastModified) }

    @Override FileTime lastAccessTime() { FileTime.from(lastModified) }

    @Override FileTime creationTime() { FileTime.from(created) }

    @Override boolean isRegularFile() { !directory }

    @Override boolean isDirectory() { directory }

    @Override boolean isSymbolicLink() { false }

    @Override boolean isOther() { false }

    @Override long size() { size }

    @Override Object fileKey() { fileKey }
}
