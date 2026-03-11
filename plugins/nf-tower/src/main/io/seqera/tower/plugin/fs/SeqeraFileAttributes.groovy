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
import io.seqera.tower.plugin.dataset.DatasetDto

/**
 * {@link BasicFileAttributes} for {@code seqera://} paths.
 * For depth &lt; 4 (directory paths): {@code isDirectory=true}, {@code size=0}.
 * For depth 4 (dataset file paths): {@code isRegularFile=true}, timestamps from {@link DatasetDto}.
 *
 * @author Seqera Labs
 */
@CompileStatic
class SeqeraFileAttributes implements BasicFileAttributes {

    private final boolean directory
    private final DatasetDto dataset

    /** Construct attributes for a virtual directory (depth 0–3). */
    SeqeraFileAttributes(boolean isDir) {
        this.directory = isDir
        this.dataset = null
    }

    /** Construct attributes for a dataset file (depth 4). */
    SeqeraFileAttributes(DatasetDto dataset) {
        this.directory = false
        this.dataset = dataset
    }

    @Override
    FileTime lastModifiedTime() {
        if (dataset?.lastUpdated) {
            try {
                return FileTime.from(Instant.parse(dataset.lastUpdated))
            } catch (Exception e) {
                // fall through to epoch
            }
        }
        return FileTime.from(Instant.EPOCH)
    }

    @Override
    FileTime lastAccessTime() { lastModifiedTime() }

    @Override
    FileTime creationTime() {
        if (dataset?.dateCreated) {
            try {
                return FileTime.from(Instant.parse(dataset.dateCreated))
            } catch (Exception e) {
                // fall through to epoch
            }
        }
        return FileTime.from(Instant.EPOCH)
    }

    @Override
    boolean isRegularFile() { !directory }

    @Override
    boolean isDirectory() { directory }

    @Override
    boolean isSymbolicLink() { false }

    @Override
    boolean isOther() { false }

    @Override
    long size() { 0L }

    @Override
    Object fileKey() { dataset?.id }
}
