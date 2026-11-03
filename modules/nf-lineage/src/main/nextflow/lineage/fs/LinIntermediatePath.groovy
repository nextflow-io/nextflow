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

package nextflow.lineage.fs

import java.nio.channels.SeekableByteChannel
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime
import java.time.Instant

class LinIntermediatePath extends LinMetadataPath{

    LinIntermediatePath(LinFileSystem fs, String path) {
        super("", FileTime.from(Instant.now()), fs, path, null)
    }

    @Override
    InputStream newInputStream() {
        throw new UnsupportedOperationException()
    }
    @Override
    SeekableByteChannel newSeekableByteChannel(){
        throw new UnsupportedOperationException()
    }
    @Override
    <A extends BasicFileAttributes> A readAttributes(Class<A> type){
        return (A) new BasicFileAttributes() {
            @Override
            long size() { return 0 }

            @Override
            FileTime lastModifiedTime() { FileTime.from(Instant.now()) }

            @Override
            FileTime lastAccessTime() { FileTime.from(Instant.now()) }

            @Override
            FileTime creationTime() { FileTime.from(Instant.now()) }

            @Override
            boolean isRegularFile() { return false }

            @Override
            boolean isDirectory() { return true }

            @Override
            boolean isSymbolicLink() { return false }

            @Override
            boolean isOther() { return false }

            @Override
            Object fileKey() { return null }
        }
    }
}
