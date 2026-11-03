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

import groovy.transform.CompileStatic

import java.nio.channels.SeekableByteChannel
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

/**
 * Class to model the metadata descriptions as a file.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class LinMetadataPath extends LinPath {
    private byte[] results
    private FileTime creationTime

    LinMetadataPath(String resultsObject, FileTime creationTime, LinFileSystem fs, String path, String fragment) {
        super(fs, "${path}${fragment ? '#'+ fragment : ''}")
        this.results = resultsObject.getBytes("UTF-8")
        this.creationTime = creationTime
    }

    InputStream newInputStream() {
        return new ByteArrayInputStream(results)
    }

    SeekableByteChannel newSeekableByteChannel(){
        return new LinMetadataSeekableByteChannel(results)
    }

    <A extends BasicFileAttributes> A readAttributes(Class<A> type){
        return (A) new BasicFileAttributes() {
            @Override
            long size() { return results.length }

            @Override
            FileTime lastModifiedTime() { return creationTime }

            @Override
            FileTime lastAccessTime() { return creationTime }

            @Override
            FileTime creationTime() { return creationTime }

            @Override
            boolean isRegularFile() { return true }

            @Override
            boolean isDirectory() { return false }

            @Override
            boolean isSymbolicLink() { return false }

            @Override
            boolean isOther() { return false }

            @Override
            Object fileKey() { return null }
        }
    }
}
