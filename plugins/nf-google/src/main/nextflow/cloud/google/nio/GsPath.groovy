/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cloud.google.nio

import java.nio.file.LinkOption
import java.nio.file.Path

import com.google.cloud.storage.Blob
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import nextflow.file.ETagAwareFile
/**
 * Implements an ETag-aware wrapper around the NIO filesystem
 * provided by the Google Cloud SDK.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class GsPath implements Path, ETagAwareFile {

    @Delegate
    CloudStoragePath target

    private GsFileSystem fileSystem

    private Blob metadata

    GsPath(CloudStoragePath target) {
        this(new GsFileSystem(target.fileSystem), target)
    }

    GsPath(GsFileSystem fileSystem, CloudStoragePath target) {
        this.target = target
        this.fileSystem = fileSystem
    }

    @Override
    int compareTo(Path other) {
        return target.compareTo(unwrap0(other))
    }

    @Override
    boolean equals(Object other) {
        return target.equals(other instanceof GsPath ? other.target : other)
    }

    @Override
    Path getFileName() {
        return wrap0(target.getFileName())
    }

    @Override
    Path getName(int index) {
        return wrap0(target.getName(index))
    }

    @Override
    Path getParent() {
        return wrap0(target.getParent())
    }

    @Override
    Path getRoot() {
        return wrap0(target.getRoot())
    }

    @Override
    GsFileSystem getFileSystem() {
        return fileSystem
    }

    @Override
    int hashCode() {
        return target.hashCode()
    }

    @Override
    Path normalize() {
        return wrap0(target.normalize())
    }

    @Override
    Path relativize(Path other) {
        return wrap0(target.relativize(unwrap0(other)))
    }

    @Override
    Path resolve(String other) {
        return wrap0(target.resolve(other))
    }

    @Override
    Path resolve(Path other) {
        return wrap0(target.resolve(unwrap0(other)))
    }

    @Override
    Path resolveSibling(String other) {
        return wrap0(target.resolveSibling(other))
    }

    @Override
    Path resolveSibling(Path other) {
        return wrap0(target.resolveSibling(unwrap0(other)))
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        return wrap0(target.subpath(beginIndex, endIndex))
    }

    @Override
    Path toAbsolutePath() {
        return wrap0(target.toAbsolutePath())
    }

    @Override
    Path toRealPath(LinkOption... options) {
        return wrap0(target.toRealPath(options))
    }

    @Override
    String toString() {
        return target.toString()
    }

    @Override
    String getETag() {
        if( metadata == null ) {
            final googleOpts = GoogleOpts.create(Global.session as Session)
            final client = StorageOptions.newBuilder()
                    .setProjectId(googleOpts.projectId)
                    .build()
                    .getService()
            metadata = client.get(target.bucket(), target.toString())
        }
        return metadata?.getMd5ToHexString()
    }

    private static Path unwrap0(Path path) {
        path instanceof GsPath ? path.target : path
    }

    private Path wrap0(Path path) {
        path instanceof CloudStoragePath ? new GsPath(fileSystem, path) : path
    }

}
