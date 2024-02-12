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

import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.LinkOption
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class GsFileSystemProvider extends FileSystemProvider {

    @Delegate
    CloudStorageFileSystemProvider delegate

    GsFileSystemProvider() {
        this(new CloudStorageFileSystemProvider())
    }

    GsFileSystemProvider(CloudStorageFileSystemProvider delegate) {
        this.delegate = delegate
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) {
        delegate.checkAccess(unwrap0(path), modes)
    }

    @Override
    void copy(Path source, Path target, CopyOption... options) {
        delegate.copy(unwrap0(source), unwrap0(target), options)
    }

    @Override
    void createDirectory(Path dir, FileAttribute<?>... attrs) {
        delegate.createDirectory(unwrap0(dir), attrs)
    }

    @Override
    void delete(Path path) {
        delegate.delete(unwrap0(path))
    }

    @Override
    boolean deleteIfExists(Path path) {
        delegate.deleteIfExists(unwrap0(path))
    }

    @Override
    <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        delegate.getFileAttributeView(unwrap0(path), type, options)
    }

    @Override
    boolean isHidden(Path path) {
        delegate.isHidden(unwrap0(path))
    }

    @Override
    boolean isSameFile(Path path, Path path2) {
        delegate.isSameFile(unwrap0(path), unwrap0(path2))
    }

    @Override
    void move(Path source, Path target, CopyOption... options) {
        delegate.move(unwrap0(source), unwrap0(target), options)
    }

    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) {
        final directoryStream = delegate.newDirectoryStream(unwrap0(dir), filter)
        final iterator = new GsPathIterator(directoryStream.iterator())
        return new DirectoryStream<Path>() {
            @Override
            Iterator<Path> iterator() { iterator }

            @Override
            void close() {}
        }
    }

    @Override
    GsFileSystem newFileSystem(URI uri, Map<String,?> env) {
        new GsFileSystem(delegate.newFileSystem(uri, env))
    }

    @Override
    InputStream newInputStream(Path path, OpenOption... options) {
        delegate.newInputStream(unwrap0(path), options)
    }

    @Override
    OutputStream newOutputStream(Path path, OpenOption... options) {
        delegate.newOutputStream(unwrap0(path), options)
    }

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) {
        delegate.readAttributes(unwrap0(path), type, options)
    }

    private static Path unwrap0(Path path) {
        path instanceof GsPath ? path.target : path
    }

    private static class GsPathIterator implements Iterator<Path> {
        @Delegate
        private Iterator<Path> delegate

        GsPathIterator(Iterator<Path> delegate) {
            this.delegate = delegate
        }

        @Override
        Path next() {
            final path = delegate.next()
            path instanceof CloudStoragePath ? new GsPath(path) : path
        }
    }

}
