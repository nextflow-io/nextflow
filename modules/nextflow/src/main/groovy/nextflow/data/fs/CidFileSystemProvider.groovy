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

import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemNotFoundException
import java.nio.file.LinkOption
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider
import java.util.function.BiFunction

import groovy.transform.CompileStatic
import nextflow.data.config.DataConfig
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CidFileSystemProvider extends FileSystemProvider {

    public static final String SCHEME = "cid"

    private CidFileSystem fileSystem

    @Override
    String getScheme() {
        return SCHEME
    }

    protected CidPath toCidPath(Path path) {
        if (path !instanceof CidPath)
            throw new ProviderMismatchException()
        return (CidPath) path
    }

    private void checkScheme(URI uri) {
        final scheme = uri.scheme.toLowerCase()
        if( scheme != getScheme() )
            throw new IllegalArgumentException("Not a valid ${getScheme().toUpperCase()} scheme: $scheme")
    }

    @Override
    synchronized FileSystem newFileSystem(URI uri, Map<String, ?> config) throws IOException {
        checkScheme(uri)
        if( !fileSystem ) {
            fileSystem = new CidFileSystem(this, new DataConfig(config))
        }
        return fileSystem
    }

    @Override
    FileSystem getFileSystem(URI uri) throws FileSystemNotFoundException {
        checkScheme(uri)
        if( fileSystem==null )
            throw new FileSystemNotFoundException("File system has not been created yet - offending URI: $uri")
        return fileSystem
    }

    synchronized FileSystem getFileSystemOrCreate(URI uri) {
        checkScheme(uri)
        if( !fileSystem ) {
            fileSystem = (CidFileSystem) newFileSystem(uri, DataConfig.asMap())
        }
        return fileSystem
    }

    @Override
    CidPath getPath(URI uri) {
        // the URI authority holds the base component of the CID path
        final base = uri.authority
        final path = uri.path
        return (CidPath) getFileSystemOrCreate(uri).getPath(base, path)
    }

    CidPath getPath(String uri) {
        getPath(CidPath.asUri(uri))
    }

    protected <R> R delegateTo(Path path, BiFunction<FileSystemProvider,Path,R> action) {
        final cid = toCidPath(path)
        action.apply(cid.getTargetSystem().provider(), cid.toRealPath())
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        return delegateTo(path, (provider, real)-> provider.newByteChannel(real, options, attrs))
    }

    @Override
    DirectoryStream<Path> newDirectoryStream(Path path, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final real = path.toRealPath()
        final stream = real
                .getFileSystem()
                .provider()
                .newDirectoryStream(real, (Path it) -> filter.accept(new CidPath(fileSystem,it)))

        return new DirectoryStream<Path>() {
            @Override
            Iterator<Path> iterator() {
                return new CidIterator(fileSystem, stream.iterator())
            }

            @Override
            void close() throws IOException {
                stream.close()
            }
        }
    }

    private static class CidIterator implements Iterator<Path> {

        private final CidFileSystem fs
        private final Iterator<Path> target

        CidIterator(CidFileSystem fs, Iterator<Path> itr) {
            this.fs = fs
            this.target = itr
        }

        @Override
        boolean hasNext() {
            return target.hasNext()
        }

        @Override
        CidPath next() {
            final path = target.next()
            return path ? new CidPath(fs,path) : null
        }
    }

    @Override
    void createDirectory(Path path, FileAttribute<?>... attrs) throws IOException {
        delegateTo(path, (provider, real)-> provider.createDirectory(real, attrs))
    }

    @Override
    void delete(Path path) throws IOException {
        delegateTo(path, (provider, real)-> provider.delete(real))
    }

    @Override
    void copy(Path source, Path target, CopyOption... options) throws IOException {
        final realTarget = toCidPath(target).toRealPath()
        delegateTo(source, (provider, realSource)-> provider.copy(realSource, realTarget, options))
    }

    @Override
    void move(Path source, Path target, CopyOption... options) throws IOException {
        final realTarget = toCidPath(target).toRealPath()
        delegateTo(source, (provider, realSource)-> provider.move(realSource, realTarget, options))
    }

    @Override
    boolean isSameFile(Path path, Path path2) throws IOException {
        return path == path2
    }

    @Override
    boolean isHidden(Path path) throws IOException {
        return toCidPath(path).toRealPath().isHidden()
    }

    @Override
    FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException("File store not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        delegateTo(path, (provider, real)-> provider.checkAccess(real, modes))
    }

    @Override
    <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        return null
    }

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        final cid = toCidPath(path)
        final attrs = delegateTo(cid, (provider, real)-> provider.readAttributes(real,type,options))
        return (A)new CidFileAttributes(cid.fileId, attrs)
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Read file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Set file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

}
