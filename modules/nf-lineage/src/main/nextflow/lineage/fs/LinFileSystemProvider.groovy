/*
 * Copyright 2013-2025, Seqera Labs
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

import java.nio.ByteBuffer
import java.nio.channels.NonWritableChannelException
import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessDeniedException
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
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import nextflow.lineage.config.LineageConfig
import nextflow.util.TestOnly
/**
 * File System Provider for LID Paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class LinFileSystemProvider extends FileSystemProvider {

    public static final String SCHEME = "lid"

    private LinFileSystem fileSystem

    @Override
    String getScheme() {
        return SCHEME
    }

    protected LinPath toLinPath(Path path) {
        if (path !instanceof LinPath)
            throw new ProviderMismatchException()
        if (path instanceof LinMetadataPath)
            return (LinMetadataPath) path
        return (LinPath) path
    }

    private void checkScheme(URI uri) {
        final scheme = uri.scheme.toLowerCase()
        if (scheme != getScheme())
            throw new IllegalArgumentException("Not a valid ${getScheme().toUpperCase()} scheme: $scheme")
    }

    @Override
    synchronized FileSystem newFileSystem(URI uri, Map<String, ?> config) throws IOException {
        checkScheme(uri)
        if (fileSystem) {
            return fileSystem
        }
        //Overwrite default values with provided configuration
        final defaultConfig = LineageConfig.asMap()
        if (config) {
            for (Map.Entry<String, ?> e : config.entrySet()) {
                defaultConfig.put(e.key, e.value)
            }
        }
        return fileSystem = new LinFileSystem(this, new LineageConfig(defaultConfig))
    }

    @Override
    FileSystem getFileSystem(URI uri) throws FileSystemNotFoundException {
        if (!fileSystem)
            throw new FileSystemNotFoundException()
        return fileSystem
    }

    synchronized FileSystem getFileSystemOrCreate(URI uri) {
        checkScheme(uri)
        if (!fileSystem) {
            fileSystem = (LinFileSystem) newFileSystem(uri, LineageConfig.asMap())
        }
        return fileSystem
    }

    @Override
    LinPath getPath(URI uri) {
        return (LinPath) ((LinFileSystem) getFileSystemOrCreate(uri)).getPath(uri)
    }

    @Override
    OutputStream newOutputStream(Path path, OpenOption... options) throws IOException {
        throw new UnsupportedOperationException("Write not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final lid = toLinPath(path)
        if (lid instanceof LinMetadataPath)
            return (lid as LinMetadataPath).newInputStream()
        return newInputStream0(lid, options)
    }

    private static InputStream newInputStream0(LinPath lid, OpenOption... options) throws IOException {
        final realPath = lid.getTargetOrMetadataPath()
        if (realPath instanceof LinMetadataPath)
            return (realPath as LinMetadataPath).newInputStream()
        return realPath.fileSystem.provider().newInputStream(realPath, options)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        final lid = toLinPath(path)
        validateOptions(options)
        return newByteChannel0(lid, options, attrs)
    }

    @CompileStatic
    private class LinPathSeekableByteChannel implements SeekableByteChannel {
        SeekableByteChannel channel

        LinPathSeekableByteChannel(SeekableByteChannel channel) {
            this.channel = channel
        }

        @Override
        int read(ByteBuffer dst) throws IOException {
            return channel.read(dst)
        }

        @Override
        int write(ByteBuffer src) throws IOException {
            throw new NonWritableChannelException(){}
        }

        @Override
        long position() throws IOException {
            return channel.position()
        }

        @Override
        SeekableByteChannel position(long newPosition) throws IOException {
            channel.position(newPosition)
            return this
        }

        @Override
        long size() throws IOException {
            return channel.size()
        }

        @Override
        SeekableByteChannel truncate(long unused) throws IOException {
            throw new NonWritableChannelException()
        }

        @Override
        boolean isOpen() {
            return channel.isOpen()
        }

        @Override
        void close() throws IOException {
            channel.close()
        }
    }

    private static void validateOptions(Set<? extends OpenOption> options) {
        if (!options || options.empty)
            return
        for (OpenOption opt : options) {
            // All OpenOption values except for APPEND and WRITE are allowed
            if (opt == StandardOpenOption.APPEND || opt == StandardOpenOption.WRITE)
                throw new UnsupportedOperationException("'$opt' not allowed");
        }

    }

    private SeekableByteChannel newByteChannel0(LinPath lid, Set<? extends OpenOption> options, FileAttribute<?>... attrs) {
        if (lid instanceof LinMetadataPath) {
            return (lid as LinMetadataPath).newSeekableByteChannel()
        }
        final realPath = lid.getTargetOrMetadataPath()
        if (realPath instanceof LinMetadataPath) {
            return (realPath as LinMetadataPath).newSeekableByteChannel()
        } else {
            SeekableByteChannel channel = realPath.fileSystem.provider().newByteChannel(realPath, options, attrs)
            return new LinPathSeekableByteChannel(channel)
        }
    }

    @Override
    DirectoryStream<Path> newDirectoryStream(Path path, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final lid = toLinPath(path)
        final real = lid.getTargetOrIntermediatePath()
        if (real instanceof LinIntermediatePath)
            return getDirectoryStreamFromSubPath(lid)
        return getDirectoryStreamFromRealPath(real, lid)
    }

    private static DirectoryStream<Path> getDirectoryStreamFromSubPath(LinPath lid){
        final paths = lid.getSubPaths()
        if( !paths )
            throw new FileNotFoundException("Sub paths for '$lid' do not exist")
        return new DirectoryStream<Path>() {
            Iterator<Path> iterator() {
                return paths.iterator()
            }

            void close() {
                paths.close()
            }
        }
    }
    private DirectoryStream<Path> getDirectoryStreamFromRealPath(Path real, LinPath lid) {
        final stream = real
            .getFileSystem()
            .provider()
            .newDirectoryStream(real, new LidFilter(fileSystem))

        return new DirectoryStream<Path>() {

            @Override
            Iterator<Path> iterator() {
                return new LidIterator(fileSystem, stream.iterator(), lid, real)
            }

            @Override
            void close() throws IOException {
                stream.close()
            }
        }
    }

    @CompileStatic
    private class LidFilter implements DirectoryStream.Filter<Path> {

        private final LinFileSystem fs

        LidFilter(LinFileSystem fs) {
            this.fs = fs
        }

        @Override
        boolean accept(Path entry) throws IOException {
            return true
        }
    }

    private static LinPath fromRealToLinPath(Path toConvert, Path realBase, LinPath lidBase) {
        if (toConvert.isAbsolute()) {
            if (toConvert.class != realBase.class) {
                throw new ProviderMismatchException()
            }
            final relative = realBase.relativize(toConvert)
            return (LinPath) lidBase.resolve(relative.toString())
        } else {
            return (LinPath) lidBase.resolve(toConvert.toString())
        }
    }

    private static class LidIterator implements Iterator<Path> {

        private final LinFileSystem fs
        private final Iterator<Path> target
        private final LinPath parent
        private final Path parentReal

        LidIterator(LinFileSystem fs, Iterator<Path> itr, LinPath parent, Path real) {
            this.fs = fs
            this.target = itr
            this.parent = parent
            this.parentReal = real
        }

        @Override
        boolean hasNext() {
            return target.hasNext()
        }

        @Override
        LinPath next() {
            final path = target.next()
            return path ? fromRealToLinPath(path, parentReal, parent) : null
        }
    }

    @Override
    void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        throw new UnsupportedOperationException("Create directory not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void delete(Path path) throws IOException {
        throw new UnsupportedOperationException("Delete not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void copy(Path source, Path target, CopyOption... options) throws IOException {
        throw new UnsupportedOperationException("Copy not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void move(Path source, Path target, CopyOption... options) throws IOException {
        throw new UnsupportedOperationException("Move not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    boolean isSameFile(Path path, Path path2) throws IOException {
        return path == path2
    }

    @Override
    boolean isHidden(Path path) throws IOException {
        return toLinPath(path).getTargetOrMetadataPath().isHidden()
    }

    @Override
    FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException("File store not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        validateAccessModes(modes)
        final lid = toLinPath(path)
        if (lid instanceof LinMetadataPath)
            return
        checkAccess0(lid, modes)
    }

    private void checkAccess0(LinPath lid, AccessMode... modes) {
        final real = lid.getTargetOrMetadataPath()
        if (real instanceof LinMetadataPath)
            return
        real.fileSystem.provider().checkAccess(real, modes)
    }

    private void validateAccessModes(AccessMode... modes) {
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE)
                throw new AccessDeniedException("Write mode not supported")
            if (m == AccessMode.EXECUTE)
                throw new AccessDeniedException("Execute mode not supported")
        }
    }

    @Override
    <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        return null
    }

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        final lid = toLinPath(path)
        if (lid instanceof LinMetadataPath)
            return (lid as LinMetadataPath).readAttributes(type)
        return readAttributes0(lid, type, options)
    }

    private <A extends BasicFileAttributes> A readAttributes0(LinPath lid, Class<A> type, LinkOption... options) throws IOException {
        final target = lid.getTargetOrIntermediatePath()
        if (target instanceof LinIntermediatePath)
            return (target as LinIntermediatePath).readAttributes(type)
        return target.fileSystem.provider().readAttributes(target, type, options)
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Read file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Set file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @TestOnly
    void reset() {
        fileSystem=null
    }
}
