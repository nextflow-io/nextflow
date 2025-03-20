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

package nextflow.data.cid.fs

import javax.swing.DefaultListSelectionModel
import java.nio.ByteBuffer
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
import nextflow.data.config.DataConfig

/**
 * File System Provider for CID Paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
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
            //Overwrite default values with provided configuration
            final defaultConfig = DataConfig.asMap()
            if (config) {
                for (Map.Entry<String,?> e : config.entrySet()) {
                    defaultConfig.put(e.key, e.value)
                }
            }
            fileSystem = new CidFileSystem(this, new DataConfig(defaultConfig))
        }
        return fileSystem
    }

    @Override
    FileSystem getFileSystem(URI uri) throws FileSystemNotFoundException {
        if (!fileSystem)
            throw new FileSystemNotFoundException()
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

    @Override
    OutputStream newOutputStream(Path path, OpenOption... options) throws IOException {
        throw new UnsupportedOperationException("Write not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final cid = toCidPath(path)
        final realPath = cid.getTargetPath(true)
        if (realPath instanceof CidResultsPath)
            return (realPath as CidResultsPath).newInputStream()
        else
            return realPath.fileSystem.provider().newInputStream(realPath, options)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        final cid = toCidPath(path)
        if (options.size() > 0) {
            for (OpenOption opt: options) {
                // All OpenOption values except for APPEND and WRITE are allowed
                if (opt == StandardOpenOption.APPEND || opt == StandardOpenOption.WRITE)
                    throw new UnsupportedOperationException("'$opt' not allowed");
            }
        }
        final realPath = cid.getTargetPath(true)
        SeekableByteChannel channel
        if (realPath instanceof CidResultsPath){
            channel = (realPath as CidResultsPath).newSeekableByteChannel()
        } else {
            channel = realPath.fileSystem.provider().newByteChannel(realPath, options, attrs)
        }

        new SeekableByteChannel() {

            @Override
            int read(ByteBuffer dst) throws IOException {
                return channel.read(dst)
            }

            @Override
            int write(ByteBuffer src) throws IOException {
                throw new UnsupportedOperationException("Write operation not supported")
            }

            @Override
            long position() throws IOException {
                return channel.position()
            }

            @Override
            SeekableByteChannel position(long newPosition) throws IOException {
                throw new UnsupportedOperationException("Position operation not supported")
            }

            @Override
            long size() throws IOException {
                return channel.size()
            }

            @Override
            SeekableByteChannel truncate(long unused) throws IOException {
                throw new UnsupportedOperationException("Truncate operation not supported")
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
    }

    @Override
    DirectoryStream<Path> newDirectoryStream(Path path, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final cid = toCidPath(path)
        final real = cid.getTargetPath(false)
        final stream = real
                .getFileSystem()
                .provider()
                .newDirectoryStream(real, new CidFilter(fileSystem))

        return new DirectoryStream<Path>() {

            @Override
            Iterator<Path> iterator() {
                return new CidIterator(fileSystem, stream.iterator(), cid, real)
            }

            @Override
            void close() throws IOException {
                stream.close()
            }
        }
    }
    private class CidFilter implements DirectoryStream.Filter<Path> {

        private final CidFileSystem fs

        CidFilter(CidFileSystem fs){
            this.fs = fs
        }

        @Override
        boolean accept(Path entry) throws IOException {
            return true
        }
    }

    private static CidPath fromRealToCidPath(Path toConvert, Path realBase, CidPath cidBase){
        if (toConvert.isAbsolute()) {
            if (toConvert.class != realBase.class){
                throw new ProviderMismatchException()
            }
            final relative = realBase.relativize(toConvert)
            return (CidPath) cidBase.resolve(relative.toString())
        } else {
            return (CidPath) cidBase.resolve(toConvert.toString())
        }
    }

    private static class CidIterator implements Iterator<Path> {

        private final CidFileSystem fs
        private final Iterator<Path> target
        private final CidPath parent
        private final Path parentReal

        CidIterator(CidFileSystem fs, Iterator<Path> itr, CidPath parent, Path real) {
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
        CidPath next() {
            final path = target.next()
            return path ? fromRealToCidPath(path, parentReal, parent) : null
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
        return toCidPath(path).getTargetPath(true).isHidden()
    }

    @Override
    FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException("File store not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        final cid = toCidPath(path)
        for( AccessMode m : modes ) {
            if( m == AccessMode.WRITE )
                throw new AccessDeniedException("Write mode not supported")
            if( m == AccessMode.EXECUTE )
                throw new AccessDeniedException("Execute mode not supported")
        }
        final real = cid.getTargetPath(true)
        if (real instanceof CidResultsPath)
            return
        real.fileSystem.provider().checkAccess(real, modes)
    }

    @Override
    <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        return null
    }

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        final cid = toCidPath(path)
        final real = cid.getTargetPath(true)
        if (real instanceof CidResultsPath)
            return (real as CidResultsPath).readAttributes(type)
        else
            return real.fileSystem.provider().readAttributes(real,type,options)
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
