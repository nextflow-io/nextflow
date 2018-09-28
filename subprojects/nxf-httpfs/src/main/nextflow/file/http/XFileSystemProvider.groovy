/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.file.http

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.attribute.FileTime
import java.nio.file.spi.FileSystemProvider
import java.text.SimpleDateFormat
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import sun.net.www.protocol.ftp.FtpURLConnection

/**
 * Implements a read-only JSR-203 complaint file system provider for http/ftp protocols
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 */
@PackageScope
@CompileStatic
abstract class XFileSystemProvider extends FileSystemProvider {

    private Map<URI, FileSystem> fileSystemMap = [:]

    static public Set<String> ALL_SCHEMES = ['ftp','http','https'] as Set

    static private URI key(String s, String a) {
        new URI("$s://$a")
    }

    static private URI key(URI uri) {
        key(uri.scheme.toLowerCase(), uri.authority.toLowerCase())
    }

    @Override
    FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        final scheme = uri.scheme.toLowerCase()

        if( scheme != this.getScheme() )
            throw new IllegalArgumentException("Not a valid ${getScheme().toUpperCase()} scheme: $scheme")

        final base = key(uri)
        if (fileSystemMap.containsKey(base))
            throw new IllegalStateException("File system `$base` already exists")

        def result = new XFileSystem(this, base)
        fileSystemMap[base] = result
        return result
    }

    /**
     * Returns an existing {@code FileSystem} created by this provider.
     *
     * <p> This method returns a reference to a {@code FileSystem} that was
     * created by invoking the {@link #newFileSystem(URI,Map) newFileSystem(URI,Map)}
     * method. File systems created the {@link #newFileSystem(Path,Map)
     * newFileSystem(Path,Map)} method are not returned by this method.
     * The file system is identified by its {@code URI}. Its exact form
     * is highly provider dependent. In the case of the default provider the URI's
     * path component is {@code "/"} and the authority, query and fragment components
     * are undefined (Undefined components are represented by {@code null}).
     *
     * @param   uri
     *          URI reference
     *
     * @return  The file system
     *
     * @throws  IllegalArgumentException
     *          If the pre-conditions for the {@code uri} parameter aren't met
     * @throws  FileSystemNotFoundException
     *          If the file system does not exist
     * @throws  SecurityException
     *          If a security manager is installed and it denies an unspecified
     *          permission.
     */
    @Override
    FileSystem getFileSystem(URI uri) {
        getFileSystem(uri,false)
    }

    FileSystem getFileSystem(URI uri, boolean canCreate) {
        assert fileSystemMap != null

        def scheme = uri.scheme.toLowerCase()

        if( scheme != this.getScheme() )
            throw new IllegalArgumentException("Not a valid ${getScheme().toUpperCase()} scheme: $scheme")

        def key = key(uri)

        FileSystem result = fileSystemMap[key]
        if( !result ) {
            if( canCreate )
                result = newFileSystem(uri,Collections.emptyMap())
            else
                throw new FileSystemNotFoundException("File system not found: $key")
        }

        return result
    }

    @Override
    Path getPath(URI uri) {
        return getFileSystem(uri,true).getPath(uri.path)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {

        if (path.class != XPath)
            throw new ProviderMismatchException()

        if (options.size() > 0) {
            for (OpenOption opt: options) {
                // All OpenOption values except for APPEND and WRITE are allowed
                if (opt == StandardOpenOption.APPEND || opt == StandardOpenOption.WRITE)
                    throw new UnsupportedOperationException("'$opt' not allowed");
            }
        }

        final conn = new URL(path.toUri().toString()).openConnection()
        final size = conn.getContentLengthLong()
        final stream = new BufferedInputStream(conn.getInputStream())

        new SeekableByteChannel() {

            private long _position

            @Override
            int read(ByteBuffer buffer) throws IOException {
                def data=0
                int len=0
                while( len<buffer.capacity() && (data=stream.read())!=-1 ) {
                    buffer.put((byte)data)
                    len++
                }
                _position += len
                return len ?: -1
            }

            @Override
            int write(ByteBuffer src) throws IOException {
                throw new UnsupportedOperationException("Write operation not supported")
            }

            @Override
            long position() throws IOException {
                return _position
            }

            @Override
            SeekableByteChannel position(long newPosition) throws IOException {
                throw new UnsupportedOperationException("Position operation not supported")
            }

            @Override
            long size() throws IOException {
                return size
            }

            @Override
            SeekableByteChannel truncate(long unused) throws IOException {
                throw new UnsupportedOperationException("Truncate operation not supported")
            }

            @Override
            boolean isOpen() {
                return true
            }

            @Override
            void close() throws IOException {
                stream.close()
            }
        }
    }

    /**
     * Opens a file, returning an input stream to read from the file. This
     * method works in exactly the manner specified by the {@link
     * Files#newInputStream} method.
     *
     * <p> The default implementation of this method opens a channel to the file
     * as if by invoking the {@link #newByteChannel} method and constructs a
     * stream that reads bytes from the channel. This method should be overridden
     * where appropriate.
     *
     * @param   path
     *          the path to the file to open
     * @param   options
     *          options specifying how the file is opened
     *
     * @return  a new input stream
     *
     * @throws  IllegalArgumentException
     *          if an invalid combination of options is specified
     * @throws  UnsupportedOperationException
     *          if an unsupported option is specified
     * @throws  IOException
     *          if an I/O error occurs
     * @throws  SecurityException
     *          In the case of the default provider, and a security manager is
     *          installed, the {@link SecurityManager#checkRead(String) checkRead}
     *          method is invoked to check read access to the file.
     */
    @Override
    public InputStream newInputStream(Path path, OpenOption... options)
            throws IOException
    {
        if (path.class != XPath)
            throw new ProviderMismatchException()

        if (options.length > 0) {
            for (OpenOption opt: options) {
                // All OpenOption values except for APPEND and WRITE are allowed
                if (opt == StandardOpenOption.APPEND ||
                        opt == StandardOpenOption.WRITE)
                    throw new UnsupportedOperationException("'$opt' not allowed");
            }
        }

        return new URL(path.toUri().toString()).openStream()
    }

    /**
     * Opens or creates a file, returning an output stream that may be used to
     * write bytes to the file. This method works in exactly the manner
     * specified by the {@link Files#newOutputStream} method.
     *
     * <p> The default implementation of this method opens a channel to the file
     * as if by invoking the {@link #newByteChannel} method and constructs a
     * stream that writes bytes to the channel. This method should be overridden
     * where appropriate.
     *
     * @param   path
     *          the path to the file to open or create
     * @param   options
     *          options specifying how the file is opened
     *
     * @return  a new output stream
     *
     * @throws  IllegalArgumentException
     *          if {@code options} contains an invalid combination of options
     * @throws  UnsupportedOperationException
     *          if an unsupported option is specified
     * @throws  IOException
     *          if an I/O error occurs
     * @throws  SecurityException
     *          In the case of the default provider, and a security manager is
     *          installed, the {@link SecurityManager#checkWrite(String) checkWrite}
     *          method is invoked to check write access to the file. The {@link
     *          SecurityManager#checkDelete(String) checkDelete} method is
     *          invoked to check delete access if the file is opened with the
     *          {@code DELETE_ON_CLOSE} option.
     */
    @Override
    public OutputStream newOutputStream(Path path, OpenOption... options)
            throws IOException
    {
        throw new UnsupportedOperationException("Write not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        throw new UnsupportedOperationException("Direcotry listing unsupported by ${getScheme().toUpperCase()} file system provider")
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
        return path.getFileName().startsWith('.')
    }

    @Override
    FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException("File store not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        readAttributes(path, XFileAttributes)

        for( AccessMode m : modes ) {
            if( m == AccessMode.WRITE )
                throw new AccessDeniedException("Write mode not supported")
            if( m == AccessMode.EXECUTE )
                throw new AccessDeniedException("Execute mode not supported")
        }
    }

    @Override
    def <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        return null
    }

    @Override
    def <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        if ( type == BasicFileAttributes || type == XFileAttributes) {
            def p = (XPath) path
            def attrs = (A)readHttpAttributes(p)
            if (attrs == null) {
                throw new IOException("Unable to access path: ${p.toString()}")
            }
            return attrs
        }
        throw new UnsupportedOperationException("Not a valid ${getScheme().toUpperCase()} file attribute type: $type")
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Read file attributes no supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Set file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

    protected XFileAttributes readHttpAttributes(XPath path) {
        final conn = path.toUri().toURL().openConnection()
        if( conn instanceof FtpURLConnection ) {
            return new XFileAttributes(null,-1)
        }
        if ( conn instanceof HttpURLConnection && conn.getResponseCode() in [200, 301, 302]) {
            def header = conn.getHeaderFields()
            return readHttpAttributes(header)
        }
        return null
    }

    protected XFileAttributes readHttpAttributes(Map<String,List<String>> header) {
        def lastMod = header.get("Last-Modified")?.get(0)
        long contentLen = header.get("Content-Length")?.get(0)?.toLong() ?: -1
        def dateFormat = new SimpleDateFormat('E, dd MMM yyyy HH:mm:ss Z')
        def modTime = lastMod ? FileTime.from(dateFormat.parse(lastMod).time, TimeUnit.MILLISECONDS) : (FileTime)null
        new XFileAttributes(modTime, contentLen)
    }

}
