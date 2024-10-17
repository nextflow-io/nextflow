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
 */

package nextflow.file.http

import nextflow.file.CopyMoveHelper

import static nextflow.file.http.XFileSystemConfig.*

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
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.util.InsensitiveMap
import sun.net.www.protocol.ftp.FtpURLConnection
/**
 * Implements a read-only JSR-203 compliant file system provider for http/ftp protocols
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 */
@Slf4j
@PackageScope
@CompileStatic
abstract class XFileSystemProvider extends FileSystemProvider {

    private Map<URI, FileSystem> fileSystemMap = new LinkedHashMap<>(20)

    private static final int[] REDIRECT_CODES = [301, 302, 307, 308]

    protected static String config(String name, def defValue) {
        return SysEnv.containsKey(name) ? SysEnv.get(name) : defValue.toString()
    }

    static private URI key(String s, String a) {
        new URI("$s://$a")
    }

    static private URI key(URI uri) {
        final base = uri.authority
        int p = base.indexOf('@')
        if( p==-1 )
            return key(uri.scheme.toLowerCase(), base.toLowerCase())
        else {
            final user = base.substring(0,p)
            final host = base.substring(p)
            return key(uri.scheme.toLowerCase(), user + host.toLowerCase())
        }
    }

    @Override
    FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        final scheme = uri.scheme.toLowerCase()

        if( scheme != this.getScheme() )
            throw new IllegalArgumentException("Not a valid ${getScheme().toUpperCase()} scheme: $scheme")

        final base = key(uri)
        if (fileSystemMap.containsKey(base))
            throw new IllegalStateException("File system `$base` already exists")

        return new XFileSystem(this, base)
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

        final scheme = uri.scheme.toLowerCase()

        if( scheme != this.getScheme() )
            throw new IllegalArgumentException("Not a valid ${getScheme().toUpperCase()} scheme: $scheme")

        final key = key(uri)

        if( !canCreate ) {
            FileSystem result = fileSystemMap[key]
            if( result==null )
                throw new FileSystemNotFoundException("File system not found: $key")
            return result
        }

        synchronized (fileSystemMap) {
            FileSystem result = fileSystemMap[key]
            if( result==null ) {
                result = newFileSystem(uri,Collections.emptyMap())
                fileSystemMap[key] = result
            }
            return result
        }
    }

    @Override
    Path getPath(URI uri) {
        def path = uri.path
        if( !path.contains('?') && uri.query )
            path += '?' + uri.query
        return getFileSystem(uri,true).getPath(path)
    }

    protected String auth(String userInfo) {
        final String BEARER = 'x-oauth-bearer:'
        int p = userInfo.indexOf(BEARER)
        if( p!=-1 ) {
            final token = userInfo.substring(BEARER.length())
            return "Bearer $token"
        }
        else {
            return "Basic ${userInfo.getBytes().encodeBase64()}"
        }
    }

    protected URLConnection toConnection(Path path) {
        final url = path.toUri().toURL()
        log.trace "File remote URL: $url"
        return toConnection0(url, 0)
    }

    protected URLConnection toConnection0(URL url, int attempt) {
        final conn = url.openConnection()
        conn.setRequestProperty("User-Agent", 'Nextflow/httpfs')
        if( conn instanceof HttpURLConnection ) {
            // by default HttpURLConnection does redirect only within the same host
            // disable the built-in to implement custom redirection logic (see below)
            conn.setInstanceFollowRedirects(false)
        }
        if( url.userInfo ) {
            conn.setRequestProperty("Authorization", auth(url.userInfo));
        }
        else {
            XAuthRegistry.instance.authorize(conn)
        }
        if ( conn instanceof HttpURLConnection && conn.getResponseCode() in REDIRECT_CODES && attempt < MAX_REDIRECT_HOPS ) {
            final header = InsensitiveMap.of(conn.getHeaderFields())
            final location = header.get("Location")?.get(0)
            log.debug "Remote redirect location: $location"
            final newUrl = new URI(absLocation(location,url)).toURL()
            if( url.protocol=='https' && newUrl.protocol=='http' )
                throw new IOException("Refuse to follow redirection from HTTPS to HTTP (unsafe) URL - origin: $url - target: $newUrl")
            return toConnection0(newUrl, attempt+1)
        }
        else if( conn instanceof HttpURLConnection && conn.getResponseCode() in config().retryCodes() && attempt < config().maxAttempts() ) {
            final delay = (Math.pow(config().backOffBase(), attempt) as long) * config().backOffDelay()
            log.debug "Got HTTP error=${conn.getResponseCode()} waiting for ${delay}ms (attempt=${attempt+1})"
            Thread.sleep(delay)
            return toConnection0(url, attempt+1)
        }
        else if( conn instanceof HttpURLConnection && conn.getResponseCode()==401 && attempt==0 ) {
            if( XAuthRegistry.instance.refreshToken(conn) ) {
                return toConnection0(url, attempt+1)
            }
        }
        return conn
    }

    protected String absLocation(String location, URL target) {
        assert location, "Missing location argument"
        assert target, "Missing target URL argument"

        final base = FileHelper.baseUrl(location)
        if( base )
            return location
        if( !location.startsWith('/') )
            location = '/' + location
        return "${target.protocol}://${target.authority}$location"
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

        final conn = toConnection(path)
        final length = conn.getContentLengthLong()
        final target = length>0
                ? new FixedInputStream(conn.getInputStream(),length)
                : conn.getInputStream()
        final stream = new BufferedInputStream(target)

        new SeekableByteChannel() {

            private long _position

            @Override
            int read(ByteBuffer buffer) throws IOException {
                def data=0
                int len=0
                while( buffer.hasRemaining() && (data=stream.read())!=-1 ) {
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
                // this value is going to be used as the buffer size
                // file related operation. See for example {@link Files#readAllBytes}
                return 8192
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
    InputStream newInputStream(Path path, OpenOption... options)
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

        final conn = toConnection(path)
        final length = conn.getContentLengthLong()
        // only apply the FixedInputStream check if staging files
        return length>0 && CopyMoveHelper.IN_FOREIGN_COPY.get()
            ? new FixedInputStream(conn.getInputStream(), length)
            : conn.getInputStream()
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
        throw new UnsupportedOperationException("Directory listing unsupported by ${getScheme().toUpperCase()} file system provider")
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
                throw new IOException("Unable to access path: ${FilesEx.toUriString(p)}")
            }
            return attrs
        }
        throw new UnsupportedOperationException("Not a valid ${getScheme().toUpperCase()} file attribute type: $type")
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Read file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Set file attributes not supported by ${getScheme().toUpperCase()} file system provider")
    }

    protected XFileAttributes readHttpAttributes(XPath path) {
        final conn = toConnection(path)
        if( conn instanceof FtpURLConnection ) {
            return new XFileAttributes(null,-1)
        }
        if ( conn instanceof HttpURLConnection && conn.getResponseCode() in [200, 301, 302, 307, 308]) {
            final header = conn.getHeaderFields()
            return readHttpAttributes(header)
        }
        return null
    }

    protected XFileAttributes readHttpAttributes(Map<String,List<String>> header) {
        final header0 = InsensitiveMap.<String,List<String>>of(header)
        def lastMod = header0.get("Last-Modified")?.get(0)
        long contentLen = header0.get("Content-Length")?.get(0)?.toLong() ?: -1L
        def dateFormat = new SimpleDateFormat('E, dd MMM yyyy HH:mm:ss Z', Locale.ENGLISH) // <-- make sure date parse is not language dependent (for the week day)
        def modTime = lastMod ? FileTime.from(dateFormat.parse(lastMod).time, TimeUnit.MILLISECONDS) : (FileTime)null
        new XFileAttributes(modTime, contentLen)
    }

}
