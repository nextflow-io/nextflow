/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.file.ggfs

import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileStore
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.GgGridFactory
import nextflow.extension.FilesExtensions
import nextflow.util.OnlyCloseChannel
import org.gridgain.grid.Grid
import org.gridgain.grid.GridGain
import org.gridgain.grid.ggfs.GridGgfs
/**
 * Implements a GridGain file system provider complaint with JSR203 a.k.a. NIO2
 *
 * @link https://jcp.org/en/jsr/detail?id=203
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GgFileSystemProvider extends FileSystemProvider {

    /**
     * The 'scheme' string used by this provider i.e. the path url prefix that is
     * the string {@code ggfs}
     */
    public static final String SCHEME = 'ggfs'

    private static final OpenOption[] CREATE_AND_TRUNCATE = [ StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.WRITE  ] as OpenOption[]

    private GgFileSystem currentFileSystem

    /**
     * The underlying {@link Grid} instance
     */
    private Grid grid

    /**
     * Returns the URI scheme that identifies this provider.
     *
     * @return  The string {@code 'ggfs'}
     */
    @Override
    String getScheme() {
        SCHEME
    }

    /**
     * @return The underlying {@link Grid} instance
     */
    Grid getGrid() {
        grid
    }

    /**
     * Create the {@code FileSystem} object for the provided URI.
     * Actually a single file system instance independently of the specified URI.
     * <p>
     * The underlying grid and file systems object have to be specified using the {@code env} map
     * object.
     * <p>
     * For example:
     * <code>
     *      def gridObj = GridGain.grid().ggfs('gg_fs_name')
     *      new GgFileSystemProvider().newFileSystem(null, [ ggfs: gridObj ] )
     * </code>
     *
     * @param uri This parameter is ignored since only a single file system instance is supported
     * @param env Provides the {@code Grid} and {@code GridGgfs} object instances. The following parameters are supported:
     *              <li>{@code ggfs} The {@code GridGgfs} object instance - or - the name of the grid FS instance, retried
     *              by using the method {@code Grid#ggfs(String)}. Default {@code 'ggfs'}
     *              <li>{@code grid} The {@code Grid} object instance - or -the name of the grid object retried by using the
     *              method {@code GridGain.grid(String)}. Note: it is required only when the {@code ggfs} parameter is specified
     *              by the name.
     *
     * @return  The newly created {@code GgFileSystem} instance
     * @throws  FileSystemAlreadyExistsException
     *          If the file system has already been created
     */
    @Override
    synchronized GgFileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {

        if( currentFileSystem )
            throw new FileSystemAlreadyExistsException("File system object already exists -- URI: $uri")

        def ggfs = getGgfsFor(env)
        return currentFileSystem = new GgFileSystem(this, ggfs)
    }

    /**
     * Short-cut for {@link #newFileSystem(java.net.URI, java.util.Map)}
     *
     * @param env
     * @return
     * @throws IOException
     */
    GgFileSystem newFileSystem(Map<String, ?> env) throws IOException {
        newFileSystem((URI)null, env)
    }


    static private Grid getGridFor( Map env ) {

        def grid = env.get('grid')
        if( !grid ) {
            if( env.containsKey('session') ) {
                return new GgGridFactory( env.session as Session ).start()
            }
            else {
                return GridGain.grid()
            }
        }

        if( grid instanceof String )
            return GridGain.grid(grid as String)

        if( grid instanceof UUID )
            return GridGain.grid(grid as UUID)

        if( grid instanceof Grid )
            return grid as Grid

        throw new IllegalStateException("Cannot access underlying GridGain instance -- not a valid grid identifier: [${grid?.class?.name} $grid]")
    }

    private GridGgfs getGgfsFor( Map env ) {

        // -- now find out the grid file system object
        def ggfs = env.get('ggfs')
        if( ggfs instanceof GridGgfs )
            return (GridGgfs)ggfs

        grid = getGridFor(env)
        if( !ggfs )
            return grid.ggfs('ggfs')

        if( ggfs instanceof String )
            return grid.ggfs(ggfs as String)


        throw new IllegalStateException("Cannot access underlying GridGain file system instance -- not a valid fs idenitifier: [${ggfs?.class?.name}] $ggfs")
    }

    @Override
    GgFileSystem getFileSystem(URI uri) {
        return currentFileSystem
    }

    /**
     * @Inheritdoc
     */
    @Override
    GgPath getPath(URI uri) {
        assert uri.getScheme() == SCHEME
        (GgPath)getFileSystem(uri).getPath(uri.getPath());
    }

    /**
     * GridGain do not support {@code SeekableByteChannel} implemented this method only
     * to have {@code Files.create} to work properly
     *
     * @param path
     * @param opts
     * @param attrs
     * @return
     * @throws IOException
     */
    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> opts, FileAttribute<?>... attrs) throws IOException {

        if( opts.size()!=2 || !opts.contains(StandardOpenOption.CREATE_NEW) || !opts.contains(StandardOpenOption.WRITE))
            throw new UnsupportedOperationException("Options not supported: ${opts.join(', ')}")

        if( Files.exists(path) )
            throw new FileAlreadyExistsException("File ggfs://${path} already exists")

        def gg = path as GgPath
        def stream = gg.fileSystem.ggfs.create(gg.toGridGgfsPath(), false)
        new OnlyCloseChannel(stream)
    }

    /**
     * Opens a file, returning an input stream to read from the file. This
     * method works in exactly the manner specified by the {@link
     * java.nio.file.Files#newInputStream} method.
     *
     * <p> The default implementation of this method opens a channel to the file
     * as if by invoking the {@link #newByteChannel} method and constructs a
     * stream that reads bytes from the channel. This method should be overridden
     * where appropriate.
     *
     * @param   path
     *          the path to the file to open
     * @param   opts
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
    public InputStream newInputStream(Path path, OpenOption... opts)
            throws IOException
    {
        for (OpenOption opt: opts) {
            if (opt != StandardOpenOption.READ)
                throw new UnsupportedOperationException("'" + opt + "' not allowed");
        }

        def gg = path as GgPath
        gg.fileSystem.ggfs.open( gg.toGridGgfsPath() )
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
    public OutputStream newOutputStream(Path path, OpenOption... opts)
            throws IOException
    {
        def gg = path as GgPath

        // set the default options
        if( !opts )
            opts = CREATE_AND_TRUNCATE

        if( StandardOpenOption.READ in opts )
            throw new IllegalArgumentException("READ not allowed");

        if( StandardOpenOption.CREATE_NEW in opts && Files.exists(path) )
            throw new FileAlreadyExistsException("File ggfs://${path} already exists")


        final append = ( StandardOpenOption.APPEND in opts )
        if( append ) {
            boolean create = StandardOpenOption.CREATE in opts
            gg.fileSystem.ggfs.append( gg.toGridGgfsPath(), create )
        }
        else {
            boolean overwrite = ( StandardOpenOption.TRUNCATE_EXISTING in opts )
            gg.fileSystem.ggfs.create( gg.toGridGgfsPath(), overwrite )
        }

    }

    /**
     * @Inheritdoc
     */
    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        new GgDirectoryStream( dir as GgPath, filter )
    }

    /**
     * @Inheritdoc
     */
    @Override
    void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        (dir as GgPath).nativeMkdirs()
    }

    /**
     * @Inheritdoc
     */
    @Override
    void delete(Path path) throws IOException {
        (path as GgPath).nativeDelete()
    }

    /**
     * Copy a file to a target file. This method works in exactly the manner
     * specified by the {@link Files#copy(Path,Path,CopyOption[])} method
     * except that both the source and target paths must be associated with
     * this provider.
     *
     * @param   source
     *          the path to the file to copy
     * @param   target
     *          the path to the target file
     * @param   options
     *          options specifying how the copy should be done
     *
     * @throws  UnsupportedOperationException
     *          if the array contains a copy option that is not supported
     * @throws  FileAlreadyExistsException
     *          if the target file exists but cannot be replaced because the
     *          {@code REPLACE_EXISTING} option is not specified <i>(optional
     *          specific exception)</i>
     * @throws  java.nio.file.DirectoryNotEmptyException
     *          the {@code REPLACE_EXISTING} option is specified but the file
     *          cannot be replaced because it is a non-empty directory
     *          <i>(optional specific exception)</i>
     * @throws  IOException
     *          if an I/O error occurs
     * @throws  SecurityException
     *          In the case of the default provider, and a security manager is
     *          installed, the {@link SecurityManager#checkRead(String) checkRead}
     *          method is invoked to check read access to the source file, the
     *          {@link SecurityManager#checkWrite(String) checkWrite} is invoked
     *          to check write access to the target file. If a symbolic link is
     *          copied the security manager is invoked to check {@link
     *          java.nio.file.LinkPermission}{@code ("symbolic")}.
     */
    @Override
    void copy(Path source, Path target, CopyOption... opts) throws IOException {

        // delete target if it exists and REPLACE_EXISTING is specified
        if ( StandardCopyOption.REPLACE_EXISTING in opts ) {
            Files.deleteIfExists(target)
        }
        else if( Files.exists(target) )
            throw new FileAlreadyExistsException(target.toString());

        // create directory or copy file
        if (Files.isDirectory(source)) {
            Files.createDirectory(target);
        }
        else {
            InputStream input = Files.newInputStream(source)
            try {
                Files.copy(input, target);
            }
            finally {
                FilesExtensions.closeQuietly(input)
            }
        }

    }

    /**
     * Move or rename a file to a target file. This method works in exactly the
     * manner specified by the {@link Files#move} method except that both the
     * source and target paths must be associated with this provider.
     *
     * @param   source
     *          the path to the file to move
     * @param   target
     *          the path to the target file
     * @param   options
     *          options specifying how the move should be done
     *
     * @throws  UnsupportedOperationException
     *          if the array contains a copy option that is not supported
     * @throws  FileAlreadyExistsException
     *          if the target file exists but cannot be replaced because the
     *          {@code REPLACE_EXISTING} option is not specified <i>(optional
     *          specific exception)</i>
     * @throws  java.nio.file.DirectoryNotEmptyException
     *          the {@code REPLACE_EXISTING} option is specified but the file
     *          cannot be replaced because it is a non-empty directory
     *          <i>(optional specific exception)</i>
     * @throws  java.nio.file.AtomicMoveNotSupportedException
     *          if the options array contains the {@code ATOMIC_MOVE} option but
     *          the file cannot be moved as an atomic file system operation.
     * @throws  IOException
     *          if an I/O error occurs
     * @throws  SecurityException
     *          In the case of the default provider, and a security manager is
     *          installed, the {@link SecurityManager#checkWrite(String) checkWrite}
     *          method is invoked to check write access to both the source and
     *          target file.
     */
    @Override
    void move(Path source, Path target, CopyOption... options) throws IOException {
        def ggSource = source as GgPath
        def ggTarget = target as GgPath

        if( Files.exists(ggTarget) && !options.contains(StandardCopyOption.REPLACE_EXISTING))
            throw new FileAlreadyExistsException("File 'ggfs://$target' already exists")

        ggSource.fileSystem.ggfs.rename( ggSource.toGridGgfsPath(), ggTarget.toGridGgfsPath() )
    }

    /**
     * Tests if two paths locate the same file. This method works in exactly the
     * manner specified by the {@link Files#isSameFile} method.
     *
     * @param   path
     *          one path to the file
     * @param   path2
     *          the other path
     *
     * @return  {@code true} if, and only if, the two paths locate the same file
     *
     * @throws  IOException
     *          if an I/O error occurs
     * @throws  SecurityException
     *          In the case of the default provider, and a security manager is
     *          installed, the {@link SecurityManager#checkRead(String) checkRead}
     *          method is invoked to check read access to both files.
     */
    @Override
    boolean isSameFile(Path path1, Path path2) throws IOException {
        return path1.toRealPath() == path2.toRealPath()
    }

    /**
     * Not supported
     *
     * @param path A path object
     * @return {@code false}
     */
    @Override
    boolean isHidden(Path path) throws IOException {
        return false
    }

    /**
     * Not supported
     *
     * @param path A path object
     * @return {@code null}
     */
    @Override
    FileStore getFileStore(Path path) throws IOException {
        return null
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        def exists = (path as GgPath).nativeExists()
        if( !exists )
            throw new NoSuchFileException("File 'ggfs://$path' do not exist")

        if( AccessMode.EXECUTE in modes )
            throw new UnsupportedOperationException("Execute access is not supported by GridGain file system")
    }

    @Override
    def <A extends BasicFileAttributes> A readAttributes(Path obj, Class<A> type, LinkOption... options) throws IOException {
        new GgFileAttributes( (obj as GgPath).nativeReadAttributes() )
    }

    /**
     * Not implemented
     *
     * @param path
     * @param type
     * @param options
     * @return
     */

    // TODO getFileAttributeView
    @Override
    def <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        throw new UnsupportedOperationException("Method 'getFileAttributeView' not supported");
    }


    /**
     * Not implemented
     *
     * @param path
     * @param attributes
     * @param options
     * @return
     * @throws IOException
     */
    // TODO readAttributes
    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException();
    }

    /**
     * Not implemented
     * @param path
     * @param attribute
     * @param value
     * @param options
     * @throws IOException
     */

    // TODO
    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException();
    }


}
