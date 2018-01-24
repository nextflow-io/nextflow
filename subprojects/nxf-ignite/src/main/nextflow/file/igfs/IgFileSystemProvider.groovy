/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.file.igfs
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
import java.nio.file.attribute.BasicFileAttributeView
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import static nextflow.Const.ROLE_MASTER

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.daemon.IgGridFactory
import nextflow.util.OnlyCloseChannel
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteFileSystem
import org.apache.ignite.Ignition
/**
 * Implements an Ignite file system provider complaint with JSR203 a.k.a. NIO2
 *
 * @link https://jcp.org/en/jsr/detail?id=203
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgFileSystemProvider extends FileSystemProvider {

    /**
     * The 'scheme' string used by this provider i.e. the path url prefix that is
     * the string {@code ggfs}
     */
    public static final String SCHEME = 'igfs'

    private static final OpenOption[] CREATE_AND_TRUNCATE = [ StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.WRITE  ] as OpenOption[]

    private IgFileSystem currentFileSystem

    /**
     * The underlying {@link Ignite} instance
     */
    private Ignite grid

    /**
     * Returns the URI scheme that identifies this provider.
     *
     * @return  The string {@code 'igfs'}
     */
    @Override
    String getScheme() {
        SCHEME
    }

    /**
     * @return The underlying {@link Ignite} instance
     */
    Ignite getGrid() {
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
    synchronized IgFileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {

        if( currentFileSystem )
            throw new FileSystemAlreadyExistsException("File system object already exists -- URI: $uri")

        def ggfs = getIgniteFileSystem(env)
        return currentFileSystem = new IgFileSystem(this, ggfs)
    }

    /**
     * Short-cut for {@link #newFileSystem(java.net.URI, java.util.Map)}
     *
     * @param env
     * @return
     * @throws IOException
     */
    IgFileSystem newFileSystem(Map<String, ?> env) throws IOException {
        newFileSystem((URI)null, env)
    }


    static private Ignite getGridFor( Map env ) {

        def grid = env.get('grid')
        if( !grid ) {
            final session = env.session as Session
            if( !session )
                throw new IllegalStateException("Missing `session` object -- Cannot instantiate Ignite grid instance")
            final factory = new IgGridFactory(ROLE_MASTER, session.config ?: [:])
            return factory.start()
        }

        if( grid instanceof String )
            return Ignition.ignite(grid as String)

        if( grid instanceof UUID )
            return Ignition.ignite(grid as UUID)

        if( grid instanceof Ignite )
            return grid as Ignite

        throw new IllegalStateException("Cannot access underlying Apache Ignite instance -- not a valid grid identifier: [${grid?.class?.name} $grid]")
    }

    private IgniteFileSystem getIgniteFileSystem( Map env ) {

        // -- now find out the grid file system object
        def igfs = env.get('igfs')
        if( igfs instanceof IgniteFileSystem )
            return (IgniteFileSystem)igfs

        grid = getGridFor(env)
        if( !igfs )
            return grid.fileSystem('igfs')

        if( igfs instanceof String )
            return grid.fileSystem(igfs as String)


        throw new IllegalStateException("Cannot access underlying Apache Ignite file system instance -- not a valid fs idenitifier: [${igfs?.class?.name}] $igfs")
    }

    /**
     * @Inheritdoc
     *
     * @param uri
     * @return
     */
    @Override
    IgFileSystem getFileSystem(URI uri) {
        return currentFileSystem
    }

    /**
     * @Inheritdoc
     */
    @Override
    IgPath getPath(URI uri) {
        assert uri.getScheme() == SCHEME
        assert uri.path, 'igfs path cannot be empty'
        (IgPath)getFileSystem(uri).getPath(uri.path);
    }

    /**
     * Ignite do not support {@link SeekableByteChannel} implemented this method only
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

// TODO it should implement options validation
//        if( opts.size()!=2 || !opts.contains(StandardOpenOption.CREATE_NEW) || !opts.contains(StandardOpenOption.WRITE))
//            throw new UnsupportedOperationException("Options not supported: ${opts.join(', ')}")
//
        if( Files.exists(path) )
            throw new FileAlreadyExistsException("File igfs://${path} already exists")

        def gg = path as IgPath
        def stream = gg.fileSystem.igfs.create(gg.toIgnitePath(), false)
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

        def gg = path as IgPath
        gg.fileSystem.igfs.open( gg.toIgnitePath() )
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
        def gg = path as IgPath

        // set the default options
        if( !opts )
            opts = CREATE_AND_TRUNCATE

        if( StandardOpenOption.READ in opts )
            throw new IllegalArgumentException("READ not allowed");

        if( StandardOpenOption.CREATE_NEW in opts && Files.exists(path) )
            throw new FileAlreadyExistsException("File igfs://${path} already exists")


        final append = ( StandardOpenOption.APPEND in opts )
        if( append ) {
            boolean create = StandardOpenOption.CREATE in opts
            gg.fileSystem.igfs.append( gg.toIgnitePath(), create )
        }
        else {
            boolean overwrite = ( StandardOpenOption.TRUNCATE_EXISTING in opts )
            gg.fileSystem.igfs.create( gg.toIgnitePath(), overwrite )
        }

    }

    /**
     * @Inheritdoc
     */
    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        new IgDirectoryStream( dir as IgPath, filter )
    }

    /**
     * @Inheritdoc
     */
    @Override
    void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        (dir as IgPath).nativeMkdirs()
    }

    /**
     * @Inheritdoc
     */
    @Override
    void delete(Path path) throws IOException {
        (path as IgPath).nativeDelete()
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
                input.closeQuietly()
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
        def ggSource = source as IgPath
        def ggTarget = target as IgPath

        if( Files.exists(ggTarget) && !options.contains(StandardCopyOption.REPLACE_EXISTING))
            throw new FileAlreadyExistsException("File 'igfs://$target' already exists")

        ggSource.fileSystem.igfs.rename( ggSource.toIgnitePath(), ggTarget.toIgnitePath() )
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

    /**
     * @Inheritdoc
     *
     * @param path
     * @param modes
     * @throws IOException
     */
    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        def exists = (path as IgPath).nativeExists()
        if( !exists )
            throw new NoSuchFileException("File 'igfs://$path' do not exist")

        if( AccessMode.EXECUTE in modes )
            throw new UnsupportedOperationException("Execute access is not supported by Apache Ignite file system")
    }

    /**
     * @Inheritdoc
     *
     * @param obj
     * @param type
     * @param options
     * @return
     * @throws IOException
     */
    @Override
    def <V extends BasicFileAttributes> V readAttributes(Path obj, Class<V> type, LinkOption... options) throws IOException {
        def attrs = (obj as IgPath).nativeReadAttributes()
        if( attrs )
            return (V)new IgFileAttributes(attrs)
        else
            throw new NoSuchFileException("Cannot read attributes for file: $obj")
    }

    /**
     * @Inheritdoc
     *
     * @param path
     * @param type
     * @param options
     * @return
     */

    @Override
    def <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {

        if( BasicFileAttributeView.isAssignableFrom(type) && path instanceof IgPath )
            return (V)new IgFileAttributeView((IgPath)path)

        else
            throw new IllegalArgumentException("FileAttributeView of type: ${type?.name} is not supported")
    }


    /**
     * @Inheritdoc
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
     * @Inheritdoc
     *
     * @param path
     * @param attribute
     * @param value
     * @param options
     * @throws IOException
     */

    // TODO  setAttribute
    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException();
    }


}
