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

package nextflow.extension
import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.NotLinkException
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.StandardCopyOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileTime
import java.nio.file.attribute.PosixFilePermission
import java.nio.file.attribute.PosixFilePermissions

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.util.ByteBufferBackedInputStream
import nextflow.util.CharsetHelper
import nextflow.util.FileHelper
/**
 * Add utility methods to standard classes {@code File} and {@code Path}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class FilesExtensions {

    private static CR = 0x0D

    private static LF = 0x0A


    /**
     * Check if a file - or - a directory is empty
     *
     * @param file The file under test
     * @return {@code true} if the file does not exist or it is empty
     */
    def static boolean empty( File file ) {
        FileHelper.empty(file)
    }


    /**
     *
     * Check if a file - or - a directory is empty
     *
     * @param file The file under test
     * @return {@code true} if the file does not exist or it is empty
     */
    def static boolean empty( Path path ) {
        FileHelper.empty(path)
    }

    /**
     * Deletes the file or directory denoted by this abstract pathname.  If
     * this pathname denotes a directory, then the directory must be empty in
     * order to be deleted.
     *
     * @see File#delete()
     */

    def static boolean delete( Path path ) {
        try {
            Files.delete(path)
            return true
        }
        catch( IOException e ) {
            return false
        }
    }

    /**
     * Deletes a directory with all contained files and subdirectories.
     * <p>The method returns
     * <ul>
     * <li>true, when deletion was successful</li>
     * <li>true, when it is called for a non existing directory</li>
     * <li>false, when it is called for a file which isn't a directory</li>
     * <li>false, when directory couldn't be deleted</li>
     * </ul>
     *
     * @param path
     * @return
     */
    def static boolean deleteDir(Path path) {
        if( !Files.exists(path) )
             return true

        if( !Files.isDirectory(path) )
            return false

        try {
            Files.walkFileTree(path, new SimpleFileVisitor<Path>() {

                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                    Files.delete(file)
                    FileVisitResult.CONTINUE
                }

                public FileVisitResult postVisitDirectory(Path dir, IOException exc) {
                    Files.delete(dir)
                    FileVisitResult.CONTINUE
                }

            })
            return true
        }
        catch( IOException e ) {
            return false
        }
    }

    /**
     * Copy or a file or a directory. It mimics the semantic of the Linux *cp* command
     *
     * @param source
     * @param target
     * @return
     */
    def static File copyTo( File source, File target ) {
        copyTo(source.toPath(), target.toPath()).toFile()
    }

    /**
     * Copy a file or a directory to the target path.
     * It mimics the semantic of Linux *cp* command.
     *
     * @param source
     * @param target
     * @return
     */
    def static Path copyTo( Path source, Path target ) {

        if( source.isDirectory() ) {
            def parent = target.getParent()
            if( parent && !parent.exists() ) parent.mkdirs()
            return copyDirectory(source,target)
        }

        if( target.isDirectory() ) {
            target = target.resolve(source.getName())
            return Files.copy(source,target, StandardCopyOption.REPLACE_EXISTING)
        }

        // create the parent directories if do not exist
        def parent = source.getParent()
        if( parent && !parent.exists() ) {
            parent.mkdirs()
        }

        return Files.copy(source,target, StandardCopyOption.REPLACE_EXISTING)
    }


    private static Path copyDirectory( Path source, Path target ) {

        def visitor = new SimpleFileVisitor<Path>() {

            public FileVisitResult preVisitDirectory(Path current, BasicFileAttributes attr)
            throws IOException
            {
                // get the *delta* path against the source path
                def delta = source.relativize(current)
                def newFolder = delta ? target.resolve(delta) : target
                FilesExtensions.log.trace "Copy DIR: $current -> $newFolder"
                if( !newFolder.exists() ) {
                    Files.copy(current, newFolder, StandardCopyOption.REPLACE_EXISTING)
                }
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFile(Path current, BasicFileAttributes attr)
            throws IOException
            {
                // get the *delta* path against the source path
                def delta = source.relativize(current)
                def newFile = delta ? target.resolve(delta) : target
                FilesExtensions.log.trace "Copy file: $current -> $newFile"
                Files.copy(current, newFile, StandardCopyOption.REPLACE_EXISTING)
                return FileVisitResult.CONTINUE;
            }

        }

        Files.walkFileTree(source, visitor)
        return target
    }

    /**
     * Copy a file  -or- a whole directory to a new file -or- a new directory
     *
     * @param source
     * @param target
     * @return
     */
    def static File copyTo( File source, String target ) {
        copyTo(source.toPath(), FileHelper.asPath(target)).toFile()
    }

    /**
     * Copy or a file or a directory. It mimics the semantic of the Linux *cp* command
     *
     * @param source
     * @param target
     * @return
     */
    def static Path copyTo( Path source, String target ) {
        copyTo(source, FileHelper.asPath(target))
    }

    /**
     * Move a file or a directory. Mimics the Linux *mv* command
     *
     * @param source
     * @param target
     * @return
     */
    def static File moveTo( File source, File target ) {
        moveTo(source.toPath(), target.toPath()).toFile()
    }

    /**
     * Move a file or a directory. Mimics the Linux *mv* command
     *
     * @param source
     * @param target
     * @return
     */
    def static File moveTo( File source, String target ) {
        moveTo(source.toPath(), FileHelper.asPath(target)).toFile()
    }

    /**
     * Move a file or a directory. Mimics the Linux *mv* command
     *
     * @param source
     * @param target
     * @return
     */
    def static Path moveTo( Path source, Path target ) {

        if( source.isDirectory() ) {
            if( target.isDirectory() ) {
                // when the target path is a directory
                // move the source path into the target folder
                return Files.move(source, target.resolve(source.getFileName()))
            }
            else {
                def parent = target.getParent()
                if( parent && !parent.exists()) { parent.mkdirs() }
                return Files.move(source, target)
            }
        }

        // when the target path is a directory
        // move the source file into the target using the same name
        if( target.isDirectory() ) {
            target = target.resolve(source.getName())
            return Files.move(source, target)
        }

        // create the parent directories if do not exist
        def parent = target.getParent()
        if( parent && !parent.exists() ) { parent.mkdirs() }

        // move the source file to the target file,
        // overwriting the target if exists
        return Files.move(source, target, StandardCopyOption.REPLACE_EXISTING)
    }

    def static Path moveTo( Path source, String target ) {
        moveTo(source, FileHelper.asPath(target) )
    }

    /**
     * Gets the base name, minus the full path and extension, from a full filename.
     *
     * This method will handle a file in either Unix or Windows format.
     * The text after the last forward or backslash and before the last dot is returned.
     * <pre>
     *   a/b/c.txt --> c
     *   a.txt     --> a
     *   a/b/c     --> c
     *   a/b/c/    --> ""
     * </pre>
     *
     * The output will be the same irrespective of the machine that the code is running on.
     * @param file The filename to query, null returns null
     * @return The name of the file without the path, or an empty string if none exists
     */
    def static String getBaseName( File file ) {
        getBaseName(file.toPath())
    }

    def static String getBaseName( Path self ) {
        assert self

        String name = self.getFileName()
        if( !name ) return ''

        int pos = name.lastIndexOf('.')
        if( pos == -1 ) return name.toString()

        return name.substring(0, pos)
    }

    /**
     * Extend {@code Path} adding a getName method semantically equivalent to {@code File#getName}
     *
     * @return The name of the path as string
     */
    def static String getName( Path self ) {
       return self.getFileName()?.toString() ?: ''
    }

    /**
     * Gets the extension of a filename.
     * This method returns the textual part of the filename after the last dot.
     * There must be no directory separator after the dot.
     * <pre>
     *   foo.txt      --> "txt"
     *   a/b/c.jpg    --> "jpg"
     *   a/b.txt/c    --> ""
     *   a/b/c        --> ""
     * </pre>
     *
     * The output will be the same irrespective of the machine that the code is running on.
     *
     *
     * @param file  The file to retrieve the extension of.
     * @return the Extension of the file or an empty string if none exists
     */
    def static String getExtension( File file ) {
        getExtension(file.toPath())
    }

    /**
     * Retrieve the file name extension
     *
     * @param file A {@code Path} referencing a file
     * @return The file name extension (not including the dot) or an empty string if the file has not extension
     */
    def static String getExtension( Path file ) {
        assert file
        String name = file.getFileName()
        if( !name ) return ''

        int pos = name.lastIndexOf('.')
        if( pos == -1 ) return ''

        return name.substring(pos+1)
    }

    def static boolean exists(Path self, LinkOption... options) {

        return options ? Files.exists(self,options) : Files.exists(self)

    }

    def static boolean mkdir(Path self, FileAttribute<?>... attr) {

        try {
            Files.createDirectory(self,attr)
            return true
        }
        catch(IOException e) {
            return false
        }

    }

    def static boolean mkdirs(Path self, FileAttribute<?>... attr) {

        try {
            Files.createDirectories(self,attr)
            return true
        }
        catch(IOException e) {
            return false
        }
    }

    /**
     * Tests whether the application can read the file denoted by this
     * abstract pathname.
     *
     * @param self
     * @return  <code>true</code> if and only if the file specified by this
     *          abstract pathname exists <em>and</em> can be read by the
     *          application; <code>false</code> otherwise
     *
     */
    def static boolean canRead(Path self) {
        try {
            Files.isReadable(self)
        }
        catch( IOException e ) {
            log.trace("Cannot get read permission on file: $self -- Cause: ${e.getMessage()}")
            return false
        }
    }

    /**
     * Tests whether the application can modify the file denoted by this
     * abstract pathname.
     *
     * @param self
     * @return  <code>true</code> if and only if the file system actually
     *          contains a file denoted by this abstract pathname <em>and</em>
     *          the application is allowed to write to the file;
     *          <code>false</code> otherwise.
     *
     */
    def static boolean canWrite(Path self) {
        try {
            Files.isWritable(self)
        }
        catch( IOException e ) {
            log.trace("Cannot get write permission on file: $self -- Cause: ${e.getMessage()}")
            return false
        }
    }

    /**
     * Tests whether the application can execute the file denoted by this
     * abstract pathname.
     *
     * @param self
     * @return  <code>true</code> if and only if the abstract pathname exists
     *          <em>and</em> the application is allowed to execute the file
     *
     * @throws  SecurityException
     *          If a security manager exists and its <code>{@link
     *          java.lang.SecurityManager#checkExec(java.lang.String)}</code>
     *          method denies execute access to the file
     *
     */
    def static boolean canExecute(Path self) {
        try {
            Files.isExecutable(self)
        }
        catch( IOException e ) {
            log.trace("Cannot get execute permission on file: $self -- Cause: ${e.getMessage()}")
            return false
        }
    }

    /**
     * Returns the time that the file denoted by this abstract pathname was
     * last modified.
     *
     * <p> Where it is required to distinguish an I/O exception from the case
     * where {@code 0L} is returned, or where several attributes of the
     * same file are required at the same time, or where the time of last
     * access or the creation time are required, then the {@link
     * java.nio.file.Files#readAttributes(Path,Class,LinkOption[])
     * Files.readAttributes} method may be used.
     *
     * @return  A <code>long</code> value representing the time the file was
     *          last modified, measured in milliseconds since the epoch
     *          (00:00:00 GMT, January 1, 1970), or <code>0L</code> if the
     *          file does not exist or if an I/O error occurs
     *
     * @param self
     * @param options
     * @return
     */
    def static long lastModified(Path self,LinkOption...options) {
        try {
            Files.getLastModifiedTime(self,options).toMillis()
        }
        catch( IOException e ) {
            log.trace "Cannot get lastModified time on file: $self -- Cause: ${e.getMessage()}"
            return 0
        }
    }

    def static boolean isHidden(Path self) {
        Files.isHidden(self)
    }

    def static boolean isDirectory(Path self,LinkOption... options) {
        Files.isDirectory(self,options)
    }

    def static boolean isFile(Path self, LinkOption...options) {
        Files.isRegularFile(self,options)
    }

    def static boolean renameTo(Path self, Path target) {
        Files.move(self,target)
    }

    def static boolean renameTo(Path self, String target) {
        renameTo( self, FileHelper.asPath(target) )
    }

    /**
     * List the content of a given path folder. The method is semantically
     * equivalent to {@code File#list()}
     *
     * @param self The folder to list
     * @return A list of strings or {@code null} if the path is not a folder
     */
    def static String[] list(Path self) {
        listFiles(self).collect { self.toString() } as String[]
    }

    /**
     * List the content of a given path folder. The method is semantically
     * equivalent to {@code File#listFiles()}
     *
     * @param self The folder to list
     * @return A list of {@code Path} or {@code null} if the path is not a folder
     */
    def static Path[] listFiles(Path self) {

        if( !self.isDirectory() )
            return null

        def result = []
        Files.walkFileTree(self, new SimpleFileVisitor<Path>() {

            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                result.add( file )
                FileVisitResult.CONTINUE
            }

            public FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) {
                if( self == dir )
                    FileVisitResult.CONTINUE

                else {
                    result.add(dir)
                    FileVisitResult.SKIP_SUBTREE
                }
            }

        } )

        return result as Path[]

    }

    /**
     * Close a file in safe manner i.e. without throwing eventual exception raised by the *close* operation
     *
     * @param self A {@code Closable} object
     */
    static void closeQuietly( Closeable self ) {
        try {
            self.close()
        }
        catch (IOException ioe) {
            log.debug "Exception closing $self -- Cause: ${ioe.getMessage() ?: ioe.toString()}"
        }
    }


    static private DEFAULT_TAIL_BLOCK_SIZE = 5 * 1024

    /**
     * Read the last 'n' lines from a {@code SeekableByteChannel} without reading the previous content.
     * <p>
     * The method does not close the channel object, thus the caller has to close it.
     *
     *
     * @param channel The channel from where read the ending lines
     * @param n The number of lines to be read
     * @param charset The charset to use, it can be specified by using a string e.g. 'UTF-8' or a {@code Charset} object.
     *      When {@code null} is specified the current default charset is used.
     * @param blockSize The size of the block read when reading the channel tail
     * @return The 'n' lines as {@code CharSequence} object
     */
    static CharSequence tail( SeekableByteChannel channel, int n, def charset = null, int blockSize = DEFAULT_TAIL_BLOCK_SIZE ) {
        assert channel != null
        assert n > 0

        final len = channel.size()
        final buffer = ByteBuffer.allocate(blockSize)
        final result = new StringBuilder()
        long count = 0

        final charsetObj = CharsetHelper.getCharset(charset)

        long pos = len
        while( pos>0 ) {
            pos = Math.max( 0, pos-blockSize )
            channel.position(pos)
            int l = channel.read(buffer)

            buffer.flip()
            boolean isHeadLF = buffer.get(0) == LF
            boolean isTailLF = buffer.get(l-1) == LF

            def stream = new BufferedReader( new InputStreamReader(new ByteBufferBackedInputStream(buffer), charsetObj) )
            def lines = stream.readLines()
            int i = lines.size()-1
            while( i >= 0 ) {
                final ln = lines[i]

                if( isTailLF && lines.size()-1 == i && result.size() ) {
                    result.insert(0, System.lineSeparator())
                    isTailLF = false
                    count++
                }

                // prepend a new line into the result buffer
                result.insert(0,ln)

                // if there are more lines (since 'i' is greater than zero) add a new line separator
                // and increment the 'count' of total added lines
                if( i ) {
                    count++

                    // when the requested number of lines have been read, return
                    if( count < n )
                        result.insert(0, System.lineSeparator() )
                    else
                        return result
                }

                // decrement the lines pointer
                i--
            }

            // reset the buffer to read the new block
            buffer.clear()
        }

        return result
    }

    /**
     * Read the last 'n' lines from a {@code Path} without reading all the file content
     *
     * @param channel The channel from where read the ending lines
     * @param n The number of lines to be read
     * @param charset The charset to use, it can be specified by using a string e.g. 'UTF-8' or a {@code Charset} object.
     *      When {@code null} is specified the current default charset is used.
     * @param blockSize The size of the block read when reading the channel tail
     * @return The 'n' lines as {@code CharSequence} object
     */
    static CharSequence tail( Path path, int n, charset = null, int blockSize = DEFAULT_TAIL_BLOCK_SIZE) {
        def channel = Files.newByteChannel(path)
        try {
            return tail( channel, n, charset, blockSize )
        }
        finally {
            channel.closeQuietly()
        }
    }

    /**
     * Read the last 'n' lines from a {@code File} without reading all the file content
     *
     * @param channel The channel from where read the ending lines
     * @param n The number of lines to be read
     * @param charset The charset to use, it can be specified by using a string e.g. 'UTF-8' or a {@code Charset} object.
     *      When {@code null} is specified the current default charset is used.
     * @param blockSize The size of the block read when reading the channel tail
     * @return The 'n' lines as {@code CharSequence} object
     */
    static CharSequence tail( File file, int n, charset = null, int blockSize = DEFAULT_TAIL_BLOCK_SIZE ) {
        tail(file.toPath(), n, charset, blockSize)
    }

    /**
     * Read 'n' lines from the beginning of a {@code Reader} object.
     * <p>
     * The method does not close the reader object, thus the caller has to call the close method on it.
     *
     * @param reader The source object from where read the first 'n' lines
     * @param n The number of lines to read
     * @return The read lines as a {@code CharSequence} object
     */
    static CharSequence head( Reader reader, int n ) {
        assert n

        final result = new StringBuilder()
        final reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)

        def line
        int count = 0
        while( (line=reader0.readLine()) != null && count<n ) {
            if( count )
                result.append( System.lineSeparator() )
            result.append(line)
            count ++
        }

        return result
    }

    /**
     * Read 'n' lines from the beginning of a {@code InputStream} object.
     * <p>
     * The method does not close the stream object, thus the caller has to call the close method on it.
     *
     * @param reader The source object from where read the first 'n' lines
     * @param n The number of lines to read
     * @param charset The charset to use, it can be specified by using a string e.g. 'UTF-8' or a {@code Charset} object.
     *      When {@code null} is specified the current default charset is used.
     * @return The read lines as a {@code CharSequence} object
     */
    static CharSequence head( InputStream stream, int n, charset ) {
        def charsetObj = CharsetHelper.getCharset(charset)
        head( new InputStreamReader(stream, charsetObj), n)
    }

    /**
     * Read the top from the beginning of a {@code File} object, stopping as soon as
     * 'n' lines has been read.
     *
     * @param file The source file from where read the first 'n' lines
     * @param n The number of lines to read
     * @param charset The charset to use, it can be specified by using a string e.g. 'UTF-8' or a {@code Charset} object.
     *      When {@code null} is specified the current default charset is used.
     * @return The read lines as a {@code CharSequence} object
     */

    static CharSequence head( File file, int n, charset = null ) {
        head(file.toPath(), n, charset)
    }

    /**
     * Read the top from the beginning of a {@code Path} object, stopping as soon as
     * 'n' lines has been read.
     *
     * @param path The source file from where read the first 'n' lines
     * @param n The number of lines to read
     * @param charset The charset to use, it can be specified by using a string e.g. 'UTF-8' or a {@code Charset} object.
     *      When {@code null} is specified the current default charset is used.
     * @return The read lines as a {@code CharSequence} object
     */

    static CharSequence head( Path path, int n, charset = null ) {
        assert path != null
        assert n

        def charsetObj = CharsetHelper.getCharset(charset)
        def reader = Files.newBufferedReader(path,charsetObj)
        try {
            head(reader,n)
        }
        finally {
            reader.closeQuietly()
        }
    }

    /**
     * Sets the last-modified time of the file or directory named by this
     * {@code Path}.
     *
     * @param self The {@code Path} to which set the last modified time
     * @param time The new last-modified time, measured in milliseconds since
     *             the epoch (00:00:00 GMT, January 1, 1970)
     * @return <code>true</code> if and only if the operation succeeded;
     *          <code>false</code> otherwise
     */
    static boolean setLastModified(Path self, long time) {
        try {
            Files.setLastModifiedTime(self, FileTime.fromMillis(time))
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set last-modified-time: $time to path: $self"
            return false
        }
    }

    /**
     * Sets the owner's or everybody's execute permission for the specified file.
     *
     * @param self  The {@code Path} for which set the permissions
     *
     * @param executable
     *          If <code>true</code>, sets the access permission to allow execute
     *          operations; if <code>false</code> to disallow execute operations
     *
     * @param ownerOnly
     *          If <code>true</code>, the execute permission applies only to the
     *          owner's execute permission; otherwise, it applies to everybody.
     *          If the underlying file system can not distinguish the owner's
     *          execute permission from that of others, then the permission will
     *          apply to everybody, regardless of this value.
     *
     * @return  <code>true</code> if and only if the operation succeeded.  The
     *          operation will fail if the user does not have permission to
     *          change the access permissions of this file.
     *
     */
    static boolean setExecutable(Path self, boolean executable, boolean ownerOnly = true) {
        def perms = null
        try {
            perms = Files.getPosixFilePermissions(self)
            if( executable ) {
                perms.add(PosixFilePermission.OWNER_EXECUTE)
                if( !ownerOnly ) {
                    perms.add(PosixFilePermission.GROUP_EXECUTE)
                    perms.add(PosixFilePermission.OTHERS_EXECUTE)
                }
            }
            else {
                perms.remove(PosixFilePermission.OWNER_EXECUTE)
                if( !ownerOnly ) {
                    perms.remove(PosixFilePermission.GROUP_EXECUTE)
                    perms.remove(PosixFilePermission.OTHERS_EXECUTE)
                }
            }
            Files.setPosixFilePermissions(self, perms)
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set executable permissions: $perms to path: $self"
            return false
        }
    }

    /**
     * A convenience method to set the owner's read permission for the specified file
     *
     * @param self The {@code Path} for which set the readable permissions
     *
     * @param readable
     *          If <code>true</code>, sets the access permission to allow read
     *          operations; if <code>false</code> to disallow read operations
     *
     * @param ownerOnly
     *          If <code>true</code>, sets the access permission to allow read
     *          operations; if <code>false</code> to disallow read operations
     *
     * @return
     *          <code>true</code> if and only if the operation succeeded.  The
     *          operation will fail if the user does not have permission to
     *          change the access permissions of this file.
     */
    static boolean setReadable(Path self, boolean readable, boolean ownerOnly = true) {
        def perms = null
        try {
            perms = Files.getPosixFilePermissions(self)
            if( readable ) {
                perms.add(PosixFilePermission.OWNER_READ)
                if( !ownerOnly ) {
                    perms.add(PosixFilePermission.GROUP_READ)
                    perms.add(PosixFilePermission.OTHERS_READ)
                }
            }
            else {
                perms.remove(PosixFilePermission.OWNER_READ)
                if( !ownerOnly ) {
                    perms.remove(PosixFilePermission.GROUP_READ)
                    perms.remove(PosixFilePermission.OTHERS_READ)
                }
            }

            Files.setPosixFilePermissions(self, perms)
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set readable permissions: $perms to path: $self"
            return false
        }
    }

    /**
     *    * Marks the file or directory named by this abstract pathname so that
     * only read operations are allowed.
     * @param self
     * @return
     */
    static boolean setReadOnly( Path self ) {
        Set<PosixFilePermission> perms = PosixFilePermissions.fromString("r--r--r--");
        try {
            Files.setPosixFilePermissions(self, perms)
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set read-only permissions: $perms to path: $self"
            return false
        }
    }

    /**
     * Sets the owner's or everybody's write permission for the specified file.
     *
     * @param self The {@code Path} for which set the permissions
     *
     * @param writable
     *          If <code>true</code>, sets the access permission to allow write
     *          operations; if <code>false</code> to disallow write operations
     *
     * @param ownerOnly
     *          If <code>true</code>, the write permission applies only to the
     *          owner's write permission; otherwise, it applies to everybody.  If
     *          the underlying file system can not distinguish the owner's write
     *          permission from that of others, then the permission will apply to
     *          everybody, regardless of this value.
     * @return  <code>true</code> if and only if the operation succeeded.
     */
    static boolean setWritable(Path self, boolean writable, boolean ownerOnly = true ) {
        def perms = null
        try {
            perms = Files.getPosixFilePermissions(self)
            if( writable ) {
                perms.add(PosixFilePermission.OWNER_WRITE)
                if( !ownerOnly ) {
                    perms.add(PosixFilePermission.GROUP_WRITE)
                    perms.add(PosixFilePermission.OTHERS_WRITE)
                }
            }
            else {
                perms.remove(PosixFilePermission.OWNER_WRITE)
                if( !ownerOnly ) {
                    perms.remove(PosixFilePermission.GROUP_WRITE)
                    perms.remove(PosixFilePermission.OTHERS_WRITE)
                }
            }

            Files.setPosixFilePermissions(self, perms)
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set permissions: $perms to path: $self"
            return false
        }

    }

    /**
     * Get the file Unix permission as a string e.g. {@code rw-r--r--}
     *
     * @param self The {@code Path} for which the permissions string
     * @return Unix permission as a string e.g. {@code rw-r--r--}
     */
    static String getPermissions(Path self) {
        def perms = Files.getPosixFilePermissions(self)
        PosixFilePermissions.toString(perms)
    }

    /**
     * Set the file Unix permissions using a string like {@code rw-r--r--}
     *
     * @param self The {@code Path} file for which set the permissions.
     * @param permissions The permissions string e.g. {@code rw-r--r--}. It must contain 9 letters.
     */
    static boolean setPermissions( Path self, String permissions ) {
        Set<PosixFilePermission> perms = PosixFilePermissions.fromString(permissions);

        try {
            Files.setPosixFilePermissions(self, perms)
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set permissions: $permissions to path: $self"
            return false
        }
    }

    /**
     * Set the file Unix permissions using a string like {@code rw-r--r--}
     *
     * @param self The {@code File} object for which set the permissions.
     * @param permissions The permissions string e.g. {@code rw-r--r--}. It must contain 9 letters.
     */

    static boolean setPermissions( File self, String permissions ) {
        setPermissions(self.toPath(),permissions)
    }

    @PackageScope
    static digitToPerm( int value, StringBuilder sb = new StringBuilder() ) {
        assert value >= 0 && value < 8

        def x = (value & 1) == 1
        def w = (value & 2) == 2
        def r = (value & 4) == 4

        if (r) {
            sb.append('r');
        } else {
            sb.append('-');
        }
        if (w) {
            sb.append('w');
        } else {
            sb.append('-');
        }
        if (x) {
            sb.append('x');
        } else {
            sb.append('-');
        }

    }

    /**
     * Set the file Unix permission using the digit representing the respectively
     * the permissions for the owner, the group and others.
     *
     * @link http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation
     *
     * @param self The {@code Path} file for which set the permissions
     * @param owner The owner permissions using a octal numeric representation.
     * @param group The group permissions using a octal numeric representation.
     * @param other The others permissions using a octal numeric representation.
     */
    static boolean setPermissions( Path self, int owner, int group, int other ) {
        def str = new StringBuilder()
        digitToPerm(owner, str)
        digitToPerm(group, str)
        digitToPerm(other, str)
        Set<PosixFilePermission> perms = PosixFilePermissions.fromString(str.toString());
        try {
            Files.setPosixFilePermissions(self, perms)
            return true
        }
        catch( IOException e ) {
            log.debug "Unable to set permissions: $perms to path: $self"
            return false
        }
    }

    /**
     * Set the file Unix permission using the digit representing the respectively
     * the permissions for the owner, the group and others.
     *
     * @link http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation
     *
     * @param self The {@code File} object for which set the permissions
     * @param owner The owner permissions using a octal numeric representation.
     * @param group The group permissions using a octal numeric representation.
     * @param other The others permissions using a octal numeric representation.
     */
    static boolean setPermissions( File self, int owner, int group, int other ) {
        setPermissions(self.toPath(), owner, group, other)
    }


    static void deleteOnExit(Path self) {
        Runtime.getRuntime().addShutdownHook { self.delete() }
    }

    /**
     * Resolve a symbolic link to the actual target file, in a similar manner to the
     * Linux {@code readlink -f} command.
     * <p>
     * It canonicalize by following every symlink in every component of the given name recursively;
     * all but the last component must exist.
     *
     * @param self
     * @return The resolved target path of the source path itself it is not a file or the file does not exist
     */
    static Path resolveSymLink(Path self) {
        try {
            def path = Files.readSymbolicLink(self)
            return resolveSymLink(path)
        }
        catch( NotLinkException | NoSuchFileException e ) {
            return self
        }
    }

    static File resolveSymLink(File self) {
        self.toPath().resolveSymLink().toFile()
    }

    /**
     * Create a folder if not already exists
     * @param target The folder to be created
     */
    static public void createDirIfNotExists( Path target ) {
        assert target

        try {
            if( !Files.readAttributes(target, BasicFileAttributes).isDirectory() )
                throw new IOException("Cannot create folder: $target -- A file with the same name already exists")
        }
        catch (IOException e) {
            if( !Files.createDirectories(target) )
                throw new IOException("Cannot create folder: $target -- Check file systeme access permission")

        }
    }

    static public void createDirIfNotExists( File target ) {
        createDirIfNotExists(target.toPath())
    }

}
