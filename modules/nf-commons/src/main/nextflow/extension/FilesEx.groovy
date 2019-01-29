/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.extension

import org.codehaus.groovy.runtime.InvokerHelper
import static java.nio.file.StandardCopyOption.REPLACE_EXISTING

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileSystemException
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.NotLinkException
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileTime
import java.nio.file.attribute.PosixFilePermission
import java.nio.file.attribute.PosixFilePermissions

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.FromString
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.io.ByteBufferBackedInputStream
import nextflow.util.CharsetHelper
import nextflow.util.CheckHelper
/**
 * Add utility methods to standard classes {@code File} and {@code Path}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class FilesEx {

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

    static boolean delete( Path path ) {
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
    static boolean deleteDir(Path path) {
        def attr = FileHelper.readAttributes(path)
        if( !attr )
             return true

        if( !attr.isDirectory() )
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
     * Copy or a file or a directory content.
     * It mimics the semantic of Linux {@code cp -r <source> <target>} command.
     *
     * @param source The source file or directory
     * @param target The target file or directory
     * @return The target {@link Path} file
     */
    static File copyTo( File source, File target ) {
        copyTo(source.toPath(), target.toPath()).toFile()
    }

    /**
     * Copy a file or a directory content.
     * It mimics the semantic of Linux {@code cp -r <source> <target>} command.
     *
     * @param source The source file or directory
     * @param target The target file or directory
     * @return The target {@link Path} file
     */
    static Path copyTo( Path source, Path target ) {
        assert source
        assert target

        if( source.isDirectory() ) {
            boolean targetExists
            boolean targetIsDirectory
            try {
                targetIsDirectory = Files.readAttributes(target, BasicFileAttributes).isDirectory()
                targetExists = true
            }
            catch( Exception e ) {
                targetExists = false
                targetIsDirectory = false
            }

            // when the target path is a directory
            // copy the source path into the target folder
            if( targetIsDirectory )
                return FileHelper.copyPath(source, target.resolve(source.getName()))

            // otherwise it is a file => we cannot write a directory with the same name
            if( targetExists )
                throw new FileAlreadyExistsException("Cannot copy directory -- A file with the same name already exists: $target")

            def parent = target.getParent()
            if( parent && !parent.exists() ) parent.mkdirs()
            return FileHelper.copyPath(source,target)
        }

        if( target.isDirectory() ) {
            target = target.resolve(source.getName())
            return FileHelper.copyPath(source, target, REPLACE_EXISTING)
        }

        // create the target parent directories if do not exist
        def parent = target.getParent()
        if( parent && !parent.exists() ) {
            parent.mkdirs()
        }

        return FileHelper.copyPath(source, target, REPLACE_EXISTING)
    }


    /**
     * Copy or a file or a directory content.
     * It mimics the semantic of Linux {@code cp -r <source> <target>} command.
     *
     * @param source The source file or directory
     * @param target The target file or directory
     * @return The target {@link Path} file
     */
    static File copyTo( File source, String target ) {
        copyTo(source.toPath(), FileHelper.asPath(target)).toFile()
    }

    /**
     * Copy or a file or a directory content.
     * It mimics the semantic of Linux {@code cp -r <source> <target>} command.
     *
     * @param source The source file or directory
     * @param target The target file or directory
     * @return The target {@link Path} file
     */
    static Path copyTo( Path source, String target ) {
        copyTo(source, FileHelper.asPath(target))
    }

    /**
     * Move a file or a directory. Mimics the Linux *mv* command
     *
     * @param source
     * @param target
     * @return
     */
    static File moveTo( File source, File target ) {
        moveTo(source.toPath(), target.toPath()).toFile()
    }

    /**
     * Move a file or a directory. Mimics the Linux *mv* command
     *
     * @param source
     * @param target
     * @return
     */
    static File moveTo( File source, String target ) {
        moveTo(source.toPath(), FileHelper.asPath(target)).toFile()
    }

    /**
     * Move a file or a directory. Mimics the Linux *mv* command
     *
     * @param source
     * @param target
     * @return
     */
    static Path moveTo( Path source, Path target ) {

        if( source.isDirectory() ) {
            if( target.isDirectory() ) {
                // when the target path is a directory
                // move the source path into the target folder
                return FileHelper.movePath(source, target.resolve(source.getName()), REPLACE_EXISTING)
            }
            else {
                def parent = target.getParent()
                if( parent && !parent.exists()) { parent.mkdirs() }
                return FileHelper.movePath(source, target, REPLACE_EXISTING)
            }
        }

        // when the target path is a directory
        // move the source file into the target using the same name
        if( target.isDirectory() ) {
            target = target.resolve(source.getName())
            return FileHelper.movePath(source, target, REPLACE_EXISTING)
        }

        // create the parent directories if do not exist
        def parent = target.getParent()
        if( parent && !parent.exists() ) { parent.mkdirs() }

        // move the source file to the target file,
        // overwriting the target if exists
        return FileHelper.movePath(source, target, REPLACE_EXISTING)
    }

    static Path moveTo( Path source, String target ) {
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
     * </pre>
     *
     * The output will be the same irrespective of the machine that the code is running on.
     * @param file The filename to query, null returns null
     * @param times The number of times it checks for the extension to be removed (useful for files with multiple extensions)
     * @return The name of the file without the path, or an empty string if none exists
     */
    static String getBaseName( File file, int times=1 ) {
        getBaseName(file.toPath(), times)
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
     * </pre>
     *
     * The output will be the same irrespective of the machine that the code is running on.
     * @param file The filename to query, null returns null
     * @param times The number of times it checks for the extension to be removed (useful for files with multiple extensions)
     * @return The name of the file without the path, or an empty string if none exists
     */
    static String getBaseName( Path self, int times=1 ) {
        assert self

        String name = self.getFileName()
        if( !name ) return ''

        while( times-- > 0 ) {
            int pos = name.lastIndexOf('.')
            if( pos == -1 )
                break
            name = name.substring(0,pos)
        }

        return name.toString()
    }

    /**
     * Extend {@code File} adding a {@code getSimpleName()} method which returns the
     * file base name removing all file extension.
     *
     * When a file start starts with a dot character it returns the first name token
     * after the dot.
     *
     * For example:
     *
     *  <pre>
     *   a.txt     --> a
     *   a.tar.gz  --> a
     *   a/b/c     --> c
     *   a/b/c.txt --> c
     *   /         --> null
     * </pre>
     *
     * @param self The file {@code File}
     * @return The file simple name string
     */
    static String getSimpleName( File self ) {
        return getSimpleName(self.toPath())
    }

    /**
     * Extend {@code Path} adding a {@code getSimpleName()} method which returns the
     * file base name removing all file extension.
     *
     * When a file start starts with a dot character it returns the first name token
     * after the dot.
     *
     * For example:
     *
     *  <pre>
     *   a.txt     --> a
     *   a.tar.gz  --> a
     *   a/b/c     --> c
     *   a/b/c.txt --> c
     *   /         --> null
     * </pre>
     *
     * @param self The file {@code Path}
     * @return The file simple name string
     */
    static String getSimpleName( Path self ) {
        assert self

        String name = self.getFileName()
        if( name?.size() <= 1 )
            return name

        int pos = name.substring(1).indexOf('.')
        pos != -1 ? name.substring(0,pos+1) : name
    }

    /**
     * Extend {@code Path} adding a getName method semantically equivalent to {@code File#getName}
     *
     * @return The name of the path as string
     */
    static String getName( Path self ) {
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
    static String getExtension( File file ) {
        getExtension(file.toPath())
    }

    /**
     * Retrieve the file name extension
     *
     * @param file A {@code Path} referencing a file
     * @return The file name extension (not including the dot) or an empty string if the file has not extension
     */
    static String getExtension( Path file ) {
        assert file
        String name = file.getFileName()
        if( !name ) return ''

        int pos = name.lastIndexOf('.')
        if( pos == -1 ) return ''

        return name.substring(pos+1)
    }

    static boolean exists(Path self, LinkOption... options) {

        return options ? Files.exists(self,options) : Files.exists(self)

    }

    /**
     * Creates the directory named by this abstract pathname
     *
     * @param self The directory to be created
     * @param attr
     *          an optional list of file attributes to set atomically when
     *          creating the directory
     * @return {@code true} if the directory is created successfully, or {@code false} otherwise
     */
    static boolean mkdir(Path self, FileAttribute<?>... attr) {

        try {
            Files.createDirectory(self,attr)
            return true
        }
        catch(IOException e) {
            return false
        }

    }

    /**
     * Creates the directory named by this abstract pathname, including any necessary but nonexistent parent directories
     *
     * @param self The path to be created
     * @param attr
     *          an optional list of file attributes to set atomically when
     *          creating the directory
     * @return {@code true} if the directory is created successfully, or {@code false} otherwise
     */
    static boolean mkdirs(Path self, FileAttribute<?>... attr) {

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
    static boolean canRead(Path self) {
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
    static boolean canWrite(Path self) {
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
    static boolean canExecute(Path self) {
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
    static long lastModified(Path self,LinkOption...options) {
        try {
            Files.getLastModifiedTime(self,options).toMillis()
        }
        catch( IOException e ) {
            log.trace "Cannot get lastModified time on file: $self -- Cause: ${e.getMessage()}"
            return 0
        }
    }

    static boolean isHidden(Path self) {
        Files.isHidden(self)
    }

    static boolean isDirectory(Path self,LinkOption... options) {
        Files.isDirectory(self,options)
    }

    static boolean isFile(Path self, LinkOption...options) {
        Files.isRegularFile(self,options)
    }

    static boolean isLink(Path self) {
        Files.isSymbolicLink(self)
    }

    static boolean renameTo(Path self, Path target) {
        Files.move(self,target)
    }

    static boolean renameTo(Path self, String target) {
        renameTo( self, FileHelper.asPath(target) )
    }

    /**
     * List the content of a given path folder. The method is semantically
     * equivalent to {@code File#list()}
     *
     * @param self The folder to list
     * @return A list of strings or {@code null} if the path is not a folder
     */
    static String[] list(Path self) {
        listFiles(self).collect { getName(it) } as String[]
    }

    /**
     * List the content of a given path folder. The method is semantically
     * equivalent to {@code File#listFiles()}
     *
     * @param self The folder to list
     * @return A list of {@code Path} or {@code null} if the path is not a folder
     */
    static Path[] listFiles(Path self, @ClosureParams(value = FromString.class, options = ["java.nio.file.Path", "java.nio.file.Path,java.nio.file.attribute.BasicFileAttributes"]) Closure<Boolean> filter=null) {

        if( !self.isDirectory() )
            return null

        def result = []
        Files.walkFileTree(self, new SimpleFileVisitor<Path>() {

            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                if( filter==null || invokeFilter(filter,file,attrs) )
                    result.add( file )
                FileVisitResult.CONTINUE
            }

            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) {
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

    @CompileStatic
    static private boolean invokeFilter(Closure<Boolean> filter, Path file, BasicFileAttributes attrs) {
        def params = filter.maximumNumberOfParameters
        if( params==1 ) {
            filter.call(file)
        }
        else if( params==2 ) {
            filter.call(file,attrs)
        }
        else
            throw new IllegalArgumentException("Path `listFiles` filter closure cannot take more than 2 parameters")
    }

    /**
     * Close a file in safe manner i.e. without throwing eventual exception raised by the *close* operation
     *
     * @param self A {@code Closable} object
     */
    static void closeQuietly( Closeable self ) {
        try {
            if(self) self.close()
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
        def channel = null
        try {
            try {
                channel = Files.newByteChannel(path)
                return tail( channel, n, charset, blockSize )
            }
            catch( UnsupportedOperationException e ) {
                log.trace "Unable to open as a byte-channel file: $path -- trying as inputstream"
                channel = Files.newInputStream(path)
                return tail( new InputStreamReader(channel), n )
            }
        }
        finally {
            channel?.closeQuietly()
        }
    }

    /**
     * Read the last 'n' lines from a {@code Path} without reading all the file content
     *
     * @param channel The channel from where read the ending lines
     * @param n The number of lines to be read
     * @return The 'n' lines as {@code CharSequence} object
     */
    static CharSequence tail( Reader reader, int n ) {
        String line
        String[] buffer = new String[n]
        int index = 0
        while( (line = reader.readLine()) != null ) {
            buffer[ (index++) % n ] = line
        }

        StringBuilder result = new StringBuilder()
        int count = Math.min(n, index)
        while( true ) {
            result.insert( 0, buffer[ (--index) % n ] )
            if( --count )
                result.insert( 0, System.lineSeparator())
            else
                break
        }

        result
    }


    static CharSequence tail( InputStream stream, int n, charset ) {
        final charsetObj = CharsetHelper.getCharset(charset)
        tail( new InputStreamReader(stream,charsetObj), n )
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
    static void createDirIfNotExists( Path target ) {
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

    static void createDirIfNotExists( File target ) {
        createDirIfNotExists(target.toPath())
    }

    static Path plus( Path path, String other ) {
        if( !other )
            return path

        def parent = path.getParent()
        if( parent ) {
            parent.resolve( path.getName() + other )
        }
        else {
            path.getFileSystem().getPath( path.toString() + other )
        }
    }

    static Path plus( Path path, Path other ) {
        plus(path, other?.toString())
    }

    static Path div( Path path, String other ) {
        if( !other )
            return path

        path.resolve(other)
    }

    static Path div( Path path, Path other ) {
        if( !other )
            return path

        path.resolve(other)
    }

    static Path minus(Path path, int i) {
        def result = path
        while( i-- > 0 && result != null )
            result = result.getParent()

        return result ?: path.getRoot()
    }
//
//    static Path minus( Path self, Path other ) {
//        self.relativize(other)
//    }
//
//    static Path minus( Path self, String other ) {
//        def x = self.getFileSystem().getPath(other)
//        println x.class
//        println self.class
//        self.relativize(x)
//    }

    static Path or( Path path, String other ) {
        path.resolveSibling(other)
    }

    static Path or( Path path, Path other ) {
        path.resolveSibling(other)
    }


    /**
     * Roll a file moving to a new path whose name ends with .1
     *
     * @param self The file itself
     * @param maxRolls Max number of time to apply the file rolling (i.e. rename to a new name)
     *
     */
    static void rollFile( Path self, int maxRolls=9 ) {

        if( FilesEx.exists(self) ) {
            def newName = self.getFileName().toString() +'.1'
            rollFile0(self, self.resolveSibling(newName), maxRolls)
        }

    }


    @PackageScope
    static void rollFile0( Path self, Path newFile, int max, int depth=0 ) {
        assert self
        assert newFile

        if( newFile.exists() ) {
            def extension = newFile.getExtension()
            if( !extension?.isInteger() ) {
                throw new IllegalArgumentException("Unexpected rolling file extension: $extension")
            }

            if( ++depth < max ) {
                def nextName = "${newFile.getBaseName()}.${ extension.toInteger() +1 }"
                rollFile0(newFile, newFile.resolveSibling(nextName), max, depth)
            }
            else{
                Files.delete(newFile)
            }
        }

        FilesEx.renameTo(self, newFile)

    }

    static private LINK_PARAMS = [hard:Boolean, overwrite: Boolean]

    /**
     * Make a symbolic or a hard file system link
     *
     * @param existing The existing file
     * @param opts
     *      Optional parameters:
     *      {@code hard}: when {@code true} create a hard link, otherwise it creates a symbolic link
     *      {@code overwrite}
     * @param link The {@link Path} of the link to create
     */
    static Path mklink( Path existing, Map opts, Path link ) {
        CheckHelper.checkParams('mklink', opts, LINK_PARAMS)

        final hard = opts?.hard == true ?: false
        final overwrite = opts?.overwrite == true ?: false

        if( existing == link )
            throw new IllegalArgumentException("Cannot create file link -- Source and target are the same file: $existing")

        try {
            createLinkImpl(existing, link, hard)
        }
        catch( FileAlreadyExistsException e ) {
            if( !overwrite ) throw e
            Files.delete(link)
            createLinkImpl(existing, link, hard)
        }

        return link
    }

    static private void createLinkImpl(Path existing, Path link, boolean hard) {
        if( hard ) {
            try {
                Files.createLink(link, existing)
            }
            catch( FileSystemException e ) {
                if( !Files.isDirectory(existing) ) {
                    throw e
                }

                // create the target link as a directory
                Files.createDirectory(link)

                // then create a link for each entry in the `existing` folder
                existing
                        .listFiles()
                        .each { file -> createLinkImpl( file, link.resolve(file.getFileName()), true) }
            }
        }
        else {
            Files.createSymbolicLink(link, existing)
        }
    }

    static Path mklink( Path existing, Map opts, File link ) {
        mklink(existing, opts, link.toPath())
    }

    static Path mklink( Path existing, Map opts, String link ) {
        mklink(existing, opts, existing.getFileSystem().getPath(link))
    }

    static Path mklink( Path existing, Path link ) {
        mklink(existing, null, link )
    }

    static Path mklink( Path existing, File link ) {
        mklink(existing, null, link.toPath())
    }

    static Path mklink( Path existing, String link ) {
        mklink(existing, null, existing.getFileSystem().getPath(link))
    }

    static Path mklink( Path existing, Map opts = null ) {
        mklink(existing,opts,existing.getFileName())
    }

    /**
     * A shortcut for {@link Path#toAbsolutePath()} and {@link Path#normalize()}
     *
     * @param self
     * @return An absolute and normalized path
     */
    static Path complete( Path self ) {
        self.toAbsolutePath().normalize()
    }

    @Deprecated
    static BasicFileAttributes readAttributes(Path path) {
        log.warn "Method `readAttributes` is deprecated"
        try {
            Files.readAttributes(path,BasicFileAttributes)
        }
        catch( IOException e ) {
            log.trace "Unable to read attributes for file: $path"
            return null
        }
    }

    static boolean matches( Path self, String pattern ) {
        FileHelper.getPathMatcherFor("glob:$pattern", self.fileSystem).matches(self)
    }

    static boolean matches( File self, String pattern ) {
        matches(self.toPath(), pattern)
    }

    static URI getUri( Path self ) {
        self.toUri()
    }

    static URI getUri( File self ) {
        self.toURI()
    }

    static String toUriString( Path path ) {
        if(path==null)
            return null
        final scheme = getScheme(path)
        if( scheme == 'file' )
            return path.toString()
        if( scheme == 's3' )
            return "$scheme:/$path".toString()
        if( scheme == 'gs' ) {
            final bucket = InvokerHelper.invokeMethod(path, 'bucket', InvokerHelper.EMPTY_ARGS)
            return "$scheme://$bucket$path".toString()
        }
        return path.toUri().toString()
    }

    static String getScheme(Path path) {
        path.getFileSystem().provider().getScheme()
    }
}
