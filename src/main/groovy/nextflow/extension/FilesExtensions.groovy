package nextflow.extension
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.SimpleFileVisitor
import java.nio.file.StandardCopyOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute

import groovy.util.logging.Slf4j
import nextflow.util.FileHelper
/**
 * Add utility methods to standard classes {@code File} and {@code Path}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class FilesExtensions {

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

    def static delete( Path path ) {
        Files.delete(path)
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
    def static copyTo( Path source, Path target ) {

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
        copyTo(source.toPath(), Paths.get(target)).toFile()
    }

    /**
     * Copy or a file or a directory. It mimics the semantic of the Linux *cp* command
     *
     * @param source
     * @param target
     * @return
     */
    def static Path copyTo( Path source, String target ) {
        copyTo(source, Paths.get(target))
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
        moveTo(source.toPath(), Paths.get(target)).toFile()
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
        moveTo(source, target as Path )
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

    def static boolean canRead(Path self) {
        Files.isReadable(self)
    }

    def static boolean canWrite(Path self) {
        Files.isWritable(self)
    }

    def static boolean canExecute(Path self) {
        Files.isExecutable(self)
    }

    def static long lastModified(Path self,LinkOption...options) {
        Files.getLastModifiedTime(self,options).toMillis()
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

    // TODO implements using the new API
    def static boolean renameTo(Path self, Path target) {
            self.toFile().renameTo(target.toFile())
    }

    def static boolean renameTo(Path self, String target) {
        renameTo( self, target as Path )
    }

    // TODO implements using the new API
    def static String[] list(Path self) {
        self.toFile().list()
    }

    // TODO implements using the new API
    def static Path[] listFiles(Path self) {
        self.toFile().listFiles().collect { it.toPath() } as Path[]
    }


    static closeQuietly( Closeable self ) {
        try {
            self.close()
        }
        catch (IOException ioe) {
            log.debug "Exception closing $self -- Cause: ${ioe.getMessage() ?: ioe.toString()}"
        }
    }

}
