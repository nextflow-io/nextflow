package nextflow.extension

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.PosixFileAttributeView
import java.nio.file.attribute.PosixFileAttributes
import java.nio.file.attribute.PosixFilePermission

import nextflow.util.FileHelper
import org.apache.commons.io.FileUtils
import org.apache.commons.io.FilenameUtils

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilesExtensions {

    /**
     * Check if a file - or - a directory is empty
     *
     * @param file The file under test
     * @return {@code true} if the file does not exist or it is empty
     */
    @Deprecated
    def static boolean isEmpty( File file ) {
        FileHelper.isEmpty(file)
    }


    /**
     * TODO This is not working because the concrete class (UnixPath) is providing its own isEmpty() method ..
     *
     * Check if a file - or - a directory is empty
     *
     * @param file The file under test
     * @return {@code true} if the file does not exist or it is empty
     */
    @Deprecated
    def static boolean isEmpty( Path path ) {
        FileHelper.isEmpty(path)
    }

    /**
     * Check if a file - or - a folder is not empty
     *
     * @param file
     * @return
     */
    @Deprecated
    def static isNotEmpty( File file ) {
        return FileHelper.isNotEmpty(file)
    }

    @Deprecated
    def static isNotEmpty( Path file ) {
        return FileHelper.isNotEmpty(file)
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
     * Copy a file  -or- a whole directory to a new file -or- a new directory
     *
     * @param source
     * @param target
     * @return
     */
    def static copyTo( File source, File target ) {
        if( source.isDirectory() ) {
            FileUtils.copyDirectory(source, target)
        }
        else {
            FileUtils.copyFile(source, target)
        }
    }

    /*
     * TODO This MUST be reimplemented without using {@code Path#toFile} adapter method
     */
    def static copyTo( Path source, Path target ) {
        copyTo(source.toFile(), target.toFile())
    }

    /**
     * Copy a file  -or- a whole directory to a new file -or- a new directory
     *
     * @param source
     * @param target
     * @return
     */
    def static copyTo( File source, String target ) {
        copyTo(source, new File(target))
    }

    def static copyTo( Path source, String target ) {
        copyTo(source, Paths.get(target))
    }


    def static moveTo( File source, File target ) {
        if( source.isDirectory() ) {
            FileUtils.moveDirectory(source, target)
        }
        else {
            FileUtils.moveFile(source, target)
        }
    }

    def static moveTo( File source, String target ) {
        moveTo(source, new File(target))
    }

    /*
     * TODO This MUST be reimplemented without using {@code Path#toFile} adapter method
     */
    def static moveTo( Path source, Path target ) {
        moveTo(source.toFile(), target.toFile())
    }

    def static moveTo( Path source, String target ) {
        moveTo(source, Paths.get(target))
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
        assert file
        FilenameUtils.getBaseName(file.name)
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
        assert file
        FilenameUtils.getExtension(file.name)
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

        return Files.exists(self,options)

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

    def static void setExecutable(Path self, boolean flag) {

        PosixFileAttributeView view = Files.getFileAttributeView(self,PosixFileAttributeView.class);
        if( view == null ) {
            throw new UnsupportedOperationException()
        }

        PosixFileAttributes attributes = view.readAttributes();
        Set<PosixFilePermission> permissions = attributes.permissions();
        if( flag )  {
            permissions.add(PosixFilePermission.GROUP_EXECUTE)
        }
        else {
            permissions.remove(PosixFilePermission.GROUP_EXECUTE);
        }
        // finally set the permission
        view.setPermissions(permissions);
    }

}
