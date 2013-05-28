package nextflow.extension

import nextflow.util.FileHelper
import org.apache.commons.io.FileUtils
import org.apache.commons.io.FilenameUtils

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OtherExtensions {

    /**
     * Check if a file - or - a directory is empty
     *
     * @param file The file under test
     * @return {@code true} if the file does not exist or it is empty
     */
    def static boolean isEmpty( File file ) {
        FileHelper.isEmpty(file)
    }

    /**
     * Check if a file - or - a folder is not empty
     *
     * @param file
     * @return
     */
    def static isNotEmpty( File file ) {
        return FileHelper.isNotEmpty(file)
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


    def static String rightTrim(String self) {
        self.replaceAll(/\s+$/,"")
    }

    def static String leftTrim( String self ) {
        self.replaceAll(/^\s+/,"")
    }
}
