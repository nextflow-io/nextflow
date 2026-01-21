package nextflow.util

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute

@Slf4j
@CompileStatic
class BasicFileHelper {
    /**
     * Check whenever a file or a directory is empty
     *
     * @param file The file path to
     */
    static boolean empty( File file ) {
        assert file
        empty(file.toPath())
    }

    static boolean empty(Path path ) {

        def attrs
        try {
            attrs = Files.readAttributes(path, BasicFileAttributes.class)
        }
        catch (IOException e) {
            return true;
        }

        if ( attrs.isDirectory() ) {
            def stream = Files.newDirectoryStream(path)
            try {
                Iterator<Path> itr = stream.iterator()
                return !itr.hasNext()
            }
            finally {
                stream.close()
            }
        }
        else {
            return attrs.size() == 0
        }

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
            log.debug "Failed to create directory '$self'", e
            return false
        }
    }
}
