package nextflow.file

import java.nio.file.Path

import org.pf4j.ExtensionPoint

/**
 * Declares the Bash helper to support third-party file storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface BashFunExt extends ExtensionPoint {

    /**
     * Bash helper library that implements the support for third-party storage
     * to be included in the job command wrapper script
     *
     * @return The Bash snippet implementing the support for third-party such as AWS S3
     * or {@code null} if not supported
     */
    String helperLib()

    /**
     * The name of a Bash helper function to upload a file to a remote file storage
     * 
     * @return The name of the upload function or {@code null} if not supported
     */
    String uploadCmd(String source, Path target)

}
