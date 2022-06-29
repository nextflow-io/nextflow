package nextflow.file

import java.nio.file.CopyOption
import java.nio.file.Path


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
interface FileSystemProviderExt {

    boolean canCopy(Path source, Path target, CopyOption... options)

    void copyToForeignTarget(Path source, Path target, CopyOption... options) throws IOException

}