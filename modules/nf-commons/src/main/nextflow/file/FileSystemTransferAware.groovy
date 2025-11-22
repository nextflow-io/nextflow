package nextflow.file

import java.nio.file.CopyOption
import java.nio.file.Path


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
interface FileSystemTransferAware {

    boolean canUpload(Path source, Path target)

    boolean canDownload(Path source, Path target)

    void download(Path source, Path target, CopyOption... options) throws IOException

    void upload(Path source, Path target, CopyOption... options) throws IOException

}
