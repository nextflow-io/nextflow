package nextflow.cloud.aws.util

import java.nio.file.Path

import com.upplication.s3fs.S3Path
import nextflow.cloud.aws.batch.AwsBatchFileCopyStrategy
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory
/**
 * Implements the a factory strategy to parse and build S3 path URIs
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3PathFactory extends FileSystemPathFactory {

    @Override
    protected Path parseUri(String str) {
        // normalise 's3' path
        if( str.startsWith('s3://') && str[5]!='/' ) {
            final path = "s3:///${str.substring(5)}"
            // note: this URI constructor parse the path parameter and extract the `scheme` and `authority` components
            final uri = new URI(null,null, path,null,null)
            return FileHelper.getOrCreateFileSystemFor(uri).provider().getPath(uri)
        }
        return null
    }

    @Override
    protected String toUriString(Path path) {
        return path instanceof S3Path ? "s3:/$path".toString() : null
    }

    @Override
    protected String getBashLib(Path target) {
        return S3BashLib.script()
    }

    @Override
    protected String getUploadCmd(String source, Path target) {
        return target instanceof S3Path
                ? AwsBatchFileCopyStrategy.uploadCmd(source,target)
                : null
    }

}
