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
            return create(path)
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

    /**
     * Creates a {@link S3Path} from a S3 formatted URI.
     *
     * @param path
     *      A S3 URI path e.g. s3:///BUCKET_NAME/some/data.
     *      NOTE it expect the s3 prefix provided with triple `/` .
     *      This is required by the underlying implementation expecting the host name in the URI to be empty
     *      and the bucket name to be the first path element
     * @return
     *      The corresponding {@link S3Path}
     */
    static S3Path create(String path) {
        if( !path ) throw new IllegalArgumentException("Missing S3 path argument")
        if( !path.startsWith('s3:///') ) throw new IllegalArgumentException("S3 path must start with s3:/// prefix -- offending value '$path'")
        // note: this URI constructor parse the path parameter and extract the `scheme` and `authority` components
        final uri = new URI(null,null, path,null,null)
        return (S3Path)FileHelper.getOrCreateFileSystemFor(uri).provider().getPath(uri)
    }
}
