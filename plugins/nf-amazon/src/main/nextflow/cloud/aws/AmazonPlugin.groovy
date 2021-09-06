package nextflow.cloud.aws


import com.upplication.s3fs.S3FileSystemProvider
import com.upplication.s3fs.S3OutputStream
import groovy.transform.CompileStatic
import nextflow.file.FileHelper
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper
/**
 * Nextflow plugin for Amazon extensions
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AmazonPlugin extends BasePlugin {

    AmazonPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @Override
    void start() {
        super.start()
        FileHelper.getOrInstallProvider(S3FileSystemProvider)
    }

    @Override
    void stop() {
        super.stop()
        // shutdown s3 uploader
        S3OutputStream.shutdownExecutor()
    }
}
