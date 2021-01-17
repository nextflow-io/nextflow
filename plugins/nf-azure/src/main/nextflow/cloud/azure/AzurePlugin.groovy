package nextflow.cloud.azure

import com.azure.storage.blob.nio.AzureFileSystemProvider
import groovy.transform.CompileStatic
import nextflow.file.FileHelper
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Azure cloud plugin for Nextflow
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzurePlugin extends BasePlugin {

    AzurePlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @Override
    void start() {
        super.start()
        // register Azure file system
        FileHelper.getOrInstallProvider(AzureFileSystemProvider)
    }
}
