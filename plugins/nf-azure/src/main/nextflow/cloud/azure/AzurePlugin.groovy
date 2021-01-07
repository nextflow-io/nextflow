package nextflow.cloud.azure

import groovy.transform.CompileStatic
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
}
