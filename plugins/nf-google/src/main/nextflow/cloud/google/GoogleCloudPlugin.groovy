package nextflow.cloud.google

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper
/**
 * Implement the plugin entry point for Google Cloud support
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GoogleCloudPlugin extends BasePlugin {

    GoogleCloudPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

}
