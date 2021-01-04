package nextflow.ga4gh

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper
/**
 * Nextflow plugin for GA4GH extensions
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Ga4ghPlugin extends BasePlugin {

    Ga4ghPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

}
