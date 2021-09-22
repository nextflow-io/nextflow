package nextflow.ignite

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper
/**
 * Implements Apache Ignite plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class IgnitePlugin extends BasePlugin {

    IgnitePlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

}
