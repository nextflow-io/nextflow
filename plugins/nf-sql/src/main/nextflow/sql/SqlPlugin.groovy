package nextflow.sql

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Nextflow SQL plugin
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@CompileStatic
class SqlPlugin extends BasePlugin {

    SqlPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

}
