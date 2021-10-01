package nextflow.sql

import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Implements SQL plugin for Nextflow
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SqlPlugin extends BasePlugin {

    SqlPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }
}
