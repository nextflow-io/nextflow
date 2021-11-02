package nextflow.sraql

import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Implements SQL plugin for Nextflow
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
class SraqlPlugin extends BasePlugin {

    SraqlPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }
}
