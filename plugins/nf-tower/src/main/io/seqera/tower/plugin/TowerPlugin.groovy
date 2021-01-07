package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper
/**
 * Nextflow Tower plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TowerPlugin extends BasePlugin {

    TowerPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

}
