package io.nextflow.gradle.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a list of Nextflow plugins
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString
@EqualsAndHashCode
@CompileStatic
class PluginsRepo {

    List<PluginMeta> plugins

}
