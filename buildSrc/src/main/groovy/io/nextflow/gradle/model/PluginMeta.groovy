package io.nextflow.gradle.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a Nextflow plugin meta-info
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
class PluginMeta {
    String id
    String name
    String provider
    String description
    List<PluginRelease> releases
}
