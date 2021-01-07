package nextflow.config

import groovy.transform.CompileStatic

/**
 * Model a mini-dsl for plugins configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PluginsDsl {

    private Set<String> plugins = []

    Set<String> getPlugins() { plugins }

    void id( String plg ) {
        if( !plg )
            throw new IllegalArgumentException("Plugin id cannot be empty or null")
        plugins << plg
    }

}
