package nextflow.plugin

import groovy.transform.CompileStatic
import org.pf4j.update.UpdateRepository

/**
 * Extension to pf4j's UpdateRepository which supports pre-fetching
 * metadata for a specified set of plugins.
 *
 * This gives the ability to avoid downloading metadata for unused
 * plugins.
 */
@CompileStatic
interface PrefetchUpdateRepository extends UpdateRepository {
    /**
     * This will be called when Nextflow starts, before
     * initialising the plugins.
     */
    void prefetch(List<PluginSpec> plugins)
}
