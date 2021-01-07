/*
 * Copyright 2020, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package nextflow.plugin

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.extension.Bolts
import nextflow.extension.FilesEx
import nextflow.util.CacheHelper
import org.pf4j.DefaultPluginManager
import org.pf4j.PluginManager
import org.pf4j.PluginStateEvent
import org.pf4j.PluginStateListener
/**
 * Manage plugins installation and configuration
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PluginsFacade implements PluginStateListener {

    private Path PLUGINS_LOCAL_ROOT = Paths.get('.nextflow/plr')

    private Map<String,String> env = new HashMap<>(System.getenv())

    private String mode
    private Path root
    private PluginUpdater updater
    private CustomPluginManager manager
    private DefaultPlugins defaultPlugins
    private String indexUrl = Plugins.DEFAULT_PLUGINS_REPO

    PluginsFacade() {
        mode = getPluginsMode()
        root = getPluginsDir()
        System.setProperty('pf4j.mode', mode)
        defaultPlugins = new DefaultPlugins()
    }

    PluginsFacade(Path root, String mode='prod') {
        this.mode = mode
        this.root = root
        System.setProperty('pf4j.mode', mode)
        defaultPlugins = new DefaultPlugins()
    }

    protected Path getPluginsDir() {
        final dir = env.get('NXF_PLUGINS_DIR')
        if( dir ) {
            log.trace "Detected NXF_PLUGINS_DIR=$dir"
            return Paths.get(dir)
        }
        else if( env.containsKey('NXF_HOME') ) {
            log.trace "Detected NXF_HOME - Using ${env.NXF_HOME}/plugins"
            return Paths.get(env.NXF_HOME, 'plugins')
        }
        else {
            log.trace "Using local plugins directory"
            return Paths.get('plugins')
        }
    }

    protected String getPluginsMode() {
        final mode = env.get('NXF_PLUGINS_MODE')
        if( mode ) {
            log.trace "Detected NXF_PLUGINS_MODE=$mode"
            return mode
        }
        else if( env.containsKey('NXF_HOME') ) {
            log.trace "Detected NXF_HOME - Using plugins mode=prod"
            return 'prod'
        }
        else {
            log.trace "Using dev plugins mode"
            return 'dev'
        }
    }

    protected boolean getPluginsDefault() {
        if( env.containsKey('NXF_PLUGINS_DEFAULT')) {
            log.trace "Detected NXF_PLUGINS_DEFAULT=$env.NXF_PLUGINS_DEFAULT"
            return env.NXF_PLUGINS_DEFAULT!='false'
        }
        else if( env.containsKey('NXF_HOME') ) {
            log.trace "Detected NXF_HOME - Using plugins defaults"
            return true
        }
        else {
            log.trace "Disabling plugins defaults"
            return false
        }
    }

    protected void init(Path root, List<PluginSpec> specs) {
        this.manager = createManager(root, specs)
        this.updater = createUpdater(root, manager)
    }

    protected Path localRoot(List<PluginSpec> specs) {
        final unique = specs ? CacheHelper.hasher(specs).hash().toString() : 'empty'
        final localRoot = PLUGINS_LOCAL_ROOT.resolve(unique)
        log.debug "Plugins local root: $localRoot"
        FilesEx.mkdirs(localRoot)
        return localRoot
    }

    protected CustomPluginManager createManager(Path root, List<PluginSpec> specs) {
        final result = mode!='dev' ? new LocalPluginManager( localRoot(specs) ) : new DevPluginManager(root)
        result.addPluginStateListener(this)
        return result
    }

    protected PluginUpdater createUpdater(Path root, CustomPluginManager manager) {
        return ( mode!='dev'
                ? new PluginUpdater(manager, root, new URL(indexUrl))
                : new DevPluginUpdater(manager) )
    }

    @Override
    void pluginStateChanged(PluginStateEvent ev) {
        final err = ev.plugin.failedException
        final dsc = ev.plugin.descriptor
        if( err ) {
            throw new IllegalStateException("Unable to start plugin id=${dsc.pluginId} version=${dsc.version} -- cause: ${err.message ?: err}", err)
        }
    }

    PluginManager getManager() { manager }

    synchronized void setup(Map config = Collections.emptyMap()) {
        if( manager )
            throw new IllegalArgumentException("Plugin system was already setup")
        else {
            log.debug "Setting up plugin manager > mode=${mode}; plugins-dir=$root"
            // make sure plugins dir exists
            if( mode!='dev' && !FilesEx.mkdirs(root) )
                throw new IOException("Unable to create plugins dir: $root")
            final specs = pluginsRequirement(config)
            init(root, specs)
            manager.loadPlugins()
            start(specs)
        }
    }

    synchronized void stop() {
        if( manager ) {
            manager.stopPlugins()
            manager = null
        }
    }

    def <T> List<T> getExtensions(Class<T> type) {
        if( manager ) {
            return manager.getExtensions(type)
        }
        else {
            // this should oly be used to load system extensions
            // i.e. included in the app class path not provided by
            // a plugin extension
            return defaultManager().getExtensions(type)
        }
    }

    @Memoized
    private PluginManager defaultManager() {
        new DefaultPluginManager()
    }

    void start( String pluginId ) {
         start( defaultPlugins.getPlugin(pluginId) )
    }

    void start(PluginSpec plugin) {
        updater.prepareAndStart(plugin.id, plugin.version)
    }

    void start(List<PluginSpec> specs) {
        for( PluginSpec it : specs ) {
            start(it)
        }
    }

    protected List<PluginSpec> pluginsRequirement(Map config) {
        def specs = parseConf(config)
        if( env.get('NXF_PACK')=='all' && specs ) {
            // custom plugins are not allowed for nextflow self-contained package
            log.warn "Nextflow self-contained distribution only allows default plugins -- User config plugins will be ignored: ${specs.join(',')}"
            return Collections.emptyList()
        }
        if( specs ) {
            log.debug "Plugins declared=$specs"
        }
        else if( getPluginsDefault() ){
            specs = defaultPluginsConf(config)
            log.debug "Plugins default=$specs"
        }

        // add tower plugin when config contains tower options
        if( config.containsKey('tower') && !specs.find {it.id == 'tower' } ) {
            specs << defaultPlugins.getPlugin('nf-tower')
        }

        return specs
    }

    protected List<PluginSpec> defaultPluginsConf(Map config) {
        // retrieve the list from the env var
        final commaSepList = env.get('NXF_PLUGINS_DEFAULT')
        if( commaSepList && commaSepList !in ['true','false'] ) {
            return commaSepList
                    .tokenize(',')
                    .collect( it-> defaultPlugins.hasPlugin(it) ? defaultPlugins.getPlugin(it) : PluginSpec.parse(it) )
        }

        // infer from app config
        final plugins = new ArrayList<PluginSpec>()
        final executor = Bolts.navigate(config, 'process.executor')

        if( executor == 'awsbatch' )
            plugins << defaultPlugins.getPlugin('nf-amazon')

        if( executor == 'google-lifesciences' )
            plugins << defaultPlugins.getPlugin('nf-google')

        if( executor == 'ignite' || System.getProperty('nxf.node.daemon')=='true') {
            plugins << defaultPlugins.getPlugin('nf-ignite')
            plugins << defaultPlugins.getPlugin('nf-amazon')
        }

        if( !plugins ) {
            // always include amazon plugin for backward compatibility
            plugins << defaultPlugins.getPlugin('nf-amazon')
        }

        return plugins
    }

    /**
     * Lookup for plugins declared in the nextflow.config using the `plugins` scope
     *
     * @param config The nextflow config as a Map object
     * @return The list of declared plugins
     */
    protected List<PluginSpec> parseConf(Map config) {
        final pluginsConf = config.plugins as List<String>
        final result = new ArrayList( pluginsConf?.size() ?: 0 )
        if(pluginsConf) for( String it : pluginsConf ) {
            result.add( PluginSpec.parse(it) )
        }
        return result
    }

    synchronized void pullPlugins(List<String> ids) {
        updater.pullPlugins(ids)
    }

}
