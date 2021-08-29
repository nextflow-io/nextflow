/*
 * Copyright 2020-2021, Seqera Labs
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
import org.pf4j.PluginState
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
    private static final String DEV_MODE = 'dev'
    private static final String PROD_MODE = 'prod'
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
        if( mode=='dev' && root.toString()=='plugins' )
            root = detectDevPluginsRoot()
        System.setProperty('pf4j.mode', mode)
        defaultPlugins = new DefaultPlugins()
    }

    PluginsFacade(Path root, String mode=PROD_MODE) {
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

    protected Path detectDevPluginsRoot() {
        def file = new File('.').absoluteFile
        do {
            if( file.name=='nextflow' && file.isDirectory() && new File(file, 'settings.gradle').isFile() ) {
                final root = file.toPath().resolve('plugins')
                log.debug "Detected dev plugins root: $root"
                return root
            }
            file = file.parentFile
        }
        while( file!=null )
        throw new IllegalStateException("Unable to detect local plugins root")
    }

    protected String getPluginsMode() {
        final mode = env.get('NXF_PLUGINS_MODE')
        if( mode ) {
            log.trace "Detected NXF_PLUGINS_MODE=$mode"
            return mode
        }
        else if( env.containsKey('NXF_HOME') ) {
            log.trace "Detected NXF_HOME - Using plugins mode=prod"
            return PROD_MODE
        }
        else {
            log.trace "Using dev plugins mode"
            return DEV_MODE
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
        final result = mode!=DEV_MODE ? new LocalPluginManager( localRoot(specs) ) : new DevPluginManager(root)
        result.addPluginStateListener(this)
        return result
    }

    protected PluginUpdater createUpdater(Path root, CustomPluginManager manager) {
        return ( mode!=DEV_MODE
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
            if( mode!=DEV_MODE && !FilesEx.mkdirs(root) )
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

    /**
     * Return a list of extension matching the requested interface type
     *
     * @param type
     *      The request extension interface
     * @return
     *      The list of extensions matching the requested interface.
     */
    def <T> List<T> getExtensions(Class<T> type) {
        if( manager ) {
            return manager.getExtensions(type)
        }
        else {
            // this should oly be used to load system extensions
            // i.e. included in the app class path not provided by
            // a plugin extension
            log.debug "Using Default plugins manager"
            return defaultManager().getExtensions(type)
        }
    }

    /**
     * Return a list of extension matching the requested type
     * ordered by a priority value. The element at the beginning
     * of the list (index 0) has higher priority
     *
     * @param type
     *      The request extension interface
     * @return
     *      The list of extensions matching the requested interface.
     *      The extension with higher priority appears first (lower index)
     */
    def <T> List<T> getPriorityExtensions(Class<T> type,String group=null) {
        def result = getExtensions(type)
        if( group )
            result = result.findAll(it -> group0(it)==group )
        return result.sort( it -> priority0(it) )
    }

    protected int priority0(Object it) {
        final annot = it.getClass().getAnnotation(Priority)
        return annot ? annot.value() : 0
    }

    protected int priority1(Object it) {
        final annot = it.getClass().getAnnotation(Scoped)
        return annot ? annot.priority() : 0
    }

    /**
     * Find out all extensions classed with with {@code @Scoped} annotation
     * return one and exactly with for each different scope value
     *
     * @param type
     * @param scope
     * @return
     */
    def <T> Set<T> getScopedExtensions(Class<T> type,String scope=null) {
        def result = getExtensions(type).sort(it->priority1(it))
        def groups = new HashMap<String,T>()
        for( T it : result ) {
            final annot = it.getClass().getAnnotation(Scoped)
            if( annot==null )
                continue
            if( !annot.value() )
                continue
            if( groups.containsKey(annot.value()) )
                continue
            if( scope==null || annot.value()==scope )
                groups.put(annot.value(), it)
        }
        return new HashSet<T>(groups.values())
    }

    protected String group0(Object it) {
        final annot = it.getClass().getAnnotation(Priority)
        return annot && annot.group() ? annot.group() : null
    }

    @Memoized
    private PluginManager defaultManager() {
        new DefaultPluginManager()
    }

    void start( String pluginId ) {
        if( isSelfContained() ) {
            log.debug "Plugin 'start' is not required in self-contained mode -- ignoring it for plugin: $pluginId"
            return
        }

        start( defaultPlugins.getPlugin(pluginId) )
    }

    void start(PluginSpec plugin) {
        if( isSelfContained() ) {
            log.debug "Plugin 'start' is not required in self-contained mode -- ignoring it for plugin: $plugin.id"
            return
        }

        updater.prepareAndStart(plugin.id, plugin.version)
    }

    void start(List<PluginSpec> specs) {
        if( isSelfContained() ) {
            log.debug "Plugin 'start' is not required in self-contained mode -- ignoring it for plugins: $specs"
            return
        }

        for( PluginSpec it : specs ) {
            start(it)
        }
    }

    boolean isStarted(String pluginId) {
        manager.getPlugin(pluginId)?.pluginState == PluginState.STARTED
    }

    /**
     * @return {@code true} when running in self-contained mode ie. the nextflow distribution
     * include also plugin libraries. When running is this mode, plugins should not be started
     * and cannot be updated. 
     */
    protected boolean isSelfContained() {
        return env.get('NXF_PACK')=='all'
    }

    protected List<PluginSpec> pluginsRequirement(Map config) {
        def specs = parseConf(config)
        if( isSelfContained() && specs ) {
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
        if( (config.containsKey('tower') || env.TOWER_ACCESS_TOKEN ) && !specs.find {it.id == 'tower' } ) {
            specs << defaultPlugins.getPlugin('nf-tower')
        }

        return specs
    }

    protected List<PluginSpec> defaultPluginsConf(Map config) {
        // retrieve the list from the env var
        final commaSepList = env.get('NXF_PLUGINS_DEFAULT')
        if( commaSepList && commaSepList !in ['true','false'] ) {
            // if the plugin id in the list does *not* contain the @version suffix, it picks the version
            // specified in the defaults list. Otherwise parse the provider id@version string to the corresponding spec
            return commaSepList
                    .tokenize(',')
                    .collect( it-> defaultPlugins.hasPlugin(it) ? defaultPlugins.getPlugin(it) : PluginSpec.parse(it) )
        }

        // infer from app config
        final plugins = new ArrayList<PluginSpec>()
        final workDir = config.workDir as String
        final executor = Bolts.navigate(config, 'process.executor')

        if( executor == 'awsbatch' || workDir?.startsWith('s3://') )
            plugins << defaultPlugins.getPlugin('nf-amazon')

        if( executor == 'google-lifesciences' || workDir?.startsWith('gs://') )
            plugins << defaultPlugins.getPlugin('nf-google')

        if( executor == 'azurebatch' || workDir?.startsWith('az://') )
            plugins << defaultPlugins.getPlugin('nf-azure')

        if( executor == 'ignite' || System.getProperty('nxf.node.daemon')=='true') {
            plugins << defaultPlugins.getPlugin('nf-ignite')
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
            result.add( PluginSpec.parse(it, defaultPlugins) )
        }
        return result
    }

    synchronized void pullPlugins(List<String> ids) {
        updater.pullPlugins(ids)
    }

    boolean startIfMissing(String pluginId) {
        if( env.NXF_PLUGINS_DEFAULT == 'false' )
            return false

        if( isStarted(pluginId) )
            return false

        synchronized (this) {
            if( isStarted(pluginId) )
                return false
            start(pluginId)
            return true
        }
    }
}
