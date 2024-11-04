/*
 * Copyright 2013-2024, Seqera Labs
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
import nextflow.SysEnv
import nextflow.extension.Bolts
import nextflow.extension.FilesEx
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

    private static final String DEV_MODE = 'dev'
    private static final String PROD_MODE = 'prod'
    private Map<String,String> env = SysEnv.get()

    private String mode
    private Path root
    private PluginUpdater updater
    private CustomPluginManager manager
    private DefaultPlugins defaultPlugins = DefaultPlugins.INSTANCE
    private String indexUrl = Plugins.DEFAULT_PLUGINS_REPO
    private boolean embedded

    PluginsFacade() {
        mode = getPluginsMode()
        root = getPluginsDir()
        if( mode==DEV_MODE && root.toString()=='plugins' && !isRunningFromDistArchive() )
            root = detectPluginsDevRoot()
        System.setProperty('pf4j.mode', mode)
    }

    PluginsFacade(Path root, String mode=PROD_MODE) {
        this.mode = mode
        this.root = root
        System.setProperty('pf4j.mode', mode)
    }

    /**
     * Determine if it's running from a JAR archive
     * @return {@code true} if the code is running from a JAR artifact, {@code false} otherwise
     */
    protected String isRunningFromDistArchive() {
        final className = this.class.name.replace('.', '/');
        final classJar = this.class.getResource("/" + className + ".class").toString();
        return classJar.startsWith("jar:")
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

    private boolean isNextflowDevRoot(File file) {
        file.name=='nextflow' && file.isDirectory() && new File(file, 'settings.gradle').isFile()
    }

    private Path pluginsDevRoot(File path) {
        // main project root
        if( isNextflowDevRoot(path) )
            return path.toPath().resolve('plugins')
        // when nextflow is included into another build, check the sibling directory
        if( new File(path,'settings.gradle').exists() && isNextflowDevRoot(path=new File(path,'../nextflow')) )
            return path.toPath().resolve('plugins')
        else
            return null
    }

    /**
     * Determine the development plugin root. This is required to
     * allow running unit tests for plugin projects importing the
     * nextflow core runtime.
     *
     * The nextflow main project is expected to be cloned into
     * a sibling directory respect to the plugin project
     *
     * @return The nextflow plugins project path in the local file system
     */
    protected Path detectPluginsDevRoot() {
        def file = new File('.').absoluteFile
        do {
            final root = pluginsDevRoot(file)
            if( root ) {
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
            log.debug "Using dev plugins mode"
            return DEV_MODE
        }
    }

    protected boolean getPluginsDefault() {
        if( env.containsKey('NXF_PLUGINS_DEFAULT')) {
            log.trace "Detected NXF_PLUGINS_DEFAULT=$env.NXF_PLUGINS_DEFAULT"
            return env.NXF_PLUGINS_DEFAULT!='false'
        }
        else if( env.containsKey('NXF_HOME') ) {
            log.trace "Detected NXF_HOME - Using plugin defaults"
            return true
        }
        else {
            log.trace "Disabling plugin defaults"
            return false
        }
    }

    private CustomPluginManager newPluginManager(Path root, boolean embedded) {
        if( mode==DEV_MODE ) {
            // plugin manage for dev purposes
            return new DevPluginManager(root)
        }
        if( embedded ) {
            // use the custom plugin manager to by-pass the creation of a local plugin repository
            return new EmbeddedPluginManager(root)
        }
        return new LocalPluginManager(root)
    }

    protected CustomPluginManager createManager(Path root, boolean embedded) {
        final result = newPluginManager(root, embedded)
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

    void init(boolean embedded=false) {
        if( manager )
            throw new IllegalArgumentException("Plugin system already setup")

        log.debug "Setting up plugin manager > mode=${mode}; embedded=$embedded; plugins-dir=$root; core-plugins: ${defaultPlugins.toSortedString()}"
        // make sure plugins dir exists
        if( mode!=DEV_MODE && !FilesEx.mkdirs(root) )
            throw new IOException("Unable to create plugins dir: $root")

        this.manager = createManager(root, embedded)
        this.updater = createUpdater(root, manager)
        manager.loadPlugins()
        if( embedded ) {
            manager.startPlugins()
            this.embedded = embedded
        }
    }

    void init(Path root, String mode, CustomPluginManager pluginManager) {
        if( manager )
            throw new IllegalArgumentException("Plugin system already setup")
        this.root = root
        this.mode = mode
        // setup plugin manager
        this.manager = pluginManager
        this.manager.addPluginStateListener(this)
        // setup the updater
        this.updater = createUpdater(root, manager)
        // load plugins
        manager.loadPlugins()
        if( embedded ) {
            manager.startPlugins()
            this.embedded = embedded
        }
    }

    void load(Map config) {
        if( !manager )
            throw new IllegalArgumentException("Plugin system has not been initialized")
        start(pluginsRequirement(config))
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
            // this should only be used to load system extensions
            // i.e. included in the app class path not provided by
            // a plugin extension
            log.debug "Using Default plugin manager"
            return defaultManager().getExtensions(type)
        }
    }

    /**
     * Return a list of extension matching the requested interface type in a plugin
     *
     * @param type
     *      The request extension interface
     * @return
     *      The list of extensions matching the requested interface.
     */
    def <T> List<T> getExtensions(Class<T> type, String pluginId) {
        if( manager ) {
            return manager.getExtensions(type, pluginId)
        }
        else {
            return List.of()
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

    protected String group0(Object it) {
        final annot = it.getClass().getAnnotation(Priority)
        return annot && annot.group() ? annot.group() : null
    }

    @Memoized
    private PluginManager defaultManager() {
        new DefaultPluginManager()
    }

    void start( String pluginId ) {
        if( isEmbedded() && defaultPlugins.hasPlugin(pluginId) ) {
            log.debug "Plugin 'start' is not required in embedded mode -- ignoring for plugin: $pluginId"
            return
        }

        start(PluginSpec.parse(pluginId, defaultPlugins))
    }

    void start(PluginSpec plugin) {
        if( isEmbedded() && defaultPlugins.hasPlugin(plugin.id) ) {
            log.debug "Plugin 'start' is not required in embedded mode -- ignoring for plugin: $plugin.id"
            return
        }

        updater.prepareAndStart(plugin.id, plugin.version)
    }

    void start(List<PluginSpec> specs) {
        for( PluginSpec it : specs ) {
            start(it)
        }
    }

    boolean isStarted(String pluginId) {
        manager.getPlugin(pluginId)?.pluginState == PluginState.STARTED
    }

    /**
     * @return {@code true} when running in embedded mode ie. the nextflow distribution
     * include also plugin libraries. When running is this mode, plugins should not be started
     * and cannot be updated. 
     */
    protected boolean isEmbedded() {
        return embedded
    }

    protected List<PluginSpec> pluginsRequirement(Map config) {
        def specs = parseConf(config)
        if( isEmbedded() && specs ) {
            // custom plugins are not allowed for nextflow self-contained package
            log.warn "Nextflow embedded mode only core plugins -- User config plugins will be ignored: ${specs.join(',')}"
            return Collections.emptyList()
        }
        if( specs ) {
            log.debug "Plugins declared=$specs"
        }
        if( getPluginsDefault() ){
            final defSpecs = defaultPluginsConf(config)
            specs = mergePluginSpecs(specs, defSpecs)
            log.debug "Plugins default=$defSpecs"
        }

        // add tower plugin when config contains tower options
        if( (Bolts.navigate(config,'tower.enabled') || env.TOWER_ACCESS_TOKEN ) && !specs.find {it.id == 'nf-tower' } ) {
            specs << defaultPlugins.getPlugin('nf-tower')
        }
        if( (Bolts.navigate(config,'wave.enabled') || Bolts.navigate(config,'fusion.enabled')) && !specs.find {it.id == 'nf-wave' } ) {
            specs << defaultPlugins.getPlugin('nf-wave')
        }

        // add cloudcache plugin when cloudcache is enabled in the config
        if( Bolts.navigate(config, 'cloudcache.enabled')==true ) {
            specs << defaultPlugins.getPlugin('nf-cloudcache')
        }

        log.debug "Plugins resolved requirement=$specs"
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
        final bucketDir = config.bucketDir as String
        final executor = Bolts.navigate(config, 'process.executor')

        if( executor == 'awsbatch' || workDir?.startsWith('s3://') || bucketDir?.startsWith('s3://') || env.containsKey('NXF_ENABLE_AWS_SES') )
            plugins << defaultPlugins.getPlugin('nf-amazon')

        if( executor == 'google-lifesciences' || executor == 'google-batch' || workDir?.startsWith('gs://') || bucketDir?.startsWith('gs://')  )
            plugins << defaultPlugins.getPlugin('nf-google')

        if( executor == 'azurebatch' || workDir?.startsWith('az://') || bucketDir?.startsWith('az://') )
            plugins << defaultPlugins.getPlugin('nf-azure')

        if( Bolts.navigate(config, 'weblog.enabled'))
            plugins << new PluginSpec('nf-weblog')
            
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
        if( isEmbedded() && defaultPlugins.hasPlugin(pluginId) )
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

    /**
     * Merge two lists of plugin requirements
     *
     * @param configPlugins
     *      The list of plugins specified via the configuration file. This has higher priority
     * @param defaultPlugins
     *      The list of plugins specified via the environment
     * @return
     *      The list of plugins resulting from merging the two lists
     */
    protected List<PluginSpec> mergePluginSpecs(List<PluginSpec> configPlugins, List<PluginSpec> defaultPlugins) {
        final map = new LinkedHashMap<String,PluginSpec>(10)
        // add all plugins in the 'configPlugins' argument
        for( PluginSpec plugin : configPlugins ) {
            map.put(plugin.id, plugin)
        }
        // add the plugin in the 'defaultPlugins' argument
        // if the map already contains the plugin,
        // override it only if it does not specify a version
        for( PluginSpec plugin : defaultPlugins ) {
            if( !map[plugin.id] || !map[plugin.id].version ) {
                map.put(plugin.id, plugin)
            }
        }
        return new ArrayList<PluginSpec>(map.values())
    }
}
