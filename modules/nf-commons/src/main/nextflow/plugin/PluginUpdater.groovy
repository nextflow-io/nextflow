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

import static java.nio.file.StandardCopyOption.*

import java.nio.file.Files
import java.nio.file.Path
import java.util.function.Predicate
import java.util.regex.Pattern

import com.github.zafarkhaja.semver.Version
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import nextflow.Const
import nextflow.SysEnv
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.file.FileMutex
import org.pf4j.PluginDependency
import org.pf4j.PluginRuntimeException
import org.pf4j.PluginState
import org.pf4j.PluginWrapper
import org.pf4j.RuntimeMode
import org.pf4j.update.DefaultUpdateRepository
import org.pf4j.update.PluginInfo
import org.pf4j.update.UpdateManager
import org.pf4j.update.UpdateRepository
import org.pf4j.util.FileUtils
/**
 * Implements the download/install/update for nextflow plugins
 *
 * Plugins are downloaded and stored uncompressed in the
 * directory defined by the variable `NXF_PLUGINS_DIR`, $NXF_HOME/nextflow
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PluginUpdater extends UpdateManager {

    static final public Pattern META_REGEX = ~/(.+)-(\d+\.\d+\.\d+\S*)-meta\.json/

    private CustomPluginManager pluginManager

    private Path pluginsStore

    private boolean pullOnly

    private DefaultPlugins defaultPlugins = DefaultPlugins.INSTANCE

    protected PluginUpdater(CustomPluginManager pluginManager) {
        super(pluginManager)
        this.pluginManager = pluginManager
    }

    PluginUpdater(CustomPluginManager pluginManager, Path pluginsRoot, URL repo) {
        super(pluginManager, wrap(repo))
        this.pluginsStore = pluginsRoot
        this.pluginManager = pluginManager
    }

    static private List<UpdateRepository> wrap(URL repo) {
        List<UpdateRepository> result = new ArrayList<>(1)
        result << new DefaultUpdateRepository('nextflow.io', repo)
        result.addAll(customRepos())
        return result
    }

    static private List<DefaultUpdateRepository> customRepos() {
        final repos = SysEnv.get('NXF_PLUGINS_TEST_REPOSITORY')
        if( !repos )
            return List.of()
        // warn the user that a custom
        final msg = """\
                        =======================================================================
                        =                                WARNING                                    =
                        = You are running this script using a un-official plugin repository.        =
                        =                                                                           =
                        = ${repos}
                        =                                                                           =
                        = This is only meant to be used for plugin testing purposes.                =
                        =============================================================================
                        """.stripIndent(true)
        log.warn(msg)
        final result = new ArrayList<DefaultUpdateRepository>(10)
        // the repos string can contain one or more plugin repository uri separated by comma
        for( String it : repos.tokenize(',') )
            result.add(customRepo(it))
        return result
    }

    static private DefaultUpdateRepository customRepo(String uri) {
        // Check if it's a plugin meta file. The name must match the pattern `<plugin id>-X.Y.Z-meta.json`
        final matcher = META_REGEX.matcher(uri.tokenize('/')[-1])
        if( matcher.matches() ) {
            try {
                final pluginId = matcher.group(1)
                final temp = File.createTempFile('nxf-','json')
                temp.deleteOnExit()
                temp.text = /[{"id":"${pluginId}", "releases":[ ${new URL(uri).text} ]}]/
                uri = 'file://' + temp.absolutePath
            }
            catch (FileNotFoundException e) {
                throw new IllegalArgumentException("Provided repository URL does not exist or cannot be accessed: $uri")
            }
        }
        // create the update repository instance
        final fileName = uri.tokenize('/')[-1]
        return new DefaultUpdateRepository('uri', new URL(uri), fileName)
    }

    /**
     * Resolve a plugin installing or updating the dependencies if necessary
     * and start the plugin
     *
     * @param pluginId The plugin Id
     * @param version The required version, can be null for the latest version
     */
    void prepareAndStart(String pluginId, String version) {
        final PluginWrapper current = pluginManager.getPlugin(pluginId)
        if( !current ) {
            log.debug "Installing plugin ${pluginId} version: ${version ?: 'latest'}"
            // install & start the plugin
            installPlugin(pluginId, version)
        }
        else if( shouldUpdate(pluginId, version, current) ) {
            log.debug "Updating plugin ${pluginId} version: ${version ?: 'latest'} [current version: $current.descriptor.version]"
            // update & start the plugin
            updatePlugin(pluginId, version)
        }
        else {
            if( !version ) version = current.descriptor.version
            log.debug "Starting plugin ${pluginId} version: ${version}"
            pluginManager.startPlugin(pluginId)
        }
    }

    void pullPlugins(List<String> plugins) {
        pullOnly=true
        try {
            final specs = plugins.collect(it -> PluginSpec.parse(it,defaultPlugins))
            for( PluginSpec spec : specs ) {
                pullPlugin0(spec.id, spec.version)
            }
        }
        finally {
            pullOnly=false
        }
    }

    void pullPlugin0(String pluginId, String version) {
        final PluginWrapper current = pluginManager.getPlugin(pluginId)
        if( !current ) {
            log.debug "Installing plugin ${pluginId} version: ${version ?: 'latest'}"
            // install & start the plugin
            installPlugin(pluginId, version)
        }
        else if( shouldUpdate(pluginId, version, current) ) {
            log.debug "Updating plugin ${pluginId} version: ${version ?: 'latest'} [current version: $current.descriptor.version]"
            // update & start the plugin
            updatePlugin(pluginId, version)
        }
    }

    /**
     * Install a new plugin downloading the artifact from the remote source if needed
     *
     * @param id The plugin Id
     * @param version The plugin version
     * @return {@code true} when the plugin is correctly download, installed and started, {@code false} otherwise
     */
    @Override
    boolean installPlugin(String id, String version) {
        if( !pluginsStore )
            throw new IllegalStateException("Missing pluginStore attribute")

        return load0(id, version)
    }

    private Path download0(String id, String version) {

        // 0. check if already exists
        final pluginPath = pluginsStore.resolve("$id-$version")
        if( FilesEx.exists(pluginPath) ) {
            return pluginPath
        }

        // 1. determine the version
        if( !version )
            version = getLastPluginRelease(id)
        log.info "Downloading plugin ${id}@${version}"

        // 2. Download to temporary location
        Path downloaded = safeDownloadPlugin(id, version);

        // 3. unzip the content and delete downloaded file
        Path dir = FileUtils.expandIfZip(downloaded)
        FileHelper.deletePath(downloaded)

        // 4. move the final destination the plugin directory
        assert pluginPath.getFileName() == dir.getFileName()

        try {
            safeMove(dir, pluginPath)
        }
        catch (IOException e) {
            throw new PluginRuntimeException(e, "Failed to write file '$pluginPath' to plugins folder")
        }

        return pluginPath
    }

    protected Path safeDownloadPlugin(String id, String version) {
        final CheckedSupplier<Path> supplier = () -> downloadPlugin(id, version)
        final policy = retryPolicy(id,version)
        return Failsafe.with(policy).get(supplier)
    }

    protected <T> RetryPolicy<T> retryPolicy(String id, String version) {
        final listener = new dev.failsafe.event.EventListener<ExecutionAttemptedEvent<T>>() {
            @Override
            void accept(ExecutionAttemptedEvent<T> event) throws Throwable {
                log.debug("Failed to download plugin: $id; version: $version - attempt: ${event.attemptCount}", event.lastFailure)
            }
        }

        final condition = new Predicate<Throwable>() {
            @Override
            boolean test(Throwable error) {
                return error?.cause instanceof ConnectException
            }
        }

        return RetryPolicy.<T>builder()
                .handleIf(condition)
                .withMaxAttempts(3)
                .onRetry(listener)
                .build()
    }

    protected void safeMove(Path source, Path target) {
        try {
            Files.move(source, target, ATOMIC_MOVE, REPLACE_EXISTING)
        }
        catch (IOException e) {
            log.debug "Failed atomic move for plugin $source -> $target - Reason: ${e.message ?: e} - Fallback on safe move"
            safeMove0(source, target)
        }
    }

    protected void safeMove0(Path source, Path target) {
        // make sure to the target path does not exist
        FileHelper.deletePath(target)
        // copy the source to the more
        FileHelper.copyPath(source, target)
        // finally remove the source
        try {
            FileHelper.deletePath(source)
        }
        catch (IOException e) {
            log.warn("Unable to delete plugin directory: $source", e)
        }
    }

    /**
     * Race condition safe plugin download. Multiple instances are synchronised
     * using a file system lock created in the tmp directory
     *
     * @param id The plugin Id
     * @param version The plugin version string
     * @return The uncompressed plugin directory path
     */
    private Path safeDownload(String id, String version) {
        final sentinel = lockFile(id,version)
        final mutex = new FileMutex(target: sentinel, timeout: '10 min')
        try {
            return mutex.lock { download0(id, version) }
        }
        finally {
            sentinel.delete()
        }
    }

    /**
     * Get the file to be used a lock to synchronise concurrent downloaded from multiple Nextflow launches
     *
     * @param id The plugin Id
     * @param version The plugin version
     * @return The lock the path
     */
    private File lockFile(String id, String version) {
        def tmp = System.getProperty('java.io.tmpdir')
        new File(tmp, "nextflow-plugin-${id}-${version}.lock")
    }

    private boolean load0(String id, String version) {
        assert id, "Missing plugin Id"

        if( version == null )
            version = getLastPluginRelease(id)?.version

        final offline = SysEnv.get('NXF_OFFLINE')=='true'
        if( !version ) {
            final msg = offline
                ? "Cannot find version for $id plugin -- plugin versions MUST be specified in offline mode"
                : "Cannot find latest version of $id plugin"
            throw new IllegalStateException(msg)
        }

        def pluginPath = pluginsStore.resolve("$id-$version")
        if( !FilesEx.exists(pluginPath) ) {
            pluginPath = safeDownload(id, version)
        }

        // verify the plugin install path contains the expected manifest path
        if( !FilesEx.exists(pluginPath.resolve('classes/META-INF/MANIFEST.MF')) ) {
            log.warn("Plugin '${pluginPath.getFileName()}' installation looks corrupted - Delete the following directory and run nextflow again: $pluginPath")
        }

        // load the plugin from the file system
        PluginWrapper wrapper = pluginManager.loadPluginFromPath(pluginPath)

        // pull all required required deps
        final deps = wrapper.descriptor.dependencies ?: Collections.<PluginDependency>emptyList()
        for( PluginDependency it : deps ) {
            // 1. check for installed version satisfying req
            final installed = checkInstalled(it.pluginId, it.pluginVersionSupport)
            if( installed ) {
                log.debug "Plugin $id requires $it.pluginId supported version: $it.pluginVersionSupport - found version: $installed"
                continue
            }

            // 2. find latest satisfying req
            // -- if it's a core nextflow plugin use the version expected by it
            def depVersion = defaultPlugins.getPlugin(it.pluginId)?.version
            // -- otherwise try to find the newest matching release
            if( !depVersion )
                depVersion = findNewestMatchingRelease(it.pluginId, it.pluginVersionSupport)?.version
            log.debug "Plugin $id requires $it.pluginId supported version: $it.pluginVersionSupport - available version: $depVersion"
            if( pullOnly )
                pullPlugin0(it.pluginId, depVersion)
            else
                prepareAndStart(it.pluginId, depVersion)
        }

        if( pullOnly )
            return false

        // resolve the plugins
        pluginManager.resolvePlugins()
        // finally start it
        PluginState state = pluginManager.startPlugin(id)
        return PluginState.STARTED == state
    }

    /**
     * Update an existing plugin, unloading the current version and installing
     * the requested one downloading the artifact if needed
     *
     * @param id The plugin Id
     * @param version The plugin version
     * @return {@code true} when the plugin is updated and started or {@code false} otherwise
     */
    @Override
    boolean updatePlugin(String id, String version) {
        if (pluginManager.getPlugin(id) == null) {
            throw new PluginRuntimeException("Plugin $id cannot be updated since it is not installed")
        }

        PluginInfo pluginInfo = getPluginsMap().get(id)
        if (pluginInfo == null) {
            throw new PluginRuntimeException("Plugin $id does not exist in any repository")
        }

        if (!pluginManager.deletePlugin(id)) {
            return false
        }

        load0(id, version)
    }

    protected boolean shouldUpdate(String pluginId, String version, PluginWrapper current) {
        if( pluginManager.runtimeMode == RuntimeMode.DEVELOPMENT ) {
            log.debug "Update not supported during development mode"
            return false
        }

        if( !version )
            version = getLastPluginRelease(pluginId)?.version
        if( !version ) {
            log.warn "Cannot find latest version for plugin $pluginId [skip update]"
            return false
        }

        return pluginManager
                .getVersionManager()
                .compareVersions(current.descriptor.version, version) < 0
    }

    protected String checkInstalled(String id, String verConstraint) {
        final versionManager = pluginManager.getVersionManager()
        // check if the installed version satisfies the requirement
        final current = pluginManager.getPlugin(id)
        if( !current )
            return null

        final found = versionManager.checkVersionConstraint(current.descriptor.version, verConstraint)
        return found ? current.descriptor.version : null
    }

    /**
     * Find newest release version matching the requested plugin constraint and nextflow runtime
     *
     * @param id The plugin id
     * @param verConstraint The version constraint
     * @return The plugin release satisfying the requested criteria or null it none is found
     */
    protected PluginInfo.PluginRelease findNewestMatchingRelease(String id, String verConstraint) {
        assert verConstraint

        final versionManager = pluginManager.getVersionManager()

        PluginInfo pluginInfo = getPluginsMap().get(id)
        if( !pluginInfo )
            throw new IllegalArgumentException("Unknown plugin id: $id")

        // note: order releases list by descending version numbers ie. latest version comes first
        def releases = pluginInfo.releases.sort(false) { a,b -> Version.valueOf(b.version) <=> Version.valueOf(a.version) }
        for (PluginInfo.PluginRelease rel : releases ) {
            if( !versionManager.checkVersionConstraint(rel.version, verConstraint) || !rel.url )
                continue

            if( versionManager.checkVersionConstraint(BuildInfo.version, rel.requires) )
                return rel
        }

        return null
    }
}
