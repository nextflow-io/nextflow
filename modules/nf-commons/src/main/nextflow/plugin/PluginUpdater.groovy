/*
 * Copyright 2013-2026, Seqera Labs
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
 */

package nextflow.plugin

import static java.nio.file.StandardCopyOption.*

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap
import java.util.regex.Pattern

import com.github.zafarkhaja.semver.Version
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedPredicate
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import nextflow.SysEnv
import nextflow.config.RegistryConfig
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.file.FileMutex
import org.pf4j.InvalidPluginDescriptorException
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

    private boolean offline

    private boolean registryReposApplied

    private DefaultPlugins defaultPlugins = DefaultPlugins.INSTANCE

    // plain name of a moved plugin to its namespaced ID
    @PackageScope
    final Map<String,String> aliases = new ConcurrentHashMap<>()

    protected PluginUpdater(CustomPluginManager pluginManager) {
        super(pluginManager)
        this.pluginManager = pluginManager
    }

    PluginUpdater(CustomPluginManager pluginManager, Path pluginsRoot, URL repo, boolean offline) {
        super(pluginManager, wrap(repo, pluginsRoot, offline))
        this.offline = offline
        this.pluginsStore = pluginsRoot
        this.pluginManager = pluginManager
        if( pluginsRoot && Files.isDirectory(pluginsRoot) )
            aliases.putAll(storeAliases(pluginsRoot))
    }

    /**
     * The namespaced ID of a moved plugin, or the given ID
     */
    String canonical(String id) {
        return id ? aliases.getOrDefault(id, id) : id
    }

    private void readCanonicalIds() {
        // the first registry wins
        final result = new HashMap<String,String>()
        for( def repo : this.@repositories ) {
            if( repo instanceof HttpPluginRepository )
                repo.canonicalIds.each { k, v -> result.putIfAbsent(k, v) }
        }
        aliases.putAll(result)
    }

    // trust only links: only Nextflow writes them
    static Map<String,String> storeAliases(Path store) {
        final List<String> ids = FilesEx.listFiles(store)
            .findAll { Path p -> Files.isSymbolicLink(p) && Files.exists(p) }
            .collect { Path p -> [p.fileName.toString(), Files.readSymbolicLink(p).fileName.toString()].findAll { String n -> n.contains('~') } }
            .findAll { List<String> names -> names.size() == 1 }
            .collect { List<String> names -> names[0].replaceFirst(/-\d+\.\d+\.\d+\S*$/, '').replace('~', '/') }
            .unique()
        // a name linked to two namespaces is ambiguous
        return ids
            .groupBy { String id -> PluginRef.name(id) }
            .findAll { String k, List<String> v -> v.size() == 1 }
            .collectEntries { String k, List<String> v -> [k, v[0]] } as Map<String,String>
    }

    private boolean isMoved(String id) {
        return id.contains('/') && aliases.get(PluginRef.name(id)) == id
    }

    // load0 needs a real folder, not a link
    private static boolean isComplete(Path dir) {
        return Files.isDirectory(dir, LinkOption.NOFOLLOW_LINKS) && Files.exists(dir.resolve('classes/META-INF/MANIFEST.MF'))
    }

    private static boolean createLink(Path link, Path target) {
        try {
            Files.createSymbolicLink(link, target.fileName)
            return true
        }
        catch( IOException | UnsupportedOperationException e ) {
            log.debug "Unable to link plugin $link -> $target - cause: ${e.message}"
            return false
        }
    }

    // takes the plain lock, which older Nextflow holds while it copies the folder
    private boolean linkToPlain(Path link, Path plainPath, String id, String version) {
        return withLock(PluginRef.name(id), version) { isComplete(plainPath) && createLink(link, plainPath) }
    }

    static private List<UpdateRepository> wrap(URL remote, Path local, boolean offline) {
        List<UpdateRepository> result = new ArrayList<>(1)
        if( offline ) {
            log.debug "Using local update repository: ${local}"
            result.add(new LocalUpdateRepository('downloaded', local))
        }
        else {
            def remoteRepo = remote.path.endsWith('.json')
                ? new DefaultUpdateRepository('nextflow.io', remote)
                : new HttpPluginRepository('registry', remote.toURI())

            log.debug "Using plugin repository: ${remoteRepo.getClass().getSimpleName()} [${remoteRepo.id}]; url=${remote}"
            result.add(remoteRepo)
            result.addAll(customRepos())
        }
        return result
    }

    /** Ids of the default registry repositories created by {@link #wrap} */
    private static final List<String> DEFAULT_REGISTRY_REPO_IDS = ['registry', 'nextflow.io']

    /**
     * Apply the plugin registry endpoints declared in the {@code registry} config scope.
     *
     * The configured registries are authoritative: the URLs provided by {@link RegistryConfig}
     * fully replace the default registry repository. Users that want to keep resolving plugins
     * from the public registry must include its URL explicitly in {@code registry.url}.
     *
     * Callers must only invoke this when {@code registry.url} is explicitly configured (see
     * {@link PluginsFacade#applyRegistryConfig}); when it is empty or unset the default registry
     * is left in place. {@link RegistryConfig#getAllUrls} always falls back to the default URL,
     * so it is never empty here.
     *
     * Safe to call after construction and before {@link #prefetchMetadata}, which initialises
     * each {@link HttpPluginRepository} with the metadata it actually needs.
     *
     * The registries are applied only once: repeated invocations are a no-op, so re-entrant
     * callers (e.g. a second {@link PluginsFacade#load}) cannot mint duplicate repositories.
     */
    void addRegistryRepos(RegistryConfig registryConfig) {
        if( offline || !registryConfig || registryReposApplied )
            return
        registryReposApplied = true
        final urls = registryConfig.getAllUrls()
        // the configured registries take over: drop the default registry repository(ies).
        // Mutate the list directly rather than removeRepository(id): pf4j-update 2.3.0 logs a
        // misleading "Repository with id X not found" warning even on a successful removal.
        this.@repositories.removeIf { it.id in DEFAULT_REGISTRY_REPO_IDS }
        final existingIds = this.@repositories*.id as Set
        int counter = 0
        for( String url : urls ) {
            String repoId = "registry-${counter++}"
            while( repoId in existingIds )
                repoId = "registry-${counter++}"
            existingIds.add(repoId)
            final repo = new HttpPluginRepository(repoId, URI.create(url))
            log.debug "Adding plugin repository: ${repo.getClass().getSimpleName()} [${repo.id}]; url=${url}"
            addRepository(repo)
        }
    }

    /**
     * Aggregate plugin metadata across all configured registries.
     *
     * The default {@link UpdateManager#getPluginsMap()} merges repositories with
     * {@code Map.putAll}, so the last repository serving a given plugin id fully overwrites
     * the earlier ones — silently discarding releases that exist only in an earlier-listed
     * registry. Here the releases of each plugin id are instead unioned across all
     * repositories. When the same {@code version} is served by more than one registry, the
     * earlier-listed registry wins, matching the priority order in which the registries are
     * declared (first listed = highest priority).
     *
     * @return a map of plugin id to its merged {@link PluginInfo}
     */
    @Override
    Map<String, PluginInfo> getPluginsMap() {
        final result = new LinkedHashMap<String, PluginInfo>()
        for( UpdateRepository repo : getRepositories() ) {
            for( Map.Entry<String, PluginInfo> entry : repo.getPlugins().entrySet() ) {
                final merged = result.get(entry.key)
                if( merged == null )
                    result.put(entry.key, copyPluginInfo(entry.value))
                else
                    mergeReleases(merged, entry.value)
            }
        }
        // a plugins.json index lists a moved plugin only under its plain name
        for( Map.Entry<String,String> it : aliases.entrySet() ) {
            if( !result.containsKey(it.value) && result.containsKey(it.key) )
                result.put(it.value, result.get(it.key))
        }
        return result
    }

    /**
     * Create a copy of a {@link PluginInfo} with its own releases list, so that merging across
     * repositories does not mutate the metadata cached by each repository.
     */
    private static PluginInfo copyPluginInfo(PluginInfo src) {
        final copy = new PluginInfo()
        copy.id = src.id
        copy.name = src.name
        copy.description = src.description
        copy.provider = src.provider
        copy.projectUrl = src.projectUrl
        copy.repositoryId = src.repositoryId
        copy.releases = src.releases != null
            ? new ArrayList<PluginInfo.PluginRelease>(src.releases)
            : new ArrayList<PluginInfo.PluginRelease>()
        return copy
    }

    /**
     * Append the releases of {@code src} to {@code target}, skipping any version already present.
     * Releases already in {@code target} come from an earlier-listed (higher priority) registry
     * and therefore take precedence on a version clash.
     */
    private static void mergeReleases(PluginInfo target, PluginInfo src) {
        if( !src.releases )
            return
        final knownVersions = new HashSet<String>()
        for( PluginInfo.PluginRelease r : target.releases )
            knownVersions.add(r.version)
        for( PluginInfo.PluginRelease r : src.releases )
            if( knownVersions.add(r.version) )
                target.releases.add(r)
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
     * Prefetch metadata for plugins. This gives an opportunity for certain
     * repository types to perform some data-loading optimisations.
     */
    void prefetchMetadata(List<PluginRef> plugins) {
        // Skip plugins that are already installed at the requested pinned version: no remote
        // metadata is needed to start them. This eliminates the registry round-trip on every
        // Nextflow invocation in the common CI case (pinned versions, plugins already cached).
        final needed = plugins.findAll { ref -> !isAlreadyInstalled(ref) || hasUnlinkedTwin(ref, plugins) }
        if( !needed ) {
            log.trace "All requested plugins are already installed - skipping registry metadata prefetch"
            return
        }
        prefetch0(needed)
    }

    private void prefetch0(List<PluginRef> needed) {
        // use direct field access because getRepositories() calls refresh() when the repository
        // list is null (e.g. in the dev updater), which loads `repositories.json` from the
        // working directory or fails if it does not exist
        for( def repo : this.@repositories ) {
            if( repo instanceof PrefetchUpdateRepository ) {
                log.trace "Prefetching plugin metadata from repository: ${repo.getClass().getSimpleName()} [${repo.id}]; plugins=${needed}"
                repo.prefetch(needed)
            }
        }
        readCanonicalIds()
    }

    // `nf-foo` and `myorg/nf-foo` may be one moved plugin; only the registry knows
    private boolean hasUnlinkedTwin(PluginRef ref, List<PluginRef> plugins) {
        return plugins.any { it -> isUnlinkedPair(ref.id, it.id) || isUnlinkedPair(it.id, ref.id) }
    }

    private boolean isUnlinkedPair(String plain, String scoped) {
        return scoped.endsWith('/' + plain) && aliases.get(plain) != scoped
    }

    /**
     * Fetch the metadata of a plugin that was not prefetched, such as a missing
     * dependency of a plugin that was already installed.
     *
     * @param id The plugin Id
     */
    protected void fetchMissingMetadata(String id) {
        if( offline )
            return
        // use direct field access for the same reason as prefetchMetadata()
        for( def repo : this.@repositories ) {
            try {
                // HttpPluginRepository fetches and caches the metadata when missing
                repo instanceof HttpPluginRepository ? repo.fetchPlugin(id) : repo.getPlugin(id)
            }
            catch( PluginRuntimeException e ) {
                log.warn "Failed to fetch metadata for plugin ${id} from ${repo.url} - cause: ${e.message}"
            }
        }
        readCanonicalIds()
    }

    protected boolean isAlreadyInstalled(PluginRef ref) {
        // in dev mode, a locally built plugin is used regardless of the requested version
        if( pluginManager.runtimeMode == RuntimeMode.DEVELOPMENT && pluginManager.getPlugin(ref.id) )
            return true
        // Only skip when the user has pinned a version: an unpinned spec needs remote metadata
        // to resolve the latest release.
        if( !ref.version )
            return false
        // Check the on-disk plugin store rather than the runtime PluginManager: at prefetch
        // time the local plugins have not been loaded into the manager yet (its per-run root
        // is still empty), so pluginManager.getPlugin() would always return null. installPlugin()
        // reuses the cached copy whenever the store directory exists (it only downloads when
        // missing), so its presence is the authoritative signal that no remote metadata is
        // required to start the plugin.
        return pluginsStore != null && FilesEx.exists(pluginsStore.resolve("${PluginRef.dirName(ref.id)}-${ref.version}"))
    }

    /**
     * Resolve a plugin installing or updating the dependencies if necessary
     * and start the plugin
     *
     * @param pluginId The plugin Id
     * @param version The required version, can be null for the latest version
     */
    void prepareAndStart(String pluginId, String version) {
        pluginId = canonical(pluginId)
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
            final specs = plugins.collect(it -> PluginRef.parse(it,defaultPlugins))
            // installed plugins too, to add missing links
            prefetch0(specs)
            for( PluginRef spec : specs ) {
                pullPlugin0(spec.id, spec.version)
            }
        }
        finally {
            pullOnly=false
        }
    }

    void pullPlugin0(String pluginId, String version) {
        pluginId = canonical(pluginId)
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
        id = canonical(id)
        // in dev mode, core plugins must be built locally rather than downloaded
        if( pluginManager.runtimeMode == RuntimeMode.DEVELOPMENT && defaultPlugins.hasPlugin(id) )
            throw new UnsupportedOperationException("Core plugin $id is not available in dev mode -- make sure it is built locally (e.g. `make compile`)")
        if( !pluginsStore )
            throw new IllegalStateException("Missing pluginStore attribute")

        return load0(id, version)
    }

    private Path download0(String id, String version) {
        // 0. check if version is specified
        if( !version )
            throw new InvalidPluginDescriptorException("Missing version for plugin $id")
        log.info "Downloading plugin ${id}@${version}"

        // 1. check if already exists
        final pluginPath = pluginsStore.resolve("${PluginRef.dirName(id)}-$version")
        if( FilesEx.exists(pluginPath) ) {
            return pluginPath
        }

        // 2. download to temporary location
        Path downloaded = safeDownloadPlugin(id, version);

        // 3. unzip the content and delete downloaded file
        Path dir = FileUtils.expandIfZip(downloaded)
        FileHelper.deletePath(downloaded)

        // 4. move the final destination the plugin directory
        // the zip of a namespaced plugin is named <name>-<version>; safeMove renames it
        assert id.contains('/') || pluginPath.getFileName() == dir.getFileName()
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
                log.debug("Failed to download plugin: $id; version: $version - attempt: ${event.attemptCount}", event.lastException)
            }
        }

        final condition = new CheckedPredicate<Throwable>() {
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
        return withLock(id, version) { download0(id, version) }
    }

    private <T> T withLock(String id, String version, Closure<T> action) {
        final sentinel = lockFile(id,version)
        final mutex = new FileMutex(target: sentinel, timeout: '10 min')
        try {
            return mutex.lock(action)
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
        new File(tmp, "nextflow-plugin-${PluginRef.dirName(id)}-${version}.lock")
    }

    private boolean load0(String id, String requestedVersion) {
        assert id, "Missing plugin Id"
        id = canonical(id)

        if( offline && !requestedVersion ) {
            throw new IllegalStateException("Cannot find version for $id plugin -- plugin versions MUST be specified in offline mode")
        }

        if( !isAlreadyInstalled(new PluginRef(id, requestedVersion)) ) {
            fetchMissingMetadata(id)
            id = canonical(id)
        }

        def version = requestedVersion
        if( !version ) {
            version = getLastPluginRelease(id)?.version
        }
        else if( !Version.isValid(version) ) {
            // a version is 'valid' if it's an exact semver version "major.minor.patch" so
            // if it's not that, treat it as a version constraint and look for matches
            version = findNewestMatchingRelease(id, version)?.version
        }

        if( !version ) {
            final msg = requestedVersion
                ? "Cannot find version of $id plugin matching $requestedVersion"
                : "Cannot find latest version of $id plugin"
            throw new IllegalStateException(msg)
        }

        if( version != requestedVersion ) {
            log.debug "Plugin $id version $requestedVersion resolved to: $version"
        }

        final folderName = "${PluginRef.dirName(id)}-$version".toString()
        def pluginPath = pluginsStore.resolve(folderName)
        final plainPath = isMoved(id) ? pluginsStore.resolve("${PluginRef.name(id)}-$version") : null
        // read-only store: load the plain folder of an older Nextflow when the link cannot be written
        if( plainPath && !FilesEx.exists(pluginPath) && !linkToPlain(pluginPath, plainPath, id, version) && isComplete(plainPath) )
            pluginPath = plainPath
        if( !FilesEx.exists(pluginPath) )
            pluginPath = safeDownload(id, version)
        if( plainPath && isComplete(pluginPath) && pluginPath != plainPath )
            createLink(plainPath, pluginPath)

        // verify the plugin install path contains the expected manifest path
        if( !FilesEx.exists(pluginPath.resolve('classes/META-INF/MANIFEST.MF')) ) {
            log.warn("Plugin '${pluginPath.getFileName()}' installation looks corrupted - Delete the following directory and run nextflow again: $pluginPath")
        }

        // load the plugin from the file system, unless it is already loaded under the namespace its zip declares
        final declared = new NamespaceDescriptorFinder().find(pluginPath, folderName).pluginId
        PluginWrapper wrapper = declared != id && PluginRef.name(declared) == id ? pluginManager.getPlugin(declared) : null
        if( wrapper && wrapper.descriptor.version != version )
            throw new PluginRuntimeException("Plugin $id@$version cannot be loaded -- its zip declares the ID $declared, already loaded with version ${wrapper.descriptor.version}")
        if( !wrapper )
            wrapper = pluginManager.loadPluginFromPath(pluginPath, folderName)
        if( wrapper.pluginId != id ) {
            if( id.contains('/') ) {
                pluginManager.unloadUnresolved(wrapper.pluginId)
                throw new PluginRuntimeException("Plugin $id cannot be loaded -- its zip declares the ID ${wrapper.pluginId}")
            }
            if( PluginRef.name(wrapper.pluginId) == id ) {
                aliases.put(id, wrapper.pluginId)
                linkToPlain(pluginsStore.resolve("${PluginRef.dirName(wrapper.pluginId)}-$version"), pluginPath, wrapper.pluginId, version)
                id = wrapper.pluginId
            }
        }

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
            def depVersion = defaultPlugins.hasPlugin(it.pluginId) ? defaultPlugins.getPlugin(it.pluginId).version : null
            // -- otherwise try to find the newest matching release
            if( !depVersion ) {
                fetchMissingMetadata(it.pluginId)
                depVersion = findNewestMatchingRelease(it.pluginId, it.pluginVersionSupport)?.version
            }
            log.debug "Plugin $id requires $it.pluginId supported version: $it.pluginVersionSupport - available version: $depVersion"
            if( pullOnly )
                pullPlugin0(it.pluginId, depVersion)
            else
                prepareAndStart(it.pluginId, depVersion)
        }

        if( pullOnly )
            return false

        // pf4j matches dependency IDs exactly; `set`, not `replaceAll`: replaceAll breaks iterators on other threads
        for( PluginWrapper it : pluginManager.plugins ) {
            final list = it.descriptor.dependencies
            for( int i=0; i<list.size(); i++ ) {
                final dep = canonical(list[i])
                if( !dep.is(list[i]) )
                    list.set(i, dep)
            }
        }
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
        id = canonical(id)
        if (pluginManager.getPlugin(id) == null) {
            throw new PluginRuntimeException("Plugin $id cannot be updated since it is not installed")
        }

        fetchMissingMetadata(id)
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
        if( offline ) {
            log.debug "Update not supported in offline mode"
            return false
        }

        pluginId = canonical(pluginId)
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
        final current = pluginManager.getPlugin(canonical(id))
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

        PluginInfo pluginInfo = getPluginsMap().get(canonical(id))
        if( !pluginInfo )
            throw new IllegalArgumentException("Unknown plugin id: $id")

        // note: order releases list by descending version numbers ie. latest version comes first
        def releases = pluginInfo.releases.sort(false) { a,b -> Version.parse(b.version) <=> Version.parse(a.version) }
        for (PluginInfo.PluginRelease rel : releases ) {
            if( !versionManager.checkVersionConstraint(rel.version, verConstraint) || !rel.url )
                continue

            if( versionManager.checkVersionConstraint(BuildInfo.version, rel.requires) )
                return rel
        }

        return null
    }

    private PluginDependency canonical(PluginDependency dep) {
        final id = canonical(dep.pluginId)
        return id == dep.pluginId ? dep : new PluginDependency("${id}${dep.optional ? '?' : ''}@${dep.pluginVersionSupport}")
    }
}
