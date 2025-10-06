package nextflow.plugin


import java.net.http.HttpRequest
import java.net.http.HttpResponse

import nextflow.serde.gson.GsonEncoder
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import io.seqera.npr.api.schema.v1.ListDependenciesResponse
import io.seqera.npr.api.schema.v1.Plugin
import nextflow.BuildInfo
import nextflow.util.RetryConfig
import org.pf4j.PluginRuntimeException
import org.pf4j.update.FileDownloader
import org.pf4j.update.FileVerifier
import org.pf4j.update.PluginInfo
import org.pf4j.update.PluginInfo.PluginRelease
import org.pf4j.update.SimpleFileDownloader
import org.pf4j.update.verifier.CompoundVerifier
/**
 * Represents an update repository served via an HTTP api.
 *
 * It implements PrefetchUpdateRepository so that all relevant
 * plugin metadata can be loaded with a single HTTP request, rather
 * than a request-per-plugin.
 *
 * Metadata is prefetched into memory when Nextflow starts and expires
 * upon termination (or when 'refresh()' is called).
 */
@Slf4j
@CompileStatic
class HttpPluginRepository implements PrefetchUpdateRepository {
    private final String id
    private final URI url
    private final HxClient httpClient

    private Map<String, PluginInfo> plugins

    HttpPluginRepository(String id, URI url) {
        this.id = id
        // ensure url ends with a slash
        this.url = !url.toString().endsWith("/")
            ? URI.create(url.toString() + "/")
            : url
        this.httpClient = HxClient.newBuilder()
                .retryConfig(RetryConfig.config())
                .build()
    }

    // NOTE ON PREFETCHING
    //
    // The prefetch mechanism is used to work around a limitation in the
    // UpdateRepository interface from pf4j.
    //
    // Specifically, p4fj expects that getPlugins() returns a Map<> of all
    // metadata about all plugins. To implement this for an HTTP repository
    // would require either downloading the entire contents of the remote
    // repository or implementing a lazy map and making an HTTP request for
    // each required plugin.
    //
    // Instead we can use the list of configured plugins to load all relevant
    // metadata in a single HTTP request at startup, and use this to populate
    // the map. Once the prefetch is complete, this repository will behave
    // like any other implementation of UpdateRepository.
    @Override
    void prefetch(List<PluginSpec> plugins) {
        if (plugins && !plugins.isEmpty()) {
            this.plugins = fetchMetadata(plugins)
        }
    }

    @Override
    String getId() {
        return id
    }

    @Override
    URL getUrl() {
        return url.toURL()
    }

    @Override
    Map<String, PluginInfo> getPlugins() {
        if (plugins==null) {
            log.warn "getPlugins() called before prefetch() - plugins map will be empty"
            return Map.of()
        }
        return Collections.unmodifiableMap(plugins)
    }

    @Override
    PluginInfo getPlugin(String id) {
        return plugins.computeIfAbsent(id) { key -> fetchMetadataByIds([key]).get(key) }
    }

    @Override
    void refresh() {
        plugins = fetchMetadataByIds(plugins.keySet())
    }

    @Override
    FileDownloader getFileDownloader() {
        return new SimpleFileDownloader()
    }

    @Override
    FileVerifier getFileVerifier() {
        return new CompoundVerifier()
    }

    // ----------------------------------------------------------------------------
    // http handling

    private Map<String, PluginInfo> fetchMetadataByIds(Collection<String> ids) {
        def specs = ids.collect(id -> new PluginSpec(id, null))
        return fetchMetadata(specs)
    }

    private Map<String, PluginInfo> fetchMetadata(Collection<PluginSpec> specs) {
        final ordered = specs.sort(false)
        return fetchMetadata0(ordered)
    }

    private Map<String, PluginInfo> fetchMetadata0(List<PluginSpec> specs) {
        def pluginsParam = specs.collect { "${it.id}${it.version ? '@' + it.version : ''}" }.join(',')
        def uri = url.resolve("v1/plugins/dependencies?plugins=${URLEncoder.encode(pluginsParam, 'UTF-8')}&nextflowVersion=${URLEncoder.encode(BuildInfo.version, 'UTF-8')}")
        def req = HttpRequest.newBuilder()
            .uri(uri)
            .GET()
            .build()
        try {
            return sendAndParse(req)
        }
        catch (PluginRuntimeException e) {
            throw e
        }
        catch (Exception e) {
            throw new PluginRuntimeException(e, "Unable to connect to ${uri} - cause: ${e.message}")
        }
    }

    private Map<String, PluginInfo> sendAndParse(HttpRequest req) {
        final encoder = new GsonEncoder<ListDependenciesResponse>() {}
        final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
        final body = resp.body()
        log.debug "Registry request: ${resp.uri()}\n- code: ${resp.statusCode()}\n- body: ${body}"
        if( resp.statusCode() != 200 ) {
            final msg = "Invalid response while fetching plugin metadata from: ${req.uri()}\n- http status: ${resp.statusCode()}\n- response   : ${body}"
            throw new PluginRuntimeException(msg)
        }
        try {
            final ListDependenciesResponse decoded = encoder.decode(body)
            if( decoded.plugins == null ) {
                throw new PluginRuntimeException("Failed to download plugin metadata: Failed to parse response body")
            }
            final result = new HashMap<String, PluginInfo>()
            for( Plugin plugin : decoded.plugins ) {
                if( plugin.releases ) {
                    final pluginInfo = mapToPluginInfo(plugin)
                    result.put(plugin.id, pluginInfo)
                }
                else
                    log.debug "Registry ${resp.uri().host} has no releases for plugin: ${plugin}"
            }
            return result
        }
        catch( Exception e ) {
            final msg = "Unexpected error while fetching plugin metadata from: ${req.uri()}\n- message : ${e.message}\n- response: ${body}"
            throw new PluginRuntimeException(msg)
        }
    }

    /**
     * Maps a Plugin object from the repository API to a PluginInfo object for pf4j compatibility.
     * Handles conversion of OffsetDateTime to Date and ensures the releases collection is never null.
     *
     * @param plugin The Plugin object from the repository API response
     * @return A PluginInfo object compatible with pf4j's update repository interface
     */
    static protected PluginInfo mapToPluginInfo(Plugin plugin) {
        assert plugin.releases, "Plugin releases cannot be empty"

        final pluginInfo = new PluginInfo()
        pluginInfo.id = plugin.id
        pluginInfo.projectUrl = plugin.projectUrl
        pluginInfo.provider = plugin.provider
        
        // Map releases to PluginInfo.PluginRelease
        pluginInfo.releases = new ArrayList<>()
        for (def release : plugin.releases) {
            final pluginRelease = new PluginRelease()
            pluginRelease.version = release.version
            pluginRelease.date = release.date ? Date.from(release.date.toInstant()) : null
            pluginRelease.url = release.url
            pluginRelease.sha512sum = release.sha512sum
            pluginRelease.requires = release.requires
            pluginInfo.releases.add(pluginRelease)
        }
        
        return pluginInfo
    }

}
