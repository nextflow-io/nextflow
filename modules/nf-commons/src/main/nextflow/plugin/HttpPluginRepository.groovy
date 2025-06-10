package nextflow.plugin

import com.google.gson.Gson
import dev.failsafe.Failsafe
import dev.failsafe.FailsafeExecutor
import dev.failsafe.Fallback
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import org.pf4j.PluginRuntimeException
import org.pf4j.update.FileDownloader
import org.pf4j.update.FileVerifier
import org.pf4j.update.PluginInfo
import org.pf4j.update.SimpleFileDownloader
import org.pf4j.update.verifier.CompoundVerifier

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse

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
    private final HttpClient client = HttpClient.newHttpClient()
    private final String id
    private final URI url

    private Map<String, PluginInfo> plugins = new HashMap<>()

    HttpPluginRepository(String id, URI url) {
        this.id = id
        // ensure url ends with a slash
        this.url = !url.toString().endsWith("/")
            ? URI.create(url.toString() + "/")
            : url
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
        if (plugins.isEmpty()) {
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
        final CheckedSupplier<Map<String, PluginInfo>> supplier = () -> fetchMetadata0(ordered)
        return retry().get(supplier)
    }

    private Map<String, PluginInfo> fetchMetadata0(List<PluginSpec> specs) {
        final gson = new Gson()

        def reqBody = gson.toJson([
            'nextflowVersion': BuildInfo.version,
            'plugins'        : specs
        ])

        def req = HttpRequest.newBuilder()
            .uri(url.resolve("plugins/collect"))
            .POST(HttpRequest.BodyPublishers.ofString(reqBody))
            .build()

        def rep = client.send(req, HttpResponse.BodyHandlers.ofString())
        if (rep.statusCode() != 200) throw new PluginRuntimeException(errorMessage(rep, gson))

        try {
            def repBody = gson.fromJson(rep.body(), FetchResponse)
            return repBody.plugins.collectEntries { p -> Map.entry(p.id, p) }
        } catch (Exception e) {
            log.info("Plugin metadata response body: '${rep.body()}'")
            throw new PluginRuntimeException("Failed to parse response body", e)
        }
    }

    // create a retry executor using failsafe
    private static FailsafeExecutor retry() {
        EventListener<ExecutionAttemptedEvent> logAttempt = (ExecutionAttemptedEvent attempt) -> {
            log.debug("Retrying download of plugins metadata - attempt ${attempt.attemptCount}, ${attempt.lastFailure.message}", attempt.lastFailure)
        }
        Fallback fallback = Fallback.ofException { e ->
            e.lastFailure instanceof ConnectException
                ? new ConnectException("Failed to download plugins metadata")
                : new PluginRuntimeException("Failed to download plugin metadata: ${e.lastFailure.message}")
        }
        final policy = RetryPolicy.builder()
            .withMaxAttempts(3)
            .handle(ConnectException)
            .onRetry(logAttempt)
            .build()
        return Failsafe.with(fallback, policy)
    }

    private static String errorMessage(HttpResponse<String> rep, Gson gson) {
        try {
            def err = gson.fromJson(rep.body(), ErrorResponse)
            return "${err.type} - ${err.message}"
        } catch (Exception e) {
            return rep.body()
        }
    }

    // ---------------------

    /**
     * Response format object expected from repository
     */
    private static class FetchResponse {
        List<PluginInfo> plugins
    }

    private static class ErrorResponse {
        String type
        String message
    }
}
