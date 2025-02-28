package nextflow.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx
import org.pf4j.ManifestPluginDescriptorFinder
import org.pf4j.PluginDescriptor
import org.pf4j.update.FileDownloader
import org.pf4j.update.FileVerifier
import org.pf4j.update.PluginInfo
import org.pf4j.update.UpdateRepository
import org.pf4j.update.verifier.CompoundVerifier

import java.nio.file.Path

/**
 * Implementation of UpdateRepository which looks in a local directory of already-downloaded
 * plugins to find available versions.
 */
@CompileStatic
@Slf4j
class LocalUpdateRepository implements UpdateRepository {
    private final String id
    private final Path dir
    private Map<String, PluginInfo> plugins

    LocalUpdateRepository(String id, Path dir) {
        this.id = id
        this.dir = dir
    }

    @Override
    String getId() {
        return id
    }

    @Override
    URL getUrl() {
        return dir.toUri().toURL()
    }

    @Override
    Map<String, PluginInfo> getPlugins() {
        if( !plugins )
            plugins = loadPlugins(dir)
        return plugins
    }

    @Override
    PluginInfo getPlugin(String id) {
        return getPlugins().get(id)
    }

    @Override
    void refresh() {
        this.plugins = null
    }

    @Override
    FileDownloader getFileDownloader() {
        // plugins in this repo are already downloaded, so treat any download url as a file path
        return (URL url) -> Path.of(url.toURI())
    }

    @Override
    FileVerifier getFileVerifier() {
        return new CompoundVerifier()
    }

    private static Map<String, PluginInfo> loadPlugins(Path dir) {
        // each plugin is stored in a dir called $id-$version; grab the descriptor from each
        final manifestReader = new ManifestPluginDescriptorFinder()
        final descriptors = FilesEx.listFiles(dir)
            .collect { plugin -> new LocalPlugin(plugin, manifestReader.find(plugin)) }

        // now group the descriptors by id, to create a PluginInfo with list of versions
        return descriptors
            .groupBy { d -> d.getPluginId() }
            .collectEntries { id, versions -> Map.entry(id, toPluginInfo(id, versions)) }
    }

    private static PluginInfo toPluginInfo(String id, List<LocalPlugin> versions) {
        final info = new PluginInfo()
        info.id = id
        info.releases = versions.collect { v ->
            final release = new PluginInfo.PluginRelease()
            release.version = v.version
            release.requires = v.requires
            release.url = v.path.toUri().toURL()
            if( !info.provider && v.provider )
                info.provider = v.provider
            return release
        }
        return info
    }

    private static class LocalPlugin {
        Path path

        @Delegate
        PluginDescriptor descriptor

        LocalPlugin(Path path, PluginDescriptor descriptor) {
            this.path = path
            this.descriptor = descriptor
        }
    }
}
