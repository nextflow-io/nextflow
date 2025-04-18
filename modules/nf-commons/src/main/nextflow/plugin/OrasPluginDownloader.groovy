package nextflow.plugin

import groovy.transform.CompileStatic
import land.oras.ContainerRef
import land.oras.Manifest
import land.oras.Registry
import org.pf4j.PluginRuntimeException

import java.nio.file.Files
import java.nio.file.Path

/**
 * Support downloading plugins from an ORAS (OCI Registry As Storage)
 * location.
 *
 * NOTE: logically this should implement pf4j's FileDownloader interface,
 * but unfortunately that interface uses java.net.URL which doesn't support
 * custom URI protocols.
 *
 * See https://oras.land/
 */
@CompileStatic
class OrasPluginDownloader {
    private static final String PROTOCOL = "oras://"
    private static final String ARTIFACT_TYPE = "application/vnd.nextflow.plugin+zip"

    static boolean canDownload(String uri) {
        return uri.startsWith(PROTOCOL)
    }

    Path downloadFile(String uri) throws IOException {
        if (!canDownload(uri)) {
            throw new PluginRuntimeException("URI protocol not supported by ORAS downloader: {}", uri);
        }
        return downloadOras(uri)
    }

    private static Path downloadOras(String uri) {
        Path destination = Files.createTempDirectory("oras-update-downloader")
        destination.toFile().deleteOnExit()

        // strip the oras:// prefix since it's just a marker
        ContainerRef source = ContainerRef.parse(uri.replaceFirst("oras://", ""))
        Registry registry = Registry.Builder.builder()
            // use http on localhost, require https everywhere else
            .withInsecure(source.registry.startsWith("localhost:"))
            .build()

        // grab oras metadata about the url
        Manifest manifest = registry.getManifest(source)
        String type = manifest.getArtifactType()
        if (type != ARTIFACT_TYPE) {
            throw new PluginRuntimeException("Not a nextflow plugin: {}", uri)
        }
        // get the filename of the plugin artifact, which should be in layer 0
        def layers = manifest.getLayers()
        if (layers.size() < 1) {
            throw new PluginRuntimeException("Unable to find nextflow plugin at: {}", uri)
        }
        String filename = layers[0].annotations['org.opencontainers.image.title']

        // download the plugin artifact
        registry.pullArtifact(source, destination, true)
        return destination.resolve(filename)
    }
}
