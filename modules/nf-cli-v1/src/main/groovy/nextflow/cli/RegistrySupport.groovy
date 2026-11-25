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

package nextflow.cli

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.client.RegistryClient
import nextflow.SysEnv
import nextflow.config.ConfigCmdAdapter
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.PipelineSource
import nextflow.module.RegistryClientFactory
import nextflow.scm.ProviderConfig

/**
 * Helpers for registry pipeline addresses in the project commands
 */
@Slf4j
@CompileStatic
class RegistrySupport {

    /**
     * Registry config from the user config, the `nextflow.config` in the launch directory, and the given config files
     */
    static RegistryConfig registryConfig(CmdBase cmd, List<Path> configFiles = []) {
        final config = new ConfigCmdAdapter()
            .setOptions(cmd.launcher.options)
            .setUserConfigFiles(configFiles)
            .setBaseDir(Path.of('.').toAbsolutePath().normalize())
            .build()
        return config.navigate('registry') as RegistryConfig ?: new RegistryConfig()
    }

    /** A registry pipeline address: the registry API URL, the pipeline name, and the API key to send */
    @Canonical
    static class RegistryAddress {
        String url
        String name
        String apiKey
    }

    /**
     * The registry addressed by a registry page URL, or by a registry host as hub ({@code -hub} or {@code NXF_HUB}).
     * A hub that is a Git provider name keeps addressing Git, as before.
     *
     * @return the registry address, or {@code null} when {@code name} and {@code hub} address Git
     */
    static RegistryAddress registryAddress(String name, String hub, Closure<RegistryConfig> registryConfig) {
        final isUrl = name.startsWith('http://') || name.startsWith('https://')
        hub = hub ?: SysEnv.get('NXF_HUB')
        // cheap checks first, so that Git addresses never read the config
        if( !isUrl && !hub )
            return null
        final uri = isUrl ? toUri(name) : null
        // a registry URL is a pipeline page
        if( isUrl && !uri?.path?.contains('/pipelines/') )
            return null
        // a Git provider host is known without the launch config, e.g. a GitLab group named `pipelines`
        if( isUrl && isGitProviderHost(uri) )
            return null
        final providers = isUrl ? null : ProviderConfig.getDefaultProviderNames()
        if( !isUrl && hub in providers )
            return null
        final config = registryConfig.call()
        // the default public registry is always recognized, even when not configured
        final registries = config.allUrls + [RegistryConfig.DEFAULT_REGISTRY_URL]
        // a URL and a hub name a registry by its host, or host and port
        final api = registryOfHost(isUrl ? uri.authority : hub, registries)
        // an address names one registry, and the API key belongs to the primary one
        final apiKey = api && api == config.url ? config.apiKey : null
        if( isUrl ) {
            if( !api )
                return null
            // a registry URL must be a pipeline page; the scheme is ignored since requests use the configured API URL
            final prefix = (URI.create(pageUrl(api)).path ?: '') + '/pipelines/'
            if( !uri.path.startsWith(prefix) )
                return null
            return newAddress(api, uri.path.substring(prefix.length()).replaceFirst(/\/$/, ''), apiKey)
        }
        if( !api ) {
            final hosts = registries.collect { URI.create(it).authority }.unique()
            throw new AbortOperationException("Unknown repository provider or registry: `$hub` -- providers: ${providers.join(', ')}; registries: ${hosts.join(', ')}")
        }
        return newAddress(api, name, apiKey)
    }

    // checked here rather than when resolving the release, so that it applies offline too
    private static RegistryAddress newAddress(String api, String name, String apiKey) {
        if( name.contains('@') )
            throw new AbortOperationException("Not a valid registry pipeline address: ${name} -- give the version with the `-r` option")
        return new RegistryAddress(api, name, apiKey)
    }

    private static boolean isGitProviderHost(URI uri) {
        final domains = ProviderConfig.getDefaultProviderDomains()
        return domains.any { String it -> it.equalsIgnoreCase(uri.host) || it.equalsIgnoreCase(uri.authority) }
    }

    /**
     * Resolve the release of a registry address, only against the registry it names.
     *
     * @param revision The registry version, or {@code null} for the latest release
     */
    static PipelineSource resolveRelease(RegistryAddress address, String revision, RegistryClient client = null) {
        final release = PipelineSource.resolve(client ?: RegistryClientFactory.forUrls([address.url], address.apiKey), address.name, revision)
        log.info "Using registry release ${address.name}@${release.version}" + (release.commitId ? " (commit ${release.commitId.take(10)})" : '')
        return release
    }

    /** Base URL of the registry web pages: the API URL without {@code /api} */
    static String pageUrl(String apiUrl) {
        return apiUrl.replaceFirst(/\/api\/?$/, '')
    }

    private static String registryOfHost(String host, Collection<String> registryUrls) {
        return registryUrls.find { final uri = URI.create(it); uri.host?.equalsIgnoreCase(host) || uri.authority?.equalsIgnoreCase(host) }
    }

    /** A URL that {@link URI} can't parse is not a registry URL, so it is left to Git */
    private static URI toUri(String url) {
        try {
            return URI.create(url)
        }
        catch( IllegalArgumentException e ) {
            return null
        }
    }
}
