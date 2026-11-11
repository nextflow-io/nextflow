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

package nextflow.cli.pipeline

import java.nio.file.Path

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.ParameterException
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.cli.CmdBase
import nextflow.config.ConfigCmdAdapter
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig

/**
 * Helpers shared by the pipeline subcommands
 */
@CompileStatic
class PipelineCmdSupport {

    static final String OUTPUT_MODE_HELP = 'Output mode: text, simple (search), json'

    /** Accepts the output modes of list/view (text) and search (simple), plus json */
    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('text', 'simple', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    static RegistryConfig registryConfig(CmdBase cmd, Path baseDir, List<Path> configFiles = []) {
        final config = new ConfigCmdAdapter()
            .setOptions(cmd.launcher.options)
            .setUserConfigFiles(configFiles)
            .setBaseDir(baseDir)
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
    static RegistryAddress registryAddress(String name, String hub, RegistryConfig config) {
        final isUrl = name.startsWith('http://') || name.startsWith('https://')
        hub = hub ?: SysEnv.get('NXF_HUB')
        // cheap checks first, so that Git addresses don't read the config
        if( isUrl ? !name.contains('/pipelines/') : !hub )
            return null
        final providers = isUrl ? null : ProviderConfig.getDefaultProviderNames()
        if( !isUrl && hub in providers )
            return null
        // the default public registry is always recognized, even when not configured
        final registries = config.allUrls + [RegistryConfig.DEFAULT_REGISTRY_URL]
        // a URL and a hub name a registry by its host
        final api = registryOfHost(isUrl ? URI.create(name).host : hub, registries)
        // an address names one registry, and the API key belongs to the primary one
        final apiKey = api && api == config.url ? config.apiKey : null
        if( isUrl ) {
            // a registry URL must be a pipeline page
            final prefix = api ? pageUrl(api) + '/pipelines/' : null
            if( !prefix || !name.startsWith(prefix) )
                return null
            return new RegistryAddress(api, name.substring(prefix.length()).replaceFirst(/\/$/, ''), apiKey)
        }
        if( !api ) {
            final hosts = registries.collect { URI.create(it).host }.unique()
            throw new AbortOperationException("Unknown repository provider or registry: `$hub` -- providers: ${providers.join(', ')}; registries: ${hosts.join(', ')}")
        }
        return new RegistryAddress(api, name, apiKey)
    }

    /** Base URL of the registry web pages: the API URL without {@code /api} */
    static String pageUrl(String apiUrl) {
        return apiUrl.replaceFirst(/\/api\/?$/, '')
    }

    private static String registryOfHost(String host, Collection<String> registryUrls) {
        return registryUrls.find { URI.create(it).host?.equalsIgnoreCase(host) }
    }
}
