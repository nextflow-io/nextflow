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
 */

package nextflow.packages

import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.ISession
import nextflow.NextflowMeta
import nextflow.plugin.Plugins

/**
 * Manages package providers and coordinates package environment creation
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Slf4j
@CompileStatic
class PackageManager {

    private final Map<String, PackageProvider> providers = new ConcurrentHashMap<>()
    private final ISession session
    
    PackageManager(ISession session) {
        this.session = session
        initializeProviders()
    }

    /**
     * Initialize available package providers from plugins
     */
    private void initializeProviders() {
        // Load package providers from plugins
        def extensions = Plugins.getExtensions(PackageProviderExtension)
        for (PackageProviderExtension extension : extensions) {
            def provider = extension.createProvider(session)
            if (provider && provider.isAvailable()) {
                providers.put(provider.getName(), provider)
                log.debug "Registered package provider: ${provider.getName()}"
            }
        }
    }

    /**
     * Get a package provider by name
     * 
     * @param name Provider name (e.g., "conda", "pixi")
     * @return The package provider or null if not found
     */
    PackageProvider getProvider(String name) {
        return providers.get(name)
    }

    /**
     * Get all available package providers
     * 
     * @return Map of provider name to provider instance
     */
    Map<String, PackageProvider> getProviders() {
        return Collections.unmodifiableMap(providers)
    }

    /**
     * Create a package environment using the appropriate provider
     * 
     * @param spec The package specification
     * @return The path to the created environment
     */
    Path createEnvironment(PackageSpec spec) {
        if (!spec.isValid()) {
            throw new IllegalArgumentException("Invalid package specification: ${spec}")
        }

        def provider = getProvider(spec.provider)
        if (!provider) {
            throw new IllegalArgumentException("Package provider not found: ${spec.provider}")
        }

        if (!provider.supportsSpec(spec)) {
            throw new IllegalArgumentException("Package specification not supported by provider ${spec.provider}: ${spec}")
        }

        return provider.createEnvironment(spec)
    }

    /**
     * Get the activation script for an environment
     * 
     * @param spec The package specification
     * @param envPath Path to the environment
     * @return Shell script snippet to activate the environment
     */
    String getActivationScript(PackageSpec spec, Path envPath) {
        def provider = getProvider(spec.provider)
        if (!provider) {
            throw new IllegalArgumentException("Package provider not found: ${spec.provider}")
        }

        return provider.getActivationScript(envPath)
    }

    /**
     * Parse a package specification from process configuration
     * 
     * @param packageDef Package definition (string or map)
     * @param provider Default provider if not specified
     * @return Parsed package specification
     */
    static PackageSpec parseSpec(Object packageDef, String provider = null) {
        // a plain value, e.g. `package "numpy pandas"` (String or GString)
        if (packageDef instanceof CharSequence) {
            return new PackageSpec(provider, [packageDef.toString()])
        }
        if (packageDef instanceof List) {
            def list = packageDef as List
            // Groovy named-args directive form: `package "numpy", provider: "uv"`
            // is delivered as [ [provider:'uv', ...], "numpy", ... ] with the
            // options map as the leading element -- merge it and recurse.
            if (list && list[0] instanceof Map) {
                def opts = new LinkedHashMap(list[0] as Map)
                def values = list.size() > 1 ? list.subList(1, list.size()) : []
                if (!opts.containsKey('packages') && !opts.containsKey('environment') && values)
                    opts.put('packages', values.size() == 1 ? values[0] : values)
                return parseSpec(opts, provider)
            }
            return new PackageSpec(provider, list.collect { it?.toString() } as List<String>)
        }
        if (packageDef instanceof Map) {
            def map = packageDef as Map
            def spec = new PackageSpec()

            if (map.containsKey('provider')) {
                spec.provider = map.provider as String
            } else if (provider) {
                spec.provider = provider
            }

            if (map.containsKey('packages')) {
                def packages = map.packages
                if (packages instanceof List) {
                    spec.entries = packages.collect { it?.toString() } as List<String>
                } else if (packages != null) {
                    spec.entries = [packages.toString()]
                }
            }

            if (map.containsKey('environment')) {
                spec.environment = map.environment?.toString()
            }

            if (map.containsKey('channels')) {
                def channels = map.channels
                if (channels instanceof List) {
                    spec.channels = channels as List<String>
                } else if (channels != null) {
                    spec.channels = [channels.toString()]
                }
            }

            if (map.containsKey('options')) {
                spec.options = map.options as Map<String, Object>
            }

            return spec
        }

        throw new IllegalArgumentException("Invalid package definition: ${packageDef}")
    }

    /**
     * Auto-detect a package specification by scanning a module directory for a
     * provider manifest file (e.g. {@code environment.yml}, {@code requirements.txt}).
     *
     * Only providers that are available on the system are considered; the manifest
     * file names are contributed by each provider via {@link PackageProvider#getManifestFileNames()}.
     *
     * @param moduleDir the process module directory to scan
     * @return a {@link PackageSpec} referencing the detected manifest file, or {@code null}
     *         if no provider manifest is found
     */
    PackageSpec detectSpec(Path moduleDir) {
        if( !moduleDir )
            return null
        final manifests = new LinkedHashMap<String, List<String>>()
        for( String name : providers.keySet() )
            manifests.put(name, providers.get(name).getManifestFileNames())
        return findManifestSpec(moduleDir, manifests)
    }

    /**
     * Pure detection helper: given a map of provider name to the manifest file
     * names it supports, return a spec for the first manifest found in {@code moduleDir}.
     * Providers are scanned in a deterministic (name-sorted) order.
     *
     * @param moduleDir the directory to scan
     * @param providerManifests map of provider name to its manifest file names
     * @return a {@link PackageSpec} for the first match, or {@code null}
     */
    @PackageScope
    static PackageSpec findManifestSpec(Path moduleDir, Map<String, List<String>> providerManifests) {
        if( !moduleDir || !providerManifests )
            return null
        for( String provider : providerManifests.keySet().sort() ) {
            final files = providerManifests.get(provider)
            if( !files )
                continue
            for( String file : files ) {
                final candidate = moduleDir.resolve(file)
                if( candidate.exists() ) {
                    log.debug "Auto-detected ${provider} manifest file: ${candidate}"
                    final spec = new PackageSpec()
                    spec.provider = provider
                    spec.environment = candidate.toAbsolutePath().toString()
                    return spec
                }
            }
        }
        return null
    }

    /**
     * Check if the package manager feature is enabled, either by the
     * {@code nextflow.preview.package} feature flag in the pipeline script
     * or by the same setting in the configuration file
     *
     * @param session The current session
     * @return True if the feature is enabled
     */
    static boolean isEnabled(ISession session) {
        if( NextflowMeta.instance.preview.getPackage() )
            return true
        return session.config.navigate('nextflow.preview.package', false) as Boolean
    }
}