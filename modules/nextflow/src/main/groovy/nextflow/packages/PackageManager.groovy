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
import groovy.util.logging.Slf4j
import nextflow.ISession
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
        if (packageDef instanceof String) {
            return new PackageSpec(provider, [packageDef])
        } else if (packageDef instanceof List) {
            return new PackageSpec(provider, packageDef as List<String>)
        } else if (packageDef instanceof Map) {
            def map = packageDef as Map
            def spec = new PackageSpec()
            
            if (map.containsKey('provider')) {
                spec.provider = map.provider as String
            } else if (provider) {
                spec.provider = provider
            }
            
            if (map.containsKey('packages')) {
                def packages = map.packages
                if (packages instanceof String) {
                    spec.entries = [packages]
                } else if (packages instanceof List) {
                    spec.entries = packages as List<String>
                }
            }
            
            if (map.containsKey('environment')) {
                spec.environment = map.environment as String
            }
            
            if (map.containsKey('channels')) {
                def channels = map.channels
                if (channels instanceof List) {
                    spec.channels = channels as List<String>
                } else if (channels instanceof String) {
                    spec.channels = [channels]
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
     * Check if the package manager feature is enabled
     * 
     * @param session The current session
     * @return True if the feature is enabled
     */
    static boolean isEnabled(ISession session) {
        return session.config.navigate('nextflow.preview.package', false) as Boolean
    }
}