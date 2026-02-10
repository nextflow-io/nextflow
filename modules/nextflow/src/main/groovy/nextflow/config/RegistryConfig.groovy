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

package nextflow.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Configuration scope for module registry settings
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@ScopeName("registry")
@Description("""
    The `registry` scope provides configuration for the Nextflow module registry.
    This includes registry URL(s) and authentication settings.
""")
@CompileStatic
class RegistryConfig implements ConfigScope {

    final static public String DEFAULT_REGISTRY_URL = 'https://registry.nextflow.io'

    @ConfigOption
    @Description("Primary registry URL")
    final private String url

    @ConfigOption
    @Description("List of registry URLs to try in order")
    final private List<String> urls

    @ConfigOption
    @Description("Authentication configuration per registry (registry URL -> token)")
    final private Map<String, String> auth

    /* required by extension point -- do not remove */
    RegistryConfig() {
        this.url = DEFAULT_REGISTRY_URL
        this.urls = []
        this.auth = [:]
    }

    RegistryConfig(Map opts) {
        this.url = opts.url ? opts.url as String : DEFAULT_REGISTRY_URL
        this.urls = opts.urls ? opts.urls as List<String> : []
        this.auth = opts.auth ? opts.auth as Map<String, String> : [:]
    }

    /**
     * Get the primary registry URL
     *
     * @return The registry URL
     */
    String getUrl() {
        return url
    }

    /**
     * Get all registry URLs (primary + fallbacks)
     *
     * @return List of registry URLs
     */
    List<String> getAllUrls() {
        List<String> result = []
        if (urls && !urls.isEmpty()) {
            result.addAll(urls)
        } else if (url) {
            result.add(url)
        } else {
            result.add(DEFAULT_REGISTRY_URL)
        }
        return result
    }

    /**
     * Get authentication token for a registry
     *
     * @param registryUrl The registry URL
     * @return The authentication token, or null if not configured
     */
    String getAuthToken(String registryUrl) {
        return auth?.get(registryUrl)
    }

    /**
     * Get authentication token from environment variable or config
     *
     * @param registryUrl The registry URL
     * @return The authentication token, or null if not found
     */
    String getAuthTokenResolved(String registryUrl) {
        // First check config
        def token = getAuthToken(registryUrl)
        if (token) {
            // Resolve environment variable references like ${NXF_REGISTRY_TOKEN}
            if (token.startsWith('${') && token.endsWith('}')) {
                def envVar = token.substring(2, token.length() - 1)
                token = System.getenv(envVar)
            }
            return token
        }

        // Fallback to NXF_REGISTRY_TOKEN environment variable
        return System.getenv('NXF_REGISTRY_TOKEN')
    }

    /**
     * Check if authentication is configured for a registry
     *
     * @param registryUrl The registry URL
     * @return true if authentication is available
     */
    boolean hasAuth(String registryUrl) {
        return getAuthTokenResolved(registryUrl) != null
    }
}
