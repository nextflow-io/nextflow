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
import nextflow.SysEnv
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
    This includes the registry URL(s) and API key for authentication.
""")
@CompileStatic
class RegistryConfig implements ConfigScope {

    final static public String DEFAULT_REGISTRY_URL = 'https://registry.nextflow.io/api'

    @ConfigOption
    @Description("Registry URL or list of registry URLs in priority order (primary URL first)")
    final private Collection<String> url

    @ConfigOption
    @Description("API key for authenticating with the primary registry")
    final private String apiKey

    /* required by extension point -- do not remove */
    RegistryConfig() {
        this.url = [DEFAULT_REGISTRY_URL]
        this.apiKey = null
    }

    RegistryConfig(Map opts) {
        final urlObject = opts.url ?: [DEFAULT_REGISTRY_URL]
        if (urlObject instanceof Collection<String>)
            this.url = urlObject as Collection<String>
        else
            this.url = [urlObject.toString()]
        this.apiKey = opts.apiKey as String
    }

    /**
     * Get the primary (first) registry URL
     *
     * @return The primary registry URL
     */
    String getUrl() {
        return this.url ? url[0] as String : DEFAULT_REGISTRY_URL
    }

    /**
     * Get all registry URLs (primary first, fallbacks after)
     *
     * @return Collection of registry URLs
     */
    Collection<String> getAllUrls() {
        return this.url ?: [DEFAULT_REGISTRY_URL]
    }

    /**
     * Get the API key for the primary registry.
     * Authentication is only supported for the primary registry.
     *
     * @return The API key, or 'NXF_REGISTRY_TOKEN' env value if not configured
     */
    String getApiKey() {
        return apiKey ?: SysEnv.get('NXF_REGISTRY_TOKEN')
    }
}
