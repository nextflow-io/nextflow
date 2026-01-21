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
 *
 */

package nextflow.scm

/**
 * Interface for repository factory loaders.
 * Implementations provide different strategies for loading and creating
 * repository providers and configurations (e.g., with or without plugin support).
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
trait RepositoryFactoryLoader {

    /**
     * Create a new RepositoryProvider instance for the given configuration and project
     *
     * @param config The provider configuration
     * @param project The project name
     * @return A new RepositoryProvider instance
     */
    abstract RepositoryProvider createRepositoryProvider(ProviderConfig config, String project)

    /**
     * Create a new ProviderConfig instance from name and attributes
     *
     * @param name The provider name
     * @param attrs The provider attributes/configuration
     * @return A new ProviderConfig instance
     */
    abstract ProviderConfig createProviderConfig(String name, Map<String,Object> attrs)

    /**
     * Get provider configuration for a given Git URL
     *
     * @param providers List of available provider configurations
     * @param url The Git URL to match
     * @return The matching ProviderConfig
     */
    abstract ProviderConfig getProviderConfig(List<ProviderConfig> providers, GitUrl url)

    /**
     * Create list of provider configurations from a configuration map
     *
     * @param config Configuration map containing provider definitions
     * @return List of ProviderConfig instances
     */
    List<ProviderConfig> createFromMap(Map<String,?> config){
        final providers = (Map<String, ?>) config?.providers

        List<ProviderConfig> result = []
        if( providers ) {
            for( String name : providers.keySet() ) {
                def attrs = (Map) providers.get(name)
                result << createProviderConfig(name, attrs)
            }
        }

        ProviderConfig.addDefaults(result)
        return result
    }

}