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
 *
 */

package nextflow.scm


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortScmOperationException
import nextflow.plugin.Plugins
import nextflow.util.StringUtils
/**
 * Nextflow-specific implementation of repository factory loader.
 * Uses Nextflow's plugin system to discover repository factories with priority-based loading
 * and dynamic plugin activation (e.g., for CodeCommit).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Singleton
@CompileStatic
class NextflowRepositoryFactoryLoader implements RepositoryFactoryLoader {

    private static List<RepositoryFactory> factories0
    private static boolean codeCommitLoaded


    static List<RepositoryFactory> factories() {
        // simple cache
        if( factories0 != null )
            return factories0
        // scan for available plugins
        final result = Plugins.getPriorityExtensions(RepositoryFactory)
        log.debug "Found Git repository factories: ${ result.collect(it->it.class.simpleName) }"
        return factories0=result
    }

    @Override
    RepositoryProvider createRepositoryProvider(ProviderConfig config, String project) {
        // check if it's needed to load new plugins
        if( (config.name=='codecommit' || config.platform=='codecommit') && !codeCommitLoaded ) {
            Plugins.startIfMissing('nf-codecommit')
            codeCommitLoaded=true
            factories0=null
        }

        // scan all installed Git repository factories and find the first
        // returning an provider instance for the specified parameters
        final provider = factories().findResult( it -> it.createProviderInstance(config, project) )
        if( !provider ) {
            throw new AbortScmOperationException("Unable to find a Git repository provider matching platform: ${config.platform}")
        }
        // return matching provider
        return provider
    }

    @Override
    ProviderConfig createProviderConfig(String name, Map<String,Object> attrs) {
        // check if it's needed to load new plugins
        if( (name=='codecommit' || attrs.platform=='codecommit') && !codeCommitLoaded ) {
            Plugins.startIfMissing('nf-codecommit')
            codeCommitLoaded=true
            factories0=null
        }

        final config = factories().findResult( it -> it.createConfigInstance(name, attrs) )
        if( !config ) {
            throw new AbortScmOperationException("Unable to create a Git repository configuration matching name=${name} and attributes=${StringUtils.stripSecrets(attrs)}")
        }
        return config
    }

    @Override
    ProviderConfig getProviderConfig(List<ProviderConfig> providers, GitUrl url) {
        if( url.domain.startsWith('git-codecommit.') && url.domain.endsWith('.amazonaws.com') && !codeCommitLoaded ) {
            Plugins.startIfMissing('nf-codecommit')
            codeCommitLoaded=true
            factories0=null
        }

        final provider = factories().findResult( it -> it.getConfig(providers, url) )
        if( !provider ) {
            throw new AbortScmOperationException("Unable to find a Git provider config for ${url}")
        }
        // return matching provider
        return provider
    }

}
