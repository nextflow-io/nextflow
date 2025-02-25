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


import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.plugin.Priority
import nextflow.util.StringUtils
import org.pf4j.ExtensionPoint
/**
 * Implements a factory for Git {@link RepositoryProvider}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Priority(0)
class RepositoryFactory implements ExtensionPoint {

    /**
     * Create a new instance of the {@link RepositoryProvider} if a matching configuration
     * is found
     *
     * @param config
     *      The {@link ProviderConfig} object holding the Git configuration
     * @param project
     *      The Git project name e.g. `nextflow-io/hello`
     * @return
     *      An instance of a {@link RepositoryProvider} that matches the specified config
     *      or {@code null} if no match is found
     */
    protected RepositoryProvider createProviderInstance(ProviderConfig config, String project) {

        switch(config.platform) {
            case 'github':
                return new GithubRepositoryProvider(project, config)

            case 'bitbucket':
                return new BitbucketRepositoryProvider(project, config)

            case 'bitbucketserver':
                return new BitbucketServerRepositoryProvider(project, config)

            case 'gitlab':
                return new GitlabRepositoryProvider(project, config)

            case 'gitea':
                return new GiteaRepositoryProvider(project, config)

            case 'azurerepos':
                return new AzureRepositoryProvider(project, config)

            case 'file':
                // remove the 'local' prefix for the file provider
                def localName = project.tokenize('/').last()
                return new LocalRepositoryProvider(localName, config)

            default:
                return null
        }
    }

    protected ProviderConfig createConfigInstance(String name, Map attrs) {
        return new ProviderConfig(name, attrs)
    }

    protected ProviderConfig getConfig(List<ProviderConfig> providers, GitUrl url) {
        final result = providers.find(it -> it.domain == url.domain)
        log.debug "Git url=$url -> config=$result"
        return result
    }

    // --==  static definitions ==--
    private static boolean codeCommitLoaded
    private static List<RepositoryFactory> factories0

    private static List<RepositoryFactory> factories() {
        // simple cache
        if( factories0 )
            return factories0
        // scan for available plugins
        final result = Plugins.getPriorityExtensions(RepositoryFactory)
        log.debug "Found Git repository result: ${ result.collect(it->it.class.simpleName) }"
        return factories0=result
    }

    static RepositoryProvider newRepositoryProvider(ProviderConfig config, String project) {
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
            throw new AbortOperationException("Unable to find a Git repository provider matching platform: ${config.platform}")
        }
        // return matching provider
        return provider
    }

    static ProviderConfig newProviderConfig(String name, Map<String,Object> attrs) {
        // check if it's needed to load new plugins
        if( (name=='codecommit' || attrs.platform=='codecommit') && !codeCommitLoaded ) {
            Plugins.startIfMissing('nf-codecommit')
            codeCommitLoaded=true
            factories0=null
        }

        final config = factories().findResult( it -> it.createConfigInstance(name, attrs) )
        if( !config ) {
            throw new AbortOperationException("Unable to create a Git repository configuration matching name=${name} and attributes=${StringUtils.stripSecrets(attrs)}")
        }
        return config
    }

    static ProviderConfig getProviderConfig(List<ProviderConfig> providers, GitUrl url) {
        if( url.domain.startsWith('git-codecommit.') && url.domain.endsWith('.amazonaws.com') && !codeCommitLoaded ) {
            Plugins.startIfMissing('nf-codecommit')
            codeCommitLoaded=true
            factories0=null
        }
        
        final provider = factories().findResult( it -> it.getConfig(providers, url) )
        if( !provider ) {
            throw new AbortOperationException("Unable to find a Git provider config for ${url}")
        }
        // return matching provider
        return provider
    }

}
