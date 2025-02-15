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

package nextflow.cloud.aws.codecommit

import groovy.util.logging.Slf4j
import nextflow.plugin.Priority
import nextflow.scm.GitUrl
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryFactory
import nextflow.scm.RepositoryProvider
/**
 * Implements a factory to create an instance of {@link AwsCodeCommitRepositoryProvider}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Priority(-10)  // <-- lower is higher, this is needed to override default provider behavior
class AwsCodeCommitFactory extends RepositoryFactory {

    @Override
    protected RepositoryProvider createProviderInstance(ProviderConfig config, String project) {
        return config.platform=='codecommit'
                ? new AwsCodeCommitRepositoryProvider(project,config)
                : null
    }

    @Override
    protected ProviderConfig getConfig(List<ProviderConfig> providers, GitUrl url) {
        // do not care about non AWS codecommit url
        if( !url.domain.startsWith('git-codecommit.') || !url.domain.endsWith('.amazonaws.com') )
            return null

        // CodeCommit hostname vary depending the AWS region
        // try to find the config for the specified region
        def config = providers.find( it -> it.domain==url.domain )
        if( config ) {
            log.debug "Git url=$url (1) -> config=$config"
            return config
        }
        // fallback on the platform name
        config = providers.find( it -> it.platform=='codecommit' && !it.server )
        if( config ) {
            config.setServer("${url.protocol}://${url.domain}")
            log.debug "Git url=$url (2) -> config=$config"
            return config
        }
        // still nothing, create a new instance
        config = new AwsCodeCommitProviderConfig(url.domain)
        if( url.user ) {
            log.debug "Git url=$url (3) -> config=$config"
            config.setUser(url.user)
        }

        return config
    }

    @Override
    protected ProviderConfig createConfigInstance(String name, Map attrs) {
        final copy = new HashMap(attrs)
        if( name == 'codecommit' ) {
            copy.platform = 'codecommit'
        }

        return copy.platform == 'codecommit'
                ? new AwsCodeCommitProviderConfig(name, copy)
                : null
    }
}
