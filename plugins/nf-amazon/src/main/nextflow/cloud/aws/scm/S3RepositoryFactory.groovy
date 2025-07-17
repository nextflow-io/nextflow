/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.cloud.aws.scm

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.scm.jgit.TransportS3
import nextflow.plugin.Priority
import nextflow.scm.GitUrl
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryFactory
import nextflow.scm.RepositoryProvider

import java.util.concurrent.atomic.AtomicBoolean

/**
 * Implements a factory to create an instance of {@link S3RepositoryProvider}
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@Priority(-10)
@CompileStatic
class S3RepositoryFactory extends RepositoryFactory{

    private static AtomicBoolean registered = new AtomicBoolean(false)

    @Override
    protected RepositoryProvider createProviderInstance(ProviderConfig config, String project) {
        if (!registered.get()) {
            registered.set(true)
            TransportS3.register()
        }

        return config.platform == 's3'
            ? new S3RepositoryProvider(project, config)
            : null
    }

    @Override
    protected ProviderConfig getConfig(List<ProviderConfig> providers, GitUrl url) {
        // do not care about non AWS codecommit url
        if( url.protocol != 's3' )
            return null

        // S3 repository config depends on the bucket name stored as domain
        def config = providers.find( it -> it.domain == url.domain )
        if( config ) {
            log.debug "Git url=$url (1) -> config=$config"
            return config
        }

        // still nothing, create a new instance
        config = new S3ProviderConfig(url.domain)


        return config
    }

    @Override
    protected ProviderConfig createConfigInstance(String name, Map attrs) {
        final copy = new HashMap(attrs)
        return copy.platform == 's3'
            ? new S3ProviderConfig(name, copy)
            : null
    }




}
